from __future__ import annotations

from pipeconcord.comparators.base import Comparator
from pipeconcord.comparators.counts import parse_count_matrix, sample_columns_kwarg
from pipeconcord.comparators.table import read_table
from pipeconcord.core.report import ConcordanceReport
from pipeconcord.core.utils import clamp01, jaccard, mean_absolute_error, numeric_similarity, pearson, spearman
from pipeconcord.detection.filetype import detect_file_type


EXPRESSION_HINTS = {"tpm", "fpkm", "cpm", "expression", "expr", "normalized", "normalised"}


class ExpressionComparator(Comparator):
    """Comparator for normalized expression matrices such as TPM, FPKM, or CPM."""

    name = "expression"
    supported_types = ("expression", "normalized-expression")

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type in self.supported_types:
            return True
        if requested_type is not None:
            return False
        type_a = detect_file_type(file_a)
        type_b = detect_file_type(file_b)
        if type_a.kind not in {"table", "csv", "tsv"} or type_b.kind not in {"table", "csv", "tsv"}:
            return False
        if not (has_expression_hint(file_a) or has_expression_hint(file_b)):
            return False
        try:
            parse_count_matrix(read_table(file_a))
            parse_count_matrix(read_table(file_b))
        except ValueError:
            return False
        return True

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        delimiter = kwargs.get("delimiter")
        if delimiter is not None and not isinstance(delimiter, str):
            raise TypeError("delimiter must be a string")
        gene_column = string_kwarg(kwargs, "gene_column")
        sample_columns = sample_columns_kwarg(kwargs.get("sample_columns"))
        top_n = optional_int(kwargs.get("top_n"), 50, "top_n")

        matrix_a = parse_count_matrix(read_table(file_a, delimiter=delimiter), gene_column=gene_column, sample_columns=sample_columns)
        matrix_b = parse_count_matrix(read_table(file_b, delimiter=delimiter), gene_column=gene_column, sample_columns=sample_columns)

        genes_a = set(matrix_a.rows_by_gene)
        genes_b = set(matrix_b.rows_by_gene)
        shared_genes = sorted(genes_a & genes_b)
        samples_a = set(matrix_a.sample_columns)
        samples_b = set(matrix_b.sample_columns)
        shared_samples = [sample for sample in matrix_a.sample_columns if sample in samples_b]

        warnings: list[str] = []
        if matrix_a.duplicate_genes:
            warnings.append(f"file_a contains {matrix_a.duplicate_genes} duplicate gene identifiers; kept the first occurrence.")
        if matrix_b.duplicate_genes:
            warnings.append(f"file_b contains {matrix_b.duplicate_genes} duplicate gene identifiers; kept the first occurrence.")
        if not shared_genes:
            warnings.append("No shared gene identifiers were found.")
        if not shared_samples:
            warnings.append("No shared sample columns were found.")

        sample_pearsons: list[float] = []
        sample_spearmans: list[float] = []
        sample_similarities: list[float] = []
        sample_maes: list[float] = []
        sample_distribution_similarities: list[float] = []
        top_gene_jaccards: list[float] = []
        metrics: dict[str, float] = {
            "gene_overlap": jaccard(genes_a, genes_b),
            "sample_overlap": jaccard(samples_a, samples_b),
        }

        for sample in shared_samples:
            values_a = [float(matrix_a.rows_by_gene[gene][sample]) for gene in shared_genes]
            values_b = [float(matrix_b.rows_by_gene[gene][sample]) for gene in shared_genes]
            if not values_a:
                continue
            sample_pearson = pearson(values_a, values_b)
            sample_spearman = spearman(values_a, values_b)
            sample_similarity = numeric_similarity(values_a, values_b)
            sample_mae = mean_absolute_error(values_a, values_b)
            distribution_similarity = expression_distribution_similarity(values_a, values_b)
            top_jaccard = jaccard(top_genes_for_sample(matrix_a, sample, top_n), top_genes_for_sample(matrix_b, sample, top_n))
            metrics[f"sample.{sample}.pearson"] = sample_pearson
            metrics[f"sample.{sample}.spearman"] = sample_spearman
            metrics[f"sample.{sample}.similarity"] = sample_similarity
            metrics[f"sample.{sample}.mae"] = sample_mae
            metrics[f"sample.{sample}.distribution_similarity"] = distribution_similarity
            metrics[f"sample.{sample}.top_genes_jaccard"] = top_jaccard
            sample_pearsons.append(sample_pearson)
            sample_spearmans.append(sample_spearman)
            sample_similarities.append(sample_similarity)
            sample_maes.append(sample_mae)
            sample_distribution_similarities.append(distribution_similarity)
            top_gene_jaccards.append(top_jaccard)

        gene_profile_pearsons: list[float] = []
        gene_profile_spearmans: list[float] = []
        gene_profile_similarities: list[float] = []
        if len(shared_samples) >= 2:
            for gene in shared_genes:
                values_a = [float(matrix_a.rows_by_gene[gene][sample]) for sample in shared_samples]
                values_b = [float(matrix_b.rows_by_gene[gene][sample]) for sample in shared_samples]
                gene_profile_pearsons.append(pearson(values_a, values_b))
                gene_profile_spearmans.append(spearman(values_a, values_b))
                gene_profile_similarities.append(numeric_similarity(values_a, values_b))

        metrics["mean_sample_pearson"] = mean(sample_pearsons)
        metrics["mean_sample_spearman"] = mean(sample_spearmans)
        metrics["mean_sample_similarity"] = mean(sample_similarities)
        metrics["mean_sample_mae"] = mean(sample_maes)
        metrics["mean_sample_distribution_similarity"] = mean(sample_distribution_similarities)
        metrics["mean_sample_top_genes_jaccard"] = mean(top_gene_jaccards)
        metrics["mean_gene_profile_pearson"] = mean(gene_profile_pearsons)
        metrics["mean_gene_profile_spearman"] = mean(gene_profile_spearmans)
        metrics["mean_gene_profile_similarity"] = mean(gene_profile_similarities)

        sample_rank_score = clamp01((metrics["mean_sample_spearman"] + 1.0) / 2.0)
        gene_rank_score = clamp01((metrics["mean_gene_profile_spearman"] + 1.0) / 2.0)
        overall = clamp01(
            0.15 * metrics["gene_overlap"]
            + 0.10 * metrics["sample_overlap"]
            + 0.25 * sample_rank_score
            + 0.20 * metrics["mean_sample_similarity"]
            + 0.15 * metrics["mean_sample_distribution_similarity"]
            + 0.10 * metrics["mean_sample_top_genes_jaccard"]
            + 0.05 * gene_rank_score
        )

        details = {
            "top_n": top_n,
            "file_a_rows": len(matrix_a.table.rows),
            "file_b_rows": len(matrix_b.table.rows),
            "file_a_gene_column": matrix_a.gene_column,
            "file_b_gene_column": matrix_b.gene_column,
            "file_a_sample_columns": matrix_a.sample_columns,
            "file_b_sample_columns": matrix_b.sample_columns,
            "shared_genes": len(shared_genes),
            "file_a_only_genes": len(genes_a - genes_b),
            "file_b_only_genes": len(genes_b - genes_a),
            "shared_samples": shared_samples,
            "file_a_only_samples": sorted(samples_a - samples_b),
            "file_b_only_samples": sorted(samples_b - samples_a),
            "compared_cells": len(shared_genes) * len(shared_samples),
        }
        return ConcordanceReport(
            comparator=self.__class__.__name__,
            file_a=str(file_a),
            file_b=str(file_b),
            overall_concordance=overall,
            metrics=metrics,
            details=details,
            warnings=warnings,
        )


def top_genes_for_sample(matrix, sample: str, top_n: int) -> set[str]:
    ranked = sorted(
        ((float(row[sample]), gene) for gene, row in matrix.rows_by_gene.items()),
        key=lambda item: (-item[0], item[1]),
    )
    return {gene for _, gene in ranked[:top_n]}


def expression_distribution_similarity(left: list[float], right: list[float]) -> float:
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0
    return (
        numeric_similarity([mean(left)], [mean(right)])
        + numeric_similarity([median(left)], [median(right)])
        + numeric_similarity([quantile(left, 0.25)], [quantile(right, 0.25)])
        + numeric_similarity([quantile(left, 0.75)], [quantile(right, 0.75)])
    ) / 4.0


def quantile(values: list[float], probability: float) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(values)
    index = probability * (len(sorted_values) - 1)
    lower = int(index)
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = index - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def median(values: list[float]) -> float:
    return quantile(values, 0.5)


def mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def optional_int(value: object, default: int, name: str) -> int:
    if value is None:
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be positive")
    return parsed


def string_kwarg(kwargs: dict[str, object], name: str) -> str | None:
    value = kwargs.get(name)
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    return value


def has_expression_hint(path: str) -> bool:
    lowered = path.lower()
    return any(hint in lowered for hint in EXPRESSION_HINTS)

