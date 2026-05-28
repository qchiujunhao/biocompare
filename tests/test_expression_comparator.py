from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.expression import ExpressionComparator
from biocompare.core.engine import ComparisonEngine


class ExpressionComparatorTests(TestCase):
    def test_identical_expression_matrices_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left_expression.tsv"
            right = Path(temp_dir) / "right_expression.tsv"
            content = "gene_id\ts1\ts2\nG1\t10.5\t5.0\nG2\t0.0\t20.2\nG3\t30.1\t15.0\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = ExpressionComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["gene_overlap"], 1.0)
        self.assertEqual(report.metrics["mean_sample_similarity"], 1.0)
        self.assertEqual(report.metrics["mean_sample_top_genes_jaccard"], 1.0)

    def test_expression_metrics_capture_overlap_and_magnitude_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text(
                "gene_id\ts1\ts2\ts3\nG1\t10.0\t5.0\t0.5\nG2\t2.0\t8.0\t0.2\nG3\t0.0\t1.0\t3.0\nG4\t30.0\t20.0\t10.0\n",
                encoding="utf-8",
            )
            right.write_text(
                "gene_id\ts1\ts2\ts3\nG1\t10.5\t5.2\t0.4\nG2\t1.8\t7.5\t0.3\nG3\t0.1\t1.2\t2.8\nG5\t40.0\t25.0\t12.0\n",
                encoding="utf-8",
            )

            report = ExpressionComparator().compare(str(left), str(right), top_n=2)

        self.assertAlmostEqual(report.metrics["gene_overlap"], 3 / 5)
        self.assertEqual(report.details["shared_samples"], ["s1", "s2", "s3"])
        self.assertLess(report.overall_concordance, 1.0)
        self.assertIn("mean_sample_distribution_similarity", report.metrics)

    def test_engine_auto_selects_expression_when_filename_has_hint(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left_tpm.tsv"
            right = Path(temp_dir) / "right_tpm.tsv"
            left.write_text("gene_id\ts1\ts2\nG1\t10.0\t5.0\nG2\t0.0\t20.0\n", encoding="utf-8")
            right.write_text("gene_id\ts1\ts2\nG1\t11.0\t5.0\nG2\t0.0\t19.0\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "ExpressionComparator")
        self.assertIn("mean_gene_profile_similarity", report.metrics)

