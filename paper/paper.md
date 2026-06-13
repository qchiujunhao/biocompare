---
title: "PipeConcord: Semantic Comparison of Bioinformatics Pipeline Outputs"
tags:
  - Python
  - bioinformatics
  - workflow testing
  - reproducibility
  - concordance
authors:
  - name: "PipeConcord contributors"
    affiliation: 1
affiliations:
  - name: "Open source contributors"
    index: 1
date: 13 June 2026
bibliography: paper.bib
---

# Summary

PipeConcord is a Python package and command-line tool for semantic comparison
of bioinformatics pipeline outputs. It produces standardized concordance
reports for common result types, including differential expression tables,
count matrices, normalized expression matrices, BED intervals, FASTA/FASTQ
records, VCF calls, `samtools` alignment summaries, and generic tabular files.

The package is designed for workflow regression testing and method comparison:
users can ask whether two outputs are scientifically similar even when they are
not byte-for-byte identical.

# Statement of Need

Bioinformatics workflows often produce outputs where exact file identity is too
strict and format validation is too weak. A pipeline change may reorder rows,
alter insignificant floating-point values, or shift interval boundaries while
preserving the scientific conclusion. Conversely, small-looking file changes
can flip differential expression direction, change variant genotypes, or alter
peak overlap.

Existing tools are strong within specific domains, such as interval operations
with BEDTools [@quinlan2010bedtools], alignment summaries with SAMtools
[@li2009sequence], or variant representation standards [@danecek2011variant].
Workflow testing tools can verify snapshots or checksums, but checksums do not
express scientific concordance. `pipeconcord` fills the gap between bitwise
testing and format-specific benchmarking by providing a unified comparator
interface and report schema across several common output types.

# State of the Field

Specialized benchmarking remains essential for many tasks. For example, variant
calling against truth sets requires dedicated tools and carefully normalized
inputs. `pipeconcord` is not intended to replace those benchmarks. Instead, it
targets routine pipeline development, where users need lightweight checks across
many heterogeneous outputs.

The closest general-purpose alternatives are ad hoc scripts, snapshot testing,
and format validators. These approaches answer different questions: whether a
file changed, whether a file is structurally valid, or whether one metric has a
specific value. `pipeconcord` instead returns a multi-metric concordance report
with a bounded summary score suitable for continuous integration thresholds.

# Software Design

The core API has three concepts:

- `Comparator`: a plugin interface with `can_handle()` and `compare()`.
- `ConcordanceReport`: a common result model with summary score, metrics,
  details, and warnings.
- `ComparisonEngine`: the dispatcher that selects a comparator and runs it.

The built-in comparators are dependency-free and operate on text formats. This
keeps installation simple in CI environments and encourages optional heavy
parsers to live behind explicit plugin boundaries. Reports can be emitted as
JSON, TSV, text, or self-contained HTML. Batch mode reads a manifest of file
pairs and exits non-zero when outputs fail a user-defined concordance threshold.

# Research Impact Statement

`pipeconcord` supports reproducible computational biology by making regression
checks more scientifically meaningful. It can help workflow authors detect
behavioral drift during refactoring, software upgrades, container rebuilds, and
parameter changes. The tool is most useful as a complement to domain-specific
benchmarks, not as a substitute for experimental validation or truth-set
evaluation.

# AI Usage Disclosure

Early code and documentation drafts were developed with assistance from an AI
coding agent. Human review is required before scientific release, including
verification of metric definitions, references, documentation, and examples.

# References
