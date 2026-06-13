"""Semantic comparison of bioinformatics pipeline outputs."""

from pipeconcord._version import __version__
from pipeconcord.core.engine import ComparisonEngine
from pipeconcord.core.report import ConcordanceReport

__all__ = ["ComparisonEngine", "ConcordanceReport", "__version__"]

