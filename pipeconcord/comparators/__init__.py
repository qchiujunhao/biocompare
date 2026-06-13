from __future__ import annotations

from pipeconcord.comparators.bam_stats import BAMStatsComparator
from pipeconcord.comparators.bed import BEDComparator
from pipeconcord.comparators.counts import CountsComparator
from pipeconcord.comparators.deg import DEGComparator
from pipeconcord.comparators.expression import ExpressionComparator
from pipeconcord.comparators.fasta import FASTAComparator
from pipeconcord.comparators.table import TableComparator
from pipeconcord.comparators.vcf import VCFComparator
from pipeconcord.core.registry import ComparatorRegistry


def register_builtin_comparators(registry: type[ComparatorRegistry] = ComparatorRegistry) -> None:
    registry.register(DEGComparator)
    registry.register(ExpressionComparator)
    registry.register(CountsComparator)
    registry.register(BEDComparator)
    registry.register(FASTAComparator)
    registry.register(VCFComparator)
    registry.register(BAMStatsComparator)
    registry.register(TableComparator)


__all__ = ["BAMStatsComparator", "BEDComparator", "CountsComparator", "DEGComparator", "ExpressionComparator", "FASTAComparator", "TableComparator", "VCFComparator", "register_builtin_comparators"]
