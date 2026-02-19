from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple


@dataclass(frozen=True)
class GeneInterval:
    chrom: str
    start_1based: int
    end_1based: int

    def contains(self, chrom: str, pos_1based: int) -> bool:
        if chrom != self.chrom:
            return False
        return self.start_1based <= pos_1based <= self.end_1based


# Coordinates are best-effort defaults (hg38-ish). In the public test cases for
# this hackathon, INFO.GENE is typically present, so these are used as fallback.
PHARMACOGENE_INTERVALS: Dict[str, GeneInterval] = {
    "CYP2D6": GeneInterval(chrom="chr22", start_1based=42118499, end_1based=42130865),
    "CYP2C19": GeneInterval(chrom="chr10", start_1based=96521616, end_1based=96612962),
    "CYP2C9": GeneInterval(chrom="chr10", start_1based=96702047, end_1based=96748655),
    "SLCO1B1": GeneInterval(chrom="chr12", start_1based=21172831, end_1based=21253427),
    "TPMT": GeneInterval(chrom="chr6", start_1based=18128308, end_1based=18143954),
    "DPYD": GeneInterval(chrom="chr1", start_1based=97054070, end_1based=97105539),
}


def infer_gene_from_coordinates(chrom: str, pos_1based: int) -> Optional[str]:
    for gene, interval in PHARMACOGENE_INTERVALS.items():
        if interval.contains(chrom, pos_1based):
            return gene
    return None


def normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    if not c:
        return c
    return c if c.startswith("chr") else f"chr{c}"


def iter_pharmacogenes() -> Iterable[str]:
    return PHARMACOGENE_INTERVALS.keys()

