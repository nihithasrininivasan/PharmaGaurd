from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

from pathlib import Path

from .gene_coordinates import infer_gene_from_coordinates, iter_pharmacogenes, normalize_chrom
from .parser import VcfHeaderInfo, VcfVariant, iter_vcf_variants, read_vcf_header


PHARMACOGENES: Set[str] = set(iter_pharmacogenes())


@dataclass(frozen=True)
class ExtractedVariant:
    gene: str
    rsid: Optional[str]
    star: Optional[str]
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: Optional[str]
    zygosity: str
    raw_info: Dict
    filter: Optional[str] = None


def extract_pharmacogenes(
    variants: Sequence[VcfVariant],
    *,
    genes: Optional[Iterable[str]] = None,
) -> Dict[str, List[ExtractedVariant]]:
    """
    Group variants by pharmacogene.

    Primary strategy:
    - Use INFO.GENE if present
    Fallback:
    - Infer gene from approximate genomic coordinates
    """
    allowed = set(genes) if genes is not None else PHARMACOGENES
    out: Dict[str, List[ExtractedVariant]] = {g: [] for g in sorted(allowed)}

    for v in variants:
        gene = (v.gene or "").strip()
        if not gene:
            gene = infer_gene_from_coordinates(normalize_chrom(v.chrom), v.pos) or ""
        if gene not in allowed:
            continue
        out[gene].append(
            ExtractedVariant(
                gene=gene,
                rsid=v.rsid,
                star=v.star,
                chrom=v.chrom,
                pos=v.pos,
                ref=v.ref,
                alt=v.alt,
                gt=v.gt,
                zygosity=v.zygosity or "Unknown",
                raw_info=dict(v.info) if isinstance(v.info, dict) else {"INFO": v.info},
                filter=v.filter,
            )
        )

    # Drop empty gene keys for cleaner downstream JSON
    return {g: vs for g, vs in out.items() if vs}


def extract_pharmacogenes_from_vcf_path(
    path: Union[str, Path],
    *,
    genes: Optional[Iterable[str]] = None,
    max_variants_kept: Optional[int] = None,
    max_scanned_variants: Optional[int] = None,
    encoding: str = "utf-8",
) -> Tuple[VcfHeaderInfo, Dict[str, List[ExtractedVariant]]]:
    """
    Streaming extractor for large VCF files.

    Reads the VCF line-by-line and only retains variants that map to the target genes.
    """
    allowed = set(genes) if genes is not None else PHARMACOGENES
    out: Dict[str, List[ExtractedVariant]] = {}

    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    open_kwargs = {"mode": "rt", "encoding": encoding, "errors": "replace"} if p.suffix == ".gz" else {"mode": "r", "encoding": encoding, "errors": "replace", "newline": ""}
    with opener(p, **open_kwargs) as f:
        header, data_lines = read_vcf_header(f)

        for v in iter_vcf_variants(data_lines, samples=header.samples, max_variants=max_scanned_variants):
            gene = (v.gene or "").strip()
            if not gene:
                gene = infer_gene_from_coordinates(normalize_chrom(v.chrom), v.pos) or ""
            if gene not in allowed:
                continue

            bucket = out.setdefault(gene, [])
            bucket.append(
                ExtractedVariant(
                    gene=gene,
                    rsid=v.rsid,
                    star=v.star,
                    chrom=v.chrom,
                    pos=v.pos,
                    ref=v.ref,
                    alt=v.alt,
                    gt=v.gt,
                    zygosity=v.zygosity or "Unknown",
                    raw_info=dict(v.info) if isinstance(v.info, dict) else {"INFO": v.info},
                    filter=v.filter,
                )
            )
            if max_variants_kept is not None and len(bucket) >= max_variants_kept:
                # Keep bounded memory per gene; still continue scanning other genes
                continue

    return header, out
