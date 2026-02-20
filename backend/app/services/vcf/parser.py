from __future__ import annotations

import datetime as _dt
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple, Union


@dataclass(frozen=True)
class VcfVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: Optional[float] = None
    filter: Optional[str] = None
    info: Mapping[str, Union[str, int, float, bool, Sequence[str]]] = field(
        default_factory=dict
    )
    format: Sequence[str] = field(default_factory=tuple)
    sample: Mapping[str, str] = field(default_factory=dict)

    rsid: Optional[str] = None
    gene: Optional[str] = None
    star: Optional[str] = None
    gt: Optional[str] = None
    zygosity: Optional[str] = None  # 'Hom-Ref' | 'Het' | 'Hom-Alt' | 'Unknown'


@dataclass
class VcfParseResult:
    patient_id: Optional[str]
    vcf_version: Optional[str]
    genome_build: Optional[str]
    header: List[str]
    samples: List[str]
    variants: List[VcfVariant]
    parsed_at_iso: str
    quality_metrics: Dict[str, Union[bool, int, float, str]]


class VcfParseError(ValueError):
    pass


@dataclass(frozen=True)
class VcfHeaderInfo:
    patient_id: Optional[str]
    vcf_version: Optional[str]
    genome_build: Optional[str]
    header_lines: List[str]
    samples: List[str]


def parse_vcf(
    content: Union[str, bytes, Iterable[str]],
    *,
    max_variants: Optional[int] = None,
) -> VcfParseResult:
    """
    Parse a VCF v4.x file and extract common pharmacogenomics annotations.

    Tolerant behavior:
    - Accepts missing INFO tags (GENE/STAR/RS) and infers rsid from ID column
    - Accepts multiple ALT alleles; emits one variant per ALT
    - If sample columns exist, extracts GT from FORMAT where possible
    """

    lines = _normalize_to_lines(content)
    header, data_lines = read_vcf_header(lines)

    variants: List[VcfVariant] = []
    for v in iter_vcf_variants(data_lines, samples=header.samples, max_variants=max_variants):
        variants.append(v)
        if max_variants is not None and len(variants) >= max_variants:
            break

    metrics: Dict[str, Union[bool, int, float, str]] = {
        "vcf_parsing_success": True,
        "variant_count": len(variants),
        "sample_count": len(header.samples),
        "has_column_header": True,
        "vcf_version": header.vcf_version or "",
        "genome_build": header.genome_build or "Unknown",
    }
    return VcfParseResult(
        patient_id=header.patient_id,
        vcf_version=header.vcf_version,
        genome_build=header.genome_build,
        header=header.header_lines,
        samples=header.samples,
        variants=variants,
        parsed_at_iso=_dt.datetime.now(tz=_dt.timezone.utc).isoformat(),
        quality_metrics=metrics,
    )


def read_vcf_header(lines: Iterable[str]) -> Tuple[VcfHeaderInfo, Iterator[str]]:
    """
    Consume an iterable of VCF lines until the first variant record.

    Returns (header_info, remaining_lines_iterator).
    """
    it = iter(lines)
    header_lines: List[str] = []
    samples: List[str] = []
    vcf_version: Optional[str] = None
    patient_id: Optional[str] = None
    column_header_seen = False

    buffered_first_record: Optional[str] = None
    for raw in it:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("##"):
            header_lines.append(line)
            if line.lower().startswith("##fileformat="):
                vcf_version = line.split("=", 1)[1].strip()
            continue
        if line.startswith("#"):
            header_lines.append(line)
            if line.startswith("#CHROM"):
                column_header_seen = True
                cols = line.lstrip("#").split("\t")
                if len(cols) >= 10:
                    samples = cols[9:]
                    patient_id = samples[0] if samples else None
            continue

        buffered_first_record = line
        break

    if not column_header_seen:
        raise VcfParseError("Invalid VCF: missing #CHROM header line.")

    genome_build = _detect_genome_build(header_lines)

    def remaining() -> Iterator[str]:
        if buffered_first_record is not None:
            yield buffered_first_record
        yield from it

    return (
        VcfHeaderInfo(
            patient_id=patient_id,
            vcf_version=vcf_version,
            genome_build=genome_build,
            header_lines=header_lines,
            samples=samples,
        ),
        remaining(),
    )


def iter_vcf_variants(
    lines: Iterable[str],
    *,
    samples: Sequence[str],
    max_variants: Optional[int] = None,
    record_predicate: Optional[Callable[[str], bool]] = None,
) -> Iterator[VcfVariant]:
    """
    Yield `VcfVariant` records from VCF data lines.

    `record_predicate`, when provided, is applied to the raw record line (after
    stripping CRLF) to cheaply skip records before parsing.
    """
    seen = 0
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        if record_predicate is not None and not record_predicate(line):
            continue
        for v in _parse_variant_line(line, samples=samples):
            yield v
            seen += 1
            if max_variants is not None and seen >= max_variants:
                return


def _normalize_to_lines(content: Union[str, bytes, Iterable[str], Path]) -> Iterator[str]:
    # Auto-detect gzip compressed files via Path
    if isinstance(content, Path):
        if content.suffix == ".gz":
            with gzip.open(content, "rt", encoding="utf-8", errors="replace") as f:
                yield from f
        else:
            with content.open("r", encoding="utf-8", errors="replace", newline="") as f:
                yield from f
        return
    if isinstance(content, bytes):
        # VCFs are typically UTF-8/ASCII; ignore odd bytes instead of failing hard.
        text = content.decode("utf-8", errors="replace")
        yield from text.splitlines(True)
        return
    if isinstance(content, str):
        yield from content.splitlines(True)
        return
    yield from content


def _parse_info_field(info: str) -> Dict[str, Union[str, int, float, bool, List[str]]]:
    out: Dict[str, Union[str, int, float, bool, List[str]]] = {}
    if info == "." or info == "":
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" not in item:
            out[item] = True
            continue
        k, v = item.split("=", 1)
        if "," in v:
            out[k] = v.split(",")
        else:
            # Best-effort scalar typing
            vv: Union[str, int, float, bool]
            if v.isdigit():
                vv = int(v)
            else:
                try:
                    vv = float(v)
                except ValueError:
                    vv = v
            out[k] = vv
    return out


def _pick_first_str(x: Union[str, int, float, bool, Sequence[str], None]) -> Optional[str]:
    if x is None:
        return None
    if isinstance(x, (int, float, bool)):
        return str(x)
    if isinstance(x, str):
        return x
    if isinstance(x, Sequence):
        return str(x[0]) if len(x) > 0 else None
    return None


def _parse_variant_line(line: str, *, samples: Sequence[str]) -> Iterator[VcfVariant]:
    cols = line.split("\t")
    if len(cols) < 8:
        raise VcfParseError(f"Invalid VCF record (expected 8+ columns): {line[:80]}")

    chrom, pos_s, vid, ref, alt_s, qual_s, flt, info_s = cols[:8]
    pos = int(pos_s)
    qual = None if qual_s in (".", "") else float(qual_s)
    flt_val = None if flt in (".", "") else flt
    info = _parse_info_field(info_s)

    format_keys: Tuple[str, ...] = ()
    sample_map: Dict[str, str] = {}
    gt: Optional[str] = None
    if len(cols) >= 9:
        format_keys = tuple(cols[8].split(":")) if cols[8] not in (".", "") else ()
    if len(cols) >= 10 and samples:
        # Use first sample only for now (hackathon spec: single patient).
        sample_fields = cols[9].split(":")
        sample_map = {k: (sample_fields[i] if i < len(sample_fields) else "") for i, k in enumerate(format_keys)}
        gt = sample_map.get("GT")

    alt_alleles = alt_s.split(",") if alt_s not in (".", "") else ["."]

    # PharmaGuard spec says INFO tags include GENE/STAR/RS. Be flexible with casing.
    gene = _pick_first_str(info.get("GENE") or info.get("Gene") or info.get("gene"))
    star = _pick_first_str(info.get("STAR") or info.get("Star") or info.get("star"))
    rsid = _pick_first_str(info.get("RS") or info.get("Rsid") or info.get("RSID") or info.get("rs"))
    if rsid is None and isinstance(vid, str) and vid.startswith("rs"):
        rsid = vid

    for alt in alt_alleles:
        yield VcfVariant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=flt_val,
            info=info,
            format=format_keys,
            sample=sample_map,
            rsid=rsid,
            gene=gene,
            star=star,
            gt=gt,
            zygosity=infer_zygosity(gt),
        )



def _detect_genome_build(header_lines: List[str]) -> Optional[str]:
    """Detect genome build from VCF header lines."""
    for line in header_lines:
        line_lower = line.lower()
        if "##reference=" in line_lower or "##assembly=" in line_lower:
            if "grch38" in line_lower or "hg38" in line_lower:
                return "GRCh38"
            if "grch37" in line_lower or "hg19" in line_lower:
                return "GRCh37"
    return None

def _normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name to chrX format."""
    if chrom.lower().startswith("chr"):
        return chrom
    # Handle MT/M
    if chrom.upper() == "MT":
        return "chrM"
    if chrom.upper() == "M":
        return "chrM"
    return f"chr{chrom}"

def _parse_gt_indices(gt: str) -> List[int]:
    """Parse GT string into list of allele indices. Returns empty list if invalid."""
    if not gt or gt in (".", "./.", ".|."):
        return []
    sep = "|" if "|" in gt else "/"
    try:
        return [int(a) for a in gt.split(sep) if a.strip().isdigit()]
    except ValueError:
        return []

def _infer_zygosity_for_allele(gt_indices: List[int], allele_idx: int) -> str:
    """
    Infer zygosity for a specific ALT allele (1-based index) given GT indices.
    
    gt_indices: [0, 1] (Het), [1, 1] (Hom-Alt 1), [1, 2] (Het 1/Het 2)
    allele_idx: 1 (for first ALT), 2 (for second), etc.
    """
    if not gt_indices:
        return "Unknown"
    
    count = gt_indices.count(allele_idx)
    if count == 0:
        return "Absent" # Should probably not yield this variant?
    
    total_ploidy = len(gt_indices)
    if count == total_ploidy:
        return "Hom-Alt"
    else:
        return "Het"

def infer_zygosity(gt: Optional[str]) -> str:
    """Legacy helper for single-allele inference (deprecated internal usage)."""
    if not gt: return "Unknown"
    indices = _parse_gt_indices(gt)
    if not indices: return "Unknown"
    unique = set(indices)
    if len(unique) == 1:
        return "Hom-Ref" if 0 in unique else "Hom-Alt"
    return "Het"


def _parse_variant_line(line: str, *, samples: Sequence[str]) -> Iterator[VcfVariant]:
    cols = line.split("\t")
    if len(cols) < 8:
        raise VcfParseError(f"Invalid VCF record (expected 8+ columns): {line[:80]}")

    chrom_raw, pos_s, vid, ref, alt_s, qual_s, flt, info_s = cols[:8]
    chrom = _normalize_chrom(chrom_raw)
    pos = int(pos_s)
    qual = None if qual_s in (".", "") else float(qual_s)
    flt_val = None if flt in (".", "") else flt
    info = _parse_info_field(info_s)

    format_keys: Tuple[str, ...] = ()
    sample_map: Dict[str, str] = {}
    gt: Optional[str] = None
    gt_indices: List[int] = []
    
    if len(cols) >= 9:
        format_keys = tuple(cols[8].split(":")) if cols[8] not in (".", "") else ()
    if len(cols) >= 10 and samples:
        # Use first sample only for now (hackathon spec: single patient).
        sample_fields = cols[9].split(":")
        sample_map = {k: (sample_fields[i] if i < len(sample_fields) else "") for i, k in enumerate(format_keys)}
        gt = sample_map.get("GT")
        if gt:
            gt_indices = _parse_gt_indices(gt)

    alt_alleles = alt_s.split(",") if alt_s not in (".", "") else ["."]

    # PharmaGuard spec says INFO tags include GENE/STAR/RS. Be flexible with casing.
    gene = _pick_first_str(info.get("GENE") or info.get("Gene") or info.get("gene"))
    star = _pick_first_str(info.get("STAR") or info.get("Star") or info.get("star"))
    rsid = _pick_first_str(info.get("RS") or info.get("Rsid") or info.get("RSID") or info.get("rs"))
    if rsid is None and isinstance(vid, str) and vid.startswith("rs"):
        rsid = vid

    # If only one ALT (or dot), simple case
    if len(alt_alleles) == 1:
        alt = alt_alleles[0]
        # Calculate zygosity based on whether *this specific allele* (index 1) is present
        # If alt is ".", assume Hom-Ref or use GT
        if alt == ".":
             zyg = "Hom-Ref" # Effectively
        else:
             # Standard VCF: ALT index is 1
             zyg = _infer_zygosity_for_allele(gt_indices, 1)
        
        # If zygosity is "Absent" (e.g. GT=0/0 but ALT=G), we typically skip 
        # listing it as a "detected variant" in PGx context unless it's explicitly wanted.
        # But parser usually yields checks. Adapter filters "Hom-Ref".
        
        # Correction: If parser yields "Absent", adapter checks zygosity.
        # Legacy behavior yielded everything.
        # Let's yield it, but with correct zygosity interpretation.
        
        # If GT is 0/0, infer_zygosity_for_allele(..., 1) returns "Absent".
        # We map "Absent" to "Hom-Ref" from the perspective of *this variant*? 
        # No, "Hom-Ref" means 0/0. "Absent" means this ALT is not in the sample.
        
        # For backward compatibility with adapter:
        # Adapter checks `mapped_zyg = _map_zygosity(v.zygosity)`
        # _map_zygosity maps "Hom-Ref", "Het", "Hom-Alt".
        # If we return "Absent", adapter handles it? Adapter maps "Hom-Ref" -> "HOM_REF".
        
        if zyg == "Absent":
            # If the ALT is not in GT, then this variant call is effectively REF for this specific ALT.
            # E.g. ALT=G, GT=0/0. The sample is REF.
            zyg = "Hom-Ref"
            
        yield VcfVariant(
            chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual, filter=flt_val,
            info=info, format=format_keys, sample=sample_map, rsid=rsid, gene=gene, star=star,
            gt=gt, zygosity=zyg
        )
    else:
        # Multi-allelic case: ALT=G,T
        for i, alt in enumerate(alt_alleles):
            alt_idx = i + 1
            zyg = _infer_zygosity_for_allele(gt_indices, alt_idx)
            
            # If Absent (e.g. GT=0/1 (REF/G), and we are looking at T (index 2)),
            # then T is not present. We should yield it as Hom-Ref (or distinct 'Absent'?)
            # If we act as if this was a biallelic record A->T, and GT is 0/0 (not T), then it is Hom-Ref.
            
            if zyg == "Absent":
                zyg = "Hom-Ref"

            yield VcfVariant(
                chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual, filter=flt_val,
                info=info, format=format_keys, sample=sample_map, rsid=rsid, gene=gene, star=star,
                gt=gt, zygosity=zyg
            )


