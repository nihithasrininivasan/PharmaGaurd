from __future__ import annotations

import datetime as _dt
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple, Union, Set

# Import gene coordinate logic
try:
    from .gene_coordinates import infer_gene_from_coordinates, normalize_chrom
except ImportError:
    from gene_coordinates import infer_gene_from_coordinates, normalize_chrom


# ----------------------------------------------------------------------
# Constants & Configuration
# ----------------------------------------------------------------------

TARGET_PHARMACOGENES: Set[str] = {
    "CYP2D6",
    "CYP2C19",
    "CYP2C9",
    "SLCO1B1",
    "TPMT",
    "DPYD",
}

# Minimum QUAL score to accept a variant (None = no filter)
DEFAULT_MIN_QUAL: Optional[float] = 20.0

# Minimum read depth (DP) to include a variant in coverage calculation
DEFAULT_MIN_DP: int = 10


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
    zygosity: Optional[str] = None   # 'Hom-Ref' | 'Het' | 'Hom-Alt' | 'Unknown'
    depth: Optional[int] = None       # DP from FORMAT or INFO
    quality_pass: bool = True         # Whether the variant passes quality thresholds


@dataclass
class VcfParseResult:
    patient_id: Optional[str]
    vcf_version: Optional[str]
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
    header_lines: List[str]
    samples: List[str]
    has_fileformat: bool
    has_chrom_header: bool


def parse_vcf(
    content: Union[str, bytes, Iterable[str]],
    *,
    max_variants: Optional[int] = None,
    min_qual: Optional[float] = DEFAULT_MIN_QUAL,
    min_dp: int = DEFAULT_MIN_DP,
    require_pass_filter: bool = False,
) -> VcfParseResult:
    """
    Parse a VCF v4.x file and extract pharmacogenomics variants for
    the 6 target genes: CYP2D6, CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD.

    Args:
        content:             File bytes, string, Path or line iterable.
        max_variants:        Optional cap on variants returned.
        min_qual:            Minimum QUAL score to accept (default 20.0). Use None to disable.
        min_dp:              Minimum read depth (DP) for coverage_check assessment.
        require_pass_filter: If True, only accept FILTER=PASS records.
    """
    lines_iter = _normalize_to_lines(content)

    try:
        header_info, data_lines = read_vcf_header(lines_iter)
    except VcfParseError:
        return _create_failed_result("Invalid VCF header or empty file")

    if not header_info.has_fileformat or not header_info.has_chrom_header:
        return _create_failed_result(
            "Missing ##fileformat or #CHROM header",
            header_lines=header_info.header_lines,
        )

    variants: List[VcfVariant] = []
    genes_found: Set[str] = set()
    qual_failed = 0
    filter_failed = 0
    depth_values: List[int] = []

    for v in iter_vcf_variants(data_lines, samples=header_info.samples, max_variants=max_variants):
        # Quality filter
        if min_qual is not None and v.qual is not None and v.qual < min_qual:
            qual_failed += 1
            # Still include the variant but flag it
            v = _replace_variant(v, quality_pass=False)

        # FILTER column check (if required)
        if require_pass_filter and v.filter is not None and v.filter.upper() not in ("PASS", ".", ""):
            filter_failed += 1
            v = _replace_variant(v, quality_pass=False)

        if v.gene:
            genes_found.add(v.gene)
        if v.depth is not None:
            depth_values.append(v.depth)
        variants.append(v)

        if max_variants is not None and len(variants) >= max_variants:
            break

    # ── Coverage heuristic ──────────────────────────────────────────────────
    # Considers gene breadth. If DP data is present, also requires adequate depth.
    count_genes = len(genes_found)
    mean_dp = (sum(depth_values) / len(depth_values)) if depth_values else None
    has_dp_data = mean_dp is not None

    # Depth only penalises when actual DP data exists in the file
    depth_ok = (not has_dp_data) or (mean_dp >= min_dp)

    if count_genes > 3 and depth_ok:
        coverage_status = "PASS"
    elif count_genes > 0:
        coverage_status = "LOW"
    else:
        coverage_status = "FAIL"

    passed_variants = [v for v in variants if v.quality_pass]

    metrics: Dict[str, Union[bool, int, float, str]] = {
        "vcf_parsing_success": True,
        "variant_count": len(variants),
        "quality_passed_variant_count": len(passed_variants),
        "qual_failed_count": qual_failed,
        "filter_failed_count": filter_failed,
        "sample_count": len(header_info.samples),
        "has_column_header": header_info.has_chrom_header,
        "vcf_version": header_info.vcf_version or "unknown",
        "genes_detected_count": count_genes,
        "coverage_check": coverage_status,
        "genes_detected_list": ",".join(sorted(genes_found)),
        "mean_read_depth": round(mean_dp, 1) if has_dp_data else 0.0,
    }

    return VcfParseResult(
        patient_id=header_info.patient_id,
        vcf_version=header_info.vcf_version,
        header=header_info.header_lines,
        samples=header_info.samples,
        variants=variants,
        parsed_at_iso=_dt.datetime.now(tz=_dt.timezone.utc).isoformat(),
        quality_metrics=metrics,
    )


def _replace_variant(v: VcfVariant, **kwargs) -> VcfVariant:
    """Return a new VcfVariant with some fields replaced (since it is frozen)."""
    d = {f: getattr(v, f) for f in v.__dataclass_fields__}
    d.update(kwargs)
    return VcfVariant(**d)


def _create_failed_result(reason: str, header_lines: List[str] = None) -> VcfParseResult:
    return VcfParseResult(
        patient_id=None,
        vcf_version=None,
        header=header_lines or [],
        samples=[],
        variants=[],
        parsed_at_iso=_dt.datetime.now(tz=_dt.timezone.utc).isoformat(),
        quality_metrics={
            "vcf_parsing_success": False,
            "error_reason": reason,
            "coverage_check": "FAIL",
            "variant_count": 0,
            "quality_passed_variant_count": 0,
            "qual_failed_count": 0,
            "filter_failed_count": 0,
            "sample_count": 0,
            "has_column_header": False,
            "vcf_version": "unknown",
            "genes_detected_count": 0,
            "genes_detected_list": "",
            "mean_read_depth": 0.0,
        },
    )


def read_vcf_header(lines: Iterator[str]) -> Tuple[VcfHeaderInfo, Iterator[str]]:
    """
    Consume VCF lines until the first variant record.
    Returns (header_info, remaining_lines_iterator).
    """
    header_lines: List[str] = []
    samples: List[str] = []
    vcf_version: Optional[str] = None
    patient_id: Optional[str] = None
    column_header_seen = False
    fileformat_seen = False
    buffered_first_record: Optional[str] = None

    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("##"):
            header_lines.append(line)
            if line.lower().startswith("##fileformat="):
                vcf_version = line.split("=", 1)[1].strip()
                fileformat_seen = True
            continue
        if line.startswith("#CHROM"):
            header_lines.append(line)
            column_header_seen = True
            cols = line.lstrip("#").split("\t")
            if len(cols) >= 10:
                samples = cols[9:]
                patient_id = samples[0] if samples else None
            continue
        if line.startswith("#"):
            header_lines.append(line)
            continue
        buffered_first_record = line
        break

    info = VcfHeaderInfo(
        patient_id=patient_id,
        vcf_version=vcf_version,
        header_lines=header_lines,
        samples=samples,
        has_fileformat=fileformat_seen,
        has_chrom_header=column_header_seen,
    )

    def remaining_iter() -> Iterator[str]:
        if buffered_first_record is not None:
            yield buffered_first_record
        yield from lines

    return info, remaining_iter()


def iter_vcf_variants(
    lines: Iterable[str],
    *,
    samples: Sequence[str],
    max_variants: Optional[int] = None,
) -> Iterator[VcfVariant]:
    """Yield VcfVariant records from VCF data lines (already filtered to target genes)."""
    seen = 0
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        for v in _parse_variant_line(line, samples=samples):
            yield v
            seen += 1
            if max_variants is not None and seen >= max_variants:
                return


def _normalize_to_lines(content: Union[str, bytes, Iterable[str], Path]) -> Iterator[str]:
    if isinstance(content, Path):
        if content.suffix == ".gz":
            with gzip.open(content, "rt", encoding="utf-8", errors="replace") as f:
                yield from f
        else:
            with content.open("r", encoding="utf-8", errors="replace", newline="") as f:
                yield from f
        return
    if isinstance(content, bytes):
        text = content.decode("utf-8", errors="replace")
        yield from text.splitlines(True)
        return
    if isinstance(content, str):
        yield from content.splitlines(True)
        return
    yield from content


def _parse_info_field(info: str) -> Dict[str, Union[str, int, float, bool, List[str]]]:
    out: Dict[str, Union[str, int, float, bool, List[str]]] = {}
    if info in (".", ""):
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
            if v.isdigit():
                out[k] = int(v)
            else:
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
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


def _extract_depth(info: Dict, sample_map: Dict[str, str]) -> Optional[int]:
    """
    Extract read depth (DP). Priority:
    1. FORMAT DP for the sample
    2. INFO DP
    """
    fmt_dp = sample_map.get("DP")
    if fmt_dp is not None:
        try:
            return int(fmt_dp)
        except ValueError:
            pass
    info_dp = info.get("DP")
    if info_dp is not None:
        try:
            return int(info_dp)
        except (ValueError, TypeError):
            pass
    return None


def _parse_variant_line(line: str, *, samples: Sequence[str]) -> Iterator[VcfVariant]:
    cols = line.split("\t")
    if len(cols) < 8:
        # Malformed line → skip gracefully
        return

    chrom, pos_s, vid, ref, alt_s, qual_s, flt, info_s = cols[:8]

    try:
        pos = int(pos_s)
    except ValueError:
        return  # bad position like 'BADPOS'

    qual = None if qual_s in (".", "") else None
    if qual_s not in (".", ""):
        try:
            qual = float(qual_s)
        except ValueError:
            pass

    flt_val = None if flt in (".", "") else flt
    info = _parse_info_field(info_s)

    format_keys: Tuple[str, ...] = ()
    sample_map: Dict[str, str] = {}
    gt: Optional[str] = None

    if len(cols) >= 9:
        format_keys = tuple(cols[8].split(":")) if cols[8] not in (".", "") else ()

    if len(cols) >= 10 and samples:
        sample_fields = cols[9].split(":")
        sample_map = {
            k: (sample_fields[i] if i < len(sample_fields) else "")
            for i, k in enumerate(format_keys)
        }
        gt = sample_map.get("GT")

    alt_alleles = alt_s.split(",") if alt_s not in (".", "") else ["."]

    # ── Gene extraction ─────────────────────────────────────────────────────
    # 1. Try INFO GENE tag (flexible casing)
    gene_raw = _pick_first_str(info.get("GENE") or info.get("Gene") or info.get("gene"))

    # 2. Fallback: infer from GRCh38 coordinates
    if not gene_raw:
        norm_chrom = normalize_chrom(chrom)
        gene_raw = infer_gene_from_coordinates(norm_chrom, pos)

    # 3. Filter: only emit variants for target genes
    if gene_raw not in TARGET_PHARMACOGENES:
        return

    # ── Extract star allele ─────────────────────────────────────────────────
    star = _pick_first_str(info.get("STAR") or info.get("Star") or info.get("star"))

    # ── Extract RSID ────────────────────────────────────────────────────────
    rsid = _pick_first_str(
        info.get("RS") or info.get("Rsid") or info.get("RSID") or info.get("rs")
    )
    if rsid is None and isinstance(vid, str) and vid.startswith("rs"):
        rsid = vid

    # ── Extract read depth ──────────────────────────────────────────────────
    depth = _extract_depth(info, sample_map)

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
            gene=gene_raw,
            star=star,
            gt=gt,
            zygosity=infer_zygosity(gt),
            depth=depth,
            quality_pass=True,  # Will be overridden in parse_vcf if needed
        )


def infer_zygosity(gt: Optional[str]) -> str:
    """
    Infer zygosity from a VCF GT string.

    Returns:
      'Hom-Ref'  : e.g. 0/0 or 0|0
      'Het'      : e.g. 0/1, 1/0, 0|1
      'Hom-Alt'  : e.g. 1/1, 2/2
      'Unknown'  : missing, ./. or unparseable
    """
    if not gt or gt in (".", "./.", ".|."):
        return "Unknown"
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    alleles = [a.strip() for a in alleles if a.strip() not in ("", ".")]
    if not alleles:
        return "Unknown"
    unique = set(alleles)
    if len(unique) == 1:
        return "Hom-Ref" if "0" in unique else "Hom-Alt"
    return "Het"
