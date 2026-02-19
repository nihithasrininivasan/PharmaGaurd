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
    }
    return VcfParseResult(
        patient_id=header.patient_id,
        vcf_version=header.vcf_version,
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

    def remaining() -> Iterator[str]:
        if buffered_first_record is not None:
            yield buffered_first_record
        yield from it

    return (
        VcfHeaderInfo(
            patient_id=patient_id,
            vcf_version=vcf_version,
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


def infer_zygosity(gt: Optional[str]) -> str:
    """
    Infer zygosity from a VCF GT (genotype) string.

    Returns one of:
    - 'Hom-Ref'  : 0/0 or 0|0
    - 'Het'      : 0/1, 1/0, 0|1, 1|0 etc.
    - 'Hom-Alt'  : 1/1, 2/2, 1|1 etc.
    - 'Unknown'  : missing (.) or unparseable
    """
    if not gt or gt in (".", "./.", ".|."):
        return "Unknown"
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    # Strip phasing chars if any
    alleles = [a.strip() for a in alleles if a.strip() not in ("", ".")]
    if not alleles:
        return "Unknown"
    unique = set(alleles)
    if len(unique) == 1:
        return "Hom-Ref" if "0" in unique else "Hom-Alt"
    return "Het"

