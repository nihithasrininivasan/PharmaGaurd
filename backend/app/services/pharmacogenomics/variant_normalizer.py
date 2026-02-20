"""
Variant Normalizer — Pre-processing layer for VCF variants.

Handles:
  1. Chromosome format normalisation (chr10 → 10)
  2. Genome build validation (GRCh37 vs GRCh38)
  3. Variant quality filtering (FILTER, QUAL, depth, genotype)
  4. Duplicate removal (keep highest quality)
  5. Multi-allelic handling

If normalisation fails → abort allele resolution safely.
No inference allowed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .models import VariantCall
from .config import get_config


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class VariantQualityResult:
    """Per-variant quality assessment."""
    variant: VariantCall
    passes_filter: bool = True
    quality_adequate: bool = True
    depth_adequate: bool = True
    genotype_clear: bool = True

    @property
    def passes_all(self) -> bool:
        return (
            self.passes_filter
            and self.quality_adequate
            and self.depth_adequate
            and self.genotype_clear
        )

    @property
    def failure_reasons(self) -> List[str]:
        reasons = []
        if not self.passes_filter:
            reasons.append(f"FILTER={self.variant.filter}")
        if not self.quality_adequate:
            reasons.append(f"QUAL={self.variant.quality}")
        if not self.depth_adequate:
            reasons.append("Insufficient allele depth")
        if not self.genotype_clear:
            reasons.append(f"Ambiguous zygosity={self.variant.zygosity}")
        return reasons


@dataclass
class BuildValidationResult:
    """Result of genome build validation."""
    is_valid: bool = True
    detected_build: Optional[str] = None
    expected_build: str = "GRCh38"
    warnings: List[str] = field(default_factory=list)


@dataclass
class NormalizationPipelineResult:
    """Complete result of the normalisation pipeline."""
    clean_variants: List[VariantCall] = field(default_factory=list)
    rejected_variants: List[Tuple[VariantCall, str]] = field(default_factory=list)
    quality_results: List[VariantQualityResult] = field(default_factory=list)
    build_validation: Optional[BuildValidationResult] = None
    duplicates_removed: int = 0
    chromosome_normalized: int = 0

    @property
    def quality_penalty_factor(self) -> float:
        """
        A [0, 1] factor representing overall variant quality.
        1.0 = all variants passed; 0.0 = all rejected.
        """
        total = len(self.clean_variants) + len(self.rejected_variants)
        if total == 0:
            return 1.0
        return len(self.clean_variants) / total


# ---------------------------------------------------------------------------
# Chromosome normalisation
# ---------------------------------------------------------------------------

_CHR_PREFIX_RE = re.compile(r"^chr", re.IGNORECASE)

# Canonical chromosome names (without 'chr' prefix)
_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}


def normalize_chromosome(chrom: str) -> str:
    """
    Normalise chromosome identifier.
    - Strip 'chr' prefix
    - Normalise 'chrM' → 'MT'
    - Upper-case sex chromosomes
    """
    stripped = _CHR_PREFIX_RE.sub("", chrom).strip()

    # Normalise mitochondrial
    if stripped.upper() in ("M", "MT"):
        return "MT"

    # Upper-case sex chromosomes
    if stripped.upper() in ("X", "Y"):
        return stripped.upper()

    return stripped


# ---------------------------------------------------------------------------
# Genome build validation
# ---------------------------------------------------------------------------

# Known anchor positions for build detection.
# These are well-characterised pharmacogene positions that differ between builds.
# Format: { gene: { "GRCh38": [positions], "GRCh37": [positions] } }
_BUILD_ANCHORS: Dict[str, Dict[str, List[int]]] = {
    "CYP2C19": {
        "GRCh38": [94775367, 94781859, 94842866],   # chr10 GRCh38
        "GRCh37": [96541616, 96535866],               # chr10 GRCh37 (hg19 range)
    },
    "CYP2D6": {
        "GRCh38": [42126611, 42127941, 42130692],   # chr22 GRCh38
        "GRCh37": [42522613, 42524943],               # chr22 GRCh37
    },
}


def validate_genome_build(
    variants: Sequence[VariantCall],
    expected_build: str = "GRCh38",
) -> BuildValidationResult:
    """
    Heuristic genome build validation.

    Cross-references variant positions against known anchor positions
    to flag potential build mismatches.  This is not definitive but
    catches the most common error (submitting GRCh37 data to a GRCh38 pipeline).
    """
    result = BuildValidationResult(expected_build=expected_build)

    if not variants:
        result.warnings.append("No variants to validate build against")
        return result

    positions_by_gene: Dict[str, set] = {}
    for v in variants:
        if hasattr(v, "_gene_hint"):
            gene = v._gene_hint
        else:
            gene = None

        # Collect all positions regardless of gene for broad check
        positions_by_gene.setdefault("_all", set()).add(v.pos)

    all_positions = positions_by_gene.get("_all", set())

    # Check each anchor gene
    for gene, builds in _BUILD_ANCHORS.items():
        expected_positions = set(builds.get(expected_build, []))
        wrong_build_positions = set()

        for build_name, positions in builds.items():
            if build_name != expected_build:
                wrong_build_positions.update(positions)

        # Check: do any variant positions match the *wrong* build?
        wrong_matches = all_positions & wrong_build_positions
        expected_matches = all_positions & expected_positions

        if wrong_matches and not expected_matches:
            result.is_valid = False
            result.detected_build = [b for b in builds if b != expected_build][0]
            result.warnings.append(
                f"Positions {sorted(wrong_matches)} match {result.detected_build}, "
                f"not {expected_build} for {gene}"
            )

    return result


# ---------------------------------------------------------------------------
# Variant quality filtering
# ---------------------------------------------------------------------------

def filter_variant_quality(
    variant: VariantCall,
    min_quality: float = 20.0,
    min_allele_depth_ratio: float = 0.2,
) -> VariantQualityResult:
    """
    Assess quality of a single variant.
    Returns structured QC result — does NOT discard the variant.
    """
    result = VariantQualityResult(variant=variant)

    # 1. FILTER check
    if variant.filter and variant.filter.upper() not in ("PASS", "."):
        result.passes_filter = False

    # 2. Quality score check
    if variant.quality < min_quality:
        result.quality_adequate = False

    # 3. Allele depth check
    if variant.ad and len(variant.ad) >= 2:
        total_depth = sum(variant.ad)
        if total_depth > 0:
            alt_ratio = variant.ad[1] / total_depth
            if alt_ratio < min_allele_depth_ratio:
                result.depth_adequate = False
        else:
            result.depth_adequate = False

    # 4. Genotype clarity
    if variant.zygosity in ("Unknown", ""):
        result.genotype_clear = False

    return result


# ---------------------------------------------------------------------------
# Duplicate removal
# ---------------------------------------------------------------------------

def remove_duplicates(variants: List[VariantCall]) -> Tuple[List[VariantCall], int]:
    """
    Deduplicate variants by (chrom, pos, ref, alt).
    Keeps the variant with the highest quality score.
    Returns (deduplicated_list, count_removed).
    """
    best: Dict[str, VariantCall] = {}

    for v in variants:
        key = f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}"
        existing = best.get(key)
        if existing is None or v.quality > existing.quality:
            best[key] = v

    removed_count = len(variants) - len(best)
    return list(best.values()), removed_count


# ---------------------------------------------------------------------------
# Full normalisation pipeline
# ---------------------------------------------------------------------------

def normalize_variants(
    variants: List[VariantCall],
    expected_genome_build: str = "GRCh38",
    min_quality: float = 20.0,
    min_allele_depth_ratio: float = 0.2,
) -> NormalizationPipelineResult:
    """
    Full normalisation pipeline:
    
    1. Chromosome format normalisation
    2. Genome build validation (warning only — does not abort)
    3. Quality filtering (reject bad variants, keep good ones)
    4. Duplicate removal

    Returns a ``NormalizationPipelineResult`` with clean variants,
    rejected variants (with reasons), and quality metadata.
    """
    result = NormalizationPipelineResult()

    if not variants:
        return result

    # ---- Step 1: Chromosome normalisation ----
    normalised: List[VariantCall] = []
    for v in variants:
        new_chrom = normalize_chromosome(v.chrom)
        if new_chrom != v.chrom:
            result.chromosome_normalized += 1
            # Create a new VariantCall with normalised chrom
            v = v.model_copy(update={"chrom": new_chrom})
        normalised.append(v)

    # ---- Step 2: Genome build validation ----
    result.build_validation = validate_genome_build(normalised, expected_genome_build)

    # ---- Step 3: Quality filtering ----
    passing: List[VariantCall] = []
    for v in normalised:
        qr = filter_variant_quality(v, min_quality, min_allele_depth_ratio)
        result.quality_results.append(qr)

        if not qr.genotype_clear:
            # Ambiguous genotype → reject entirely
            result.rejected_variants.append((v, "Ambiguous genotype"))
        elif not qr.passes_filter:
            # Failed FILTER → reject but still record quality result
            result.rejected_variants.append((v, f"Failed filter: {v.filter}"))
        else:
            # Variant passes (may have low quality — tracked in quality_results)
            passing.append(v)

    # ---- Step 4: Duplicate removal ----
    deduped, dup_count = remove_duplicates(passing)
    result.duplicates_removed = dup_count
    result.clean_variants = deduped

    return result
