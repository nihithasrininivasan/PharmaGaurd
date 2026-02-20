"""
Confidence Scoring Framework — Deterministic, Honest, Never Boosted.

Implements a mathematically defensible confidence scoring system for
pharmacogenomic diplotype calling and risk assessment. Follows the
principle: "Final confidence = min(component confidences)."

Rules:
  - Never boost confidence beyond data support
  - Never inflate score beyond weakest component
  - Fail safely when uncertain
  - Clamp all scores to [0.0, 1.0]
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional, Sequence

from .config import get_config


# ---------------------------------------------------------------------------
# Confidence Breakdown — every component tracked independently
# ---------------------------------------------------------------------------

@dataclass
class ConfidenceBreakdown:
    """
    Confidence Architecture — two-axis model.

    Separates:
      - Phenotype determination confidence (can we determine the phenotype?)
      - Classification confidence (how confident are we in the LABEL we assigned?)

    Genotype Confidence (weighted formula):
      0.35 * allele_coverage + 0.25 * cnv_evaluation +
      0.25 * variant_quality + 0.15 * genome_build_validity

    Phenotype Confidence:
      genotype_confidence × diplotype_determinism

    Classification Confidence (prevents zero-collapse):
      If phenotype resolved:   0.6 * phenotype_confidence + 0.4 * knowledge_confidence
      If phenotype unresolved: 0.6 * (1 - phenotype_confidence) + 0.4 * knowledge_confidence
      → "We are confident the result IS inconclusive" ≠ "we know nothing"

    Automation requires: resolved phenotype + strong evidence + confirmed pair.
    """

    # --- Genotype components (patient-specific VCF quality) ---
    variant_quality: float = 1.0
    allele_coverage: float = 1.0
    cnv_evaluation: float = 1.0
    genome_build_validity: float = 1.0

    # --- Phenotype component ---
    diplotype_determinism: float = 1.0

    # --- Guideline match ---
    cpic_applicability: float = 1.0

    # --- Knowledge (external evidence, NOT patient-specific) ---
    knowledge_confidence: float = 1.0

    # --- Audit fields ---
    penalties_applied: List[str] = field(default_factory=list)
    automation_blocked_reasons: List[str] = field(default_factory=list)

    # --- gene-drug confirmation flag ---
    gene_drug_confirmed: bool = True

    # --------------- derived axes ------------------

    @property
    def genotype_confidence(self) -> float:
        """
        Weighted genotype quality score.

        Formula:
          0.35 * allele_coverage +
          0.25 * cnv_evaluation +
          0.25 * variant_quality +
          0.15 * genome_build_validity

        Unresolved diplotype → components are capped by risk engine,
        ensuring genotype_confidence < 1.0.
        """
        raw = (
            0.35 * self.allele_coverage +
            0.25 * self.cnv_evaluation +
            0.25 * self.variant_quality +
            0.15 * self.genome_build_validity
        )
        return max(0.0, min(1.0, round(raw, 4)))

    @property
    def phenotype_confidence(self) -> float:
        """
        Phenotype determination confidence.

        phenotype_confidence = genotype_confidence × diplotype_determinism

        When diplotype is unresolved (diplotype_determinism = 0),
        phenotype_confidence = 0 regardless of genotype quality.
        When resolved, phenotype confidence is bounded by genotype quality.
        """
        return max(0.0, round(
            self.genotype_confidence * self.diplotype_determinism, 4
        ))

    @property
    def classification_confidence(self) -> float:
        """
        Confidence in the assigned risk LABEL (not in phenotype determination).

        Prevents zero-collapse: when phenotype is unresolved, the label
        "Inconclusive" can still be HIGH-confidence (we're confident
        we can't determine the phenotype).

        Formula:
          If resolved:   0.6 * phenotype_confidence + 0.4 * knowledge_confidence
          If unresolved: 0.6 * (1 − phenotype_confidence) + 0.4 * knowledge_confidence
        """
        pheno = self.phenotype_confidence
        knowledge = max(0.0, min(1.0, self.knowledge_confidence))

        if pheno > 0.0:
            # Resolved: confidence tracks phenotype + knowledge
            raw = 0.6 * pheno + 0.4 * knowledge
        else:
            # Unresolved: we're confident it IS inconclusive
            raw = 0.6 * (1.0 - pheno) + 0.4 * knowledge

        return max(0.0, min(1.0, round(raw, 4)))

    @property
    def final(self) -> float:
        """
        Final confidence score with deterministic caps.

        This is the value used as confidence_score on RiskAssessment.
        It represents confidence in the LABEL, not in the phenotype.

        Caps applied:
        - phenotype_confidence = 0 → cap at 0.50 (cannot be fully confident when phenotype unresolved)
        - automation_status.allowed = false → cap at 0.70 (blocked automation reduces confidence)

        These caps enforce mathematical consistency:
        - Overall confidence cannot be maximal (1.0) when phenotype_confidence = 0
        - Overall confidence cannot be maximal when automation gates fail
        """
        base_confidence = self.classification_confidence

        # Apply phenotype cap: if phenotype is unresolved, confidence cannot exceed 0.50
        phenotype_cap = 0.50 if self.phenotype_confidence == 0 else 1.0

        # Apply automation cap: if automation is blocked, confidence cannot exceed 0.70
        automation_status = self.get_automation_status()
        automation_cap = 0.70 if not automation_status.get("allowed", True) else 1.0

        # Final confidence = min(base, phenotype_cap, automation_cap)
        final_confidence = min(base_confidence, phenotype_cap, automation_cap)

        return max(0.0, min(1.0, round(final_confidence, 4)))

    def get_automation_status(self) -> Dict[str, object]:
        """
        4-gate automation check. ALL must pass:
          1. phenotype resolved (phenotype_confidence > 0)
          2. evidence sufficient (knowledge_confidence >= 0.80)
          3. genotype quality adequate (genotype_confidence >= 0.50)
          4. gene-drug pair confirmed

        Returns dict with 'allowed' and 'blocked_reasons'.
        """
        reasons = []

        if self.phenotype_confidence <= 0.0:
            reasons.append(
                f"Phenotype unresolved (phenotype_confidence = {self.phenotype_confidence:.2f})"
            )
        if self.knowledge_confidence < 0.80:
            reasons.append(
                f"Evidence insufficient (knowledge_confidence = {self.knowledge_confidence:.2f}, requires ≥ 0.80)"
            )
        if self.genotype_confidence < 0.50:
            reasons.append(
                f"Genotype quality too low (genotype_confidence = {self.genotype_confidence:.2f}, requires ≥ 0.50)"
            )
        if not self.gene_drug_confirmed:
            reasons.append(
                "Gene-drug pair not confirmed in PharmGKB"
            )

        return {
            "allowed": len(reasons) == 0,
            "blocked_reasons": reasons,
        }

    def all_scores(self) -> List[float]:
        """Return all component scores (backward compat)."""
        return [
            self.variant_quality,
            self.allele_coverage,
            self.cnv_evaluation,
            self.genome_build_validity,
            self.diplotype_determinism,
            self.cpic_applicability,
        ]

    def to_dict(self) -> Dict[str, object]:
        """
        Serialise confidence breakdown.

        NOTE: automation_status is NOT included here — it lives
        exclusively on RiskAssessment.automation_status to prevent
        duplication.
        """
        return {
            # Knowledge (external evidence)
            "knowledge_confidence": round(self.knowledge_confidence, 4),
            # Genotype (patient-specific VCF quality)
            "genotype_confidence": round(self.genotype_confidence, 4),
            "genotype_components": {
                "variant_quality": round(self.variant_quality, 4),
                "allele_coverage": round(self.allele_coverage, 4),
                "cnv_evaluation": round(self.cnv_evaluation, 4),
                "genome_build_validity": round(self.genome_build_validity, 4),
            },
            # Phenotype determination
            "phenotype_confidence": round(self.phenotype_confidence, 4),
            # Classification (label correctness)
            "classification_confidence": round(self.classification_confidence, 4),
            # Audit trail
            "penalties_applied": list(self.penalties_applied),
        }


# ---------------------------------------------------------------------------
# Penalty constants (additive deductions from 1.0 base)
# ---------------------------------------------------------------------------

# These are *deductions*, not multipliers.
# Variant quality
PENALTY_FAILED_FILTER = 0.30          # FILTER ≠ PASS
PENALTY_LOW_QUALITY = 0.15            # QUAL < threshold
PENALTY_LOW_DEPTH = 0.10              # AD ratio below threshold
PENALTY_AMBIGUOUS_GENOTYPE = 0.20     # GT = ./. or Unknown

# Allele coverage
PENALTY_MISSING_KEY_POSITION = 0.05   # Per missing key SNP
PENALTY_NO_COVERAGE_DATA = 0.15       # No coverage information provided
PENALTY_INCOMPLETE_ALLELE_DEF = 0.20  # Not all defining SNPs evaluated

# Phase & CNV
PENALTY_UNPHASED_COMPOUND_HET = 0.10  # Compound het without phasing
PENALTY_CNV_NOT_EVALUATED = 0.20      # CYP2D6 without CNV calling

# Diplotype determinism
PENALTY_INDETERMINATE = 0.50          # Diplotype is Indeterminate
PENALTY_PARTIAL_MATCH = 0.15          # Partial allele definition match
PENALTY_WILDTYPE_NO_COVERAGE = 0.30   # *1/*1 assumed but positions not verified

# CPIC applicability
PENALTY_NO_CPIC_RULE = 0.20           # No CPIC recommendation found
PENALTY_PHENOTYPE_INDETERMINATE = 0.30  # Phenotype could not be determined
PENALTY_UNKNOWN_GENOME_BUILD = 0.20   # Genome build unknown or missing
PENALTY_GRCH37_GENOME_BUILD = 0.15    # GRCh37 detected, GRCh38 expected

# Genes requiring CNV evaluation
CNV_REQUIRED_GENES = frozenset({"CYP2D6"})

# ---------------------------------------------------------------------------
# Diplotype Quality Categories
# ---------------------------------------------------------------------------
# Instead of hardcoded magic confidence numbers in _select_best_diplotype,
# we use named quality categories with documented, defensible penalties.

DIPLOTYPE_QUALITY_PENALTIES: Dict[str, float] = {
    "EXACT_HOM":               0.05,   # All defining SNPs match + HOM_ALT
    "EXACT_HET_WILDTYPE":      0.10,   # All SNPs match, paired with *1
    "COMPOUND_HET_PHASED":     0.10,   # Two alleles confirmed by phasing
    "COMPOUND_HET_UNPHASED":   0.20,   # Two alleles, phase unknown
    "VCF_ANNOTATION":          0.20,   # Used VCF STAR tag, not positional match
    "PARTIAL":                 0.35,   # Not all defining variants present
    "AMBIGUOUS":               0.50,   # >2 candidates without phasing
    "INDETERMINATE":           0.60,   # Cannot resolve diplotype
}

# Variant quality thresholds
VARIANT_QUAL_THRESHOLD = 30.0            # QUAL below this is low quality
VARIANT_DEPTH_THRESHOLD = 20             # Total depth below this is low
VARIANT_ALLELE_BALANCE_THRESHOLD = 0.15  # AD ratio below this is suspicious


# ---------------------------------------------------------------------------
# Calculator
# ---------------------------------------------------------------------------

class ConfidenceCalculator:
    """
    Deterministic confidence scoring.

    Usage::

        calc = ConfidenceCalculator()
        bd = ConfidenceBreakdown()

        calc.apply_variant_quality_penalties(bd, quality_results)
        calc.apply_allele_coverage_penalties(bd, gene, ...)
        ...

        final = bd.final   #  min(components)
    """

    # -- Variant Quality ---------------------------------------------------

    @staticmethod
    def apply_variant_quality_penalties(
        bd: ConfidenceBreakdown,
        quality_results: Sequence["VariantQualityResult"],
    ) -> None:
        """
        Penalise variant_quality component based on individual variant QC.

        Each failed criterion deducts from the base 1.0 score.
        Multiple low-quality variants compound the penalty.
        """
        if not quality_results:
            return

        total_variants = len(quality_results)
        if total_variants == 0:
            return

        # Count failures across all variants
        failed_filter_count = sum(1 for q in quality_results if not q.passes_filter)
        low_quality_count = sum(1 for q in quality_results if not q.quality_adequate)
        low_depth_count = sum(1 for q in quality_results if not q.depth_adequate)
        ambiguous_gt_count = sum(1 for q in quality_results if not q.genotype_clear)

        # Proportional penalties (fraction of variants affected)
        deduction = 0.0

        if failed_filter_count:
            frac = failed_filter_count / total_variants
            penalty = PENALTY_FAILED_FILTER * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"FILTER≠PASS on {failed_filter_count}/{total_variants} variants (−{penalty:.2f})"
            )

        if low_quality_count:
            frac = low_quality_count / total_variants
            penalty = PENALTY_LOW_QUALITY * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Low QUAL on {low_quality_count}/{total_variants} variants (−{penalty:.2f})"
            )

        if low_depth_count:
            frac = low_depth_count / total_variants
            penalty = PENALTY_LOW_DEPTH * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Low depth on {low_depth_count}/{total_variants} variants (−{penalty:.2f})"
            )

        if ambiguous_gt_count:
            frac = ambiguous_gt_count / total_variants
            penalty = PENALTY_AMBIGUOUS_GENOTYPE * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Ambiguous GT on {ambiguous_gt_count}/{total_variants} variants (−{penalty:.2f})"
            )

        bd.variant_quality = max(0.0, 1.0 - deduction)

    @staticmethod
    def apply_variant_quality_from_vcf(
        bd: ConfidenceBreakdown,
        variants: Sequence,
    ) -> None:
        """
        Assess variant quality directly from VCF VariantCall objects.

        Inspects each variant's FILTER, QUAL, and AD fields to produce
        honest variant_quality penalties based on actual sequencing data.

        This is the PRIMARY method that should be called in the main flow.
        ``apply_variant_quality_penalties`` is for pre-computed QC results.
        """
        if not variants:
            return

        total = len(variants)
        failed_filter = 0
        low_qual = 0
        low_depth = 0
        ambiguous_gt = 0

        for v in variants:
            # 1. FILTER check
            filt = getattr(v, 'filter', None)
            if filt and filt not in ('PASS', '.', None):
                failed_filter += 1

            # 2. QUAL check
            qual = getattr(v, 'quality', 0.0) or 0.0
            if qual < VARIANT_QUAL_THRESHOLD:
                low_qual += 1

            # 3. Depth check (from AD field)
            ad = getattr(v, 'ad', None)
            if ad and len(ad) >= 2:
                total_depth = sum(ad)
                if total_depth < VARIANT_DEPTH_THRESHOLD:
                    low_depth += 1
                else:
                    # Allele balance check
                    alt_ratio = ad[1] / total_depth if total_depth > 0 else 0
                    if alt_ratio < VARIANT_ALLELE_BALANCE_THRESHOLD:
                        low_depth += 1  # Suspicious allele balance

            # 4. Genotype check
            zyg = getattr(v, 'zygosity', '')
            if zyg in ('./.', 'UNKNOWN', '', None):
                ambiguous_gt += 1

        # Apply proportional penalties
        deduction = 0.0

        if failed_filter:
            frac = failed_filter / total
            penalty = PENALTY_FAILED_FILTER * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"FILTER≠PASS on {failed_filter}/{total} variants (−{penalty:.2f})"
            )

        if low_qual:
            frac = low_qual / total
            penalty = PENALTY_LOW_QUALITY * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Low QUAL (<{VARIANT_QUAL_THRESHOLD}) on {low_qual}/{total} variants (−{penalty:.2f})"
            )

        if low_depth:
            frac = low_depth / total
            penalty = PENALTY_LOW_DEPTH * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Low depth/allele-balance on {low_depth}/{total} variants (−{penalty:.2f})"
            )

        if ambiguous_gt:
            frac = ambiguous_gt / total
            penalty = PENALTY_AMBIGUOUS_GENOTYPE * frac
            deduction += penalty
            bd.penalties_applied.append(
                f"Ambiguous genotype on {ambiguous_gt}/{total} variants (−{penalty:.2f})"
            )

        if deduction > 0:
            bd.variant_quality = max(0.0, 1.0 - deduction)

    # -- Genome Build ------------------------------------------------------

    @staticmethod
    def apply_genome_build_penalty(
        bd: ConfidenceBreakdown,
        genome_build: Optional[str],
    ) -> None:
        """Penalise allele_coverage when genome build is unknown or mismatched."""
        if not genome_build or genome_build.lower() in ("unknown", ""):
            penalty = PENALTY_UNKNOWN_GENOME_BUILD
            bd.allele_coverage = min(bd.allele_coverage, max(0.0, 1.0 - penalty))
            bd.penalties_applied.append(
                f"Unknown genome build (−{penalty:.2f})"
            )
        elif genome_build in ("GRCh37", "hg19"):
            penalty = PENALTY_GRCH37_GENOME_BUILD
            bd.allele_coverage = min(bd.allele_coverage, max(0.0, 1.0 - penalty))
            bd.penalties_applied.append(
                f"GRCh37 build detected, GRCh38 expected (−{penalty:.2f})"
            )

    # -- Allele Coverage ---------------------------------------------------

    @staticmethod
    def apply_allele_coverage_penalties(
        bd: ConfidenceBreakdown,
        gene: str,
        observed_positions: Sequence[int],
        key_positions: Sequence[int],
        has_coverage_data: bool,
    ) -> None:
        """Penalise allele_coverage based on how many key positions were evaluated."""
        if not has_coverage_data:
            bd.allele_coverage = max(0.0, 1.0 - PENALTY_NO_COVERAGE_DATA)
            bd.penalties_applied.append(
                f"No coverage data provided for {gene} (−{PENALTY_NO_COVERAGE_DATA:.2f})"
            )
            return

        if not key_positions:
            return  # No key positions defined — cannot penalise

        key_set = set(key_positions)
        observed_set = set(observed_positions)
        missing = key_set - observed_set

        if missing:
            n = len(missing)
            penalty = min(PENALTY_MISSING_KEY_POSITION * n, 0.50)  # cap
            bd.allele_coverage = max(0.0, 1.0 - penalty)
            bd.penalties_applied.append(
                f"{n} key position(s) missing for {gene} (−{penalty:.2f})"
            )

    # -- Phase Resolution --------------------------------------------------

    @staticmethod
    def apply_phase_penalties(
        bd: ConfidenceBreakdown,
        is_compound_het: bool,
        has_phasing: bool,
    ) -> None:
        """Penalise allele_coverage for unphased compound heterozygotes."""
        if is_compound_het and not has_phasing:
            bd.allele_coverage = max(0.0, bd.allele_coverage - PENALTY_UNPHASED_COMPOUND_HET)
            bd.penalties_applied.append(
                f"Unphased compound heterozygote (−{PENALTY_UNPHASED_COMPOUND_HET:.2f})"
            )

    # -- CNV Evaluation ----------------------------------------------------

    @staticmethod
    def apply_cnv_penalties(
        bd: ConfidenceBreakdown,
        gene: str,
        cnv_evaluated: bool = False,
    ) -> None:
        """Penalise cnv_evaluation for genes where CNV is relevant but not assessed."""
        if gene in CNV_REQUIRED_GENES and not cnv_evaluated:
            bd.cnv_evaluation = max(0.0, 1.0 - PENALTY_CNV_NOT_EVALUATED)
            bd.penalties_applied.append(
                f"CNV not evaluated for {gene} (−{PENALTY_CNV_NOT_EVALUATED:.2f})"
            )

    # -- Diplotype Determinism ---------------------------------------------

    @staticmethod
    def apply_diplotype_penalties(
        bd: ConfidenceBreakdown,
        diplotype: str,
        is_partial_match: bool = False,
        is_wildtype_unverified: bool = False,
        quality_category: Optional[str] = None,
    ) -> None:
        """
        Penalise diplotype_determinism based on resolution quality.

        Args:
            quality_category: One of the DIPLOTYPE_QUALITY_PENALTIES keys.
                If provided, uses the standardized penalty instead of
                ad-hoc partial_match/wildtype flags.
        """
        deduction = 0.0

        # New path: use quality category if provided
        if quality_category and quality_category in DIPLOTYPE_QUALITY_PENALTIES:
            deduction = DIPLOTYPE_QUALITY_PENALTIES[quality_category]
            bd.penalties_applied.append(
                f"Diplotype quality: {quality_category} (−{deduction:.2f})"
            )
        else:
            # Legacy path: flags-based penalties
            if diplotype in ("Indeterminate", "Unknown"):
                deduction += PENALTY_INDETERMINATE
                bd.penalties_applied.append(
                    f"Indeterminate diplotype (−{PENALTY_INDETERMINATE:.2f})"
                )

            if is_partial_match:
                deduction += PENALTY_PARTIAL_MATCH
                bd.penalties_applied.append(
                    f"Partial allele definition match (−{PENALTY_PARTIAL_MATCH:.2f})"
                )

            if is_wildtype_unverified:
                deduction += PENALTY_WILDTYPE_NO_COVERAGE
                bd.penalties_applied.append(
                    f"Wildtype assumed without verifying key positions (−{PENALTY_WILDTYPE_NO_COVERAGE:.2f})"
                )

        if deduction > 0:
            bd.diplotype_determinism = max(0.0, 1.0 - deduction)

    # -- CPIC Applicability ------------------------------------------------

    @staticmethod
    def apply_cpic_penalties(
        bd: ConfidenceBreakdown,
        has_cpic_rule: bool,
        phenotype_is_indeterminate: bool,
    ) -> None:
        """Penalise cpic_applicability when CPIC guidance is incomplete."""
        deduction = 0.0

        if not has_cpic_rule:
            deduction += PENALTY_NO_CPIC_RULE
            bd.penalties_applied.append(
                f"No specific CPIC recommendation found (−{PENALTY_NO_CPIC_RULE:.2f})"
            )

        if phenotype_is_indeterminate:
            deduction += PENALTY_PHENOTYPE_INDETERMINATE
            bd.penalties_applied.append(
                f"Phenotype is Indeterminate (−{PENALTY_PHENOTYPE_INDETERMINATE:.2f})"
            )

        if deduction > 0:
            bd.cpic_applicability = max(0.0, 1.0 - deduction)


    # -- Weighted Confidence (Alternative to Min-Based) -------------------

    @staticmethod
    def calculate_weighted_confidence(
        bd: ConfidenceBreakdown,
        weights: Optional[Dict[str, float]] = None
    ) -> float:
        """
        Calculate confidence using weighted geometric mean.

        Alternative to min-based approach — less conservative.
        Allows tuning component importance via weights.

        Args:
            bd: ConfidenceBreakdown with component scores
            weights: Component weights (default: equal weighting)

        Returns:
            Weighted confidence score [0, 1]
        """
        import math

        # Default weights (can be tuned per gene)
        if weights is None:
            weights = {
                "variant_quality": 0.25,
                "allele_coverage": 0.20,
                "genome_build_validity": 0.10,
                "cnv_evaluation": 0.15,
                "diplotype_determinism": 0.25,
                "cpic_applicability": 0.05,
            }

        # Weighted geometric mean
        log_sum = 0.0
        for component, weight in weights.items():
            score = getattr(bd, component)
            # Prevent log(0) by clamping to minimum 0.01
            log_sum += weight * math.log(max(0.01, score))

        confidence = math.exp(log_sum)

        # Clamp to [0, 1]
        return max(0.0, min(1.0, confidence))


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

def create_confidence_calculator() -> ConfidenceCalculator:
    return ConfidenceCalculator()


def calculate_weighted_confidence(
    bd: ConfidenceBreakdown,
    weights: Optional[Dict[str, float]] = None
) -> float:
    """
    Convenience function for weighted confidence calculation.

    See ConfidenceCalculator.calculate_weighted_confidence for details.
    """
    return ConfidenceCalculator.calculate_weighted_confidence(bd, weights)
