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
    Itemised confidence components.

    Each component is a penalty-adjusted score in [0, 1].
    ``final`` is always ``min(all non-None components)`` and is computed
    lazily so that callers can inspect or override individual fields.
    """

    variant_quality: float = 1.0
    allele_coverage: float = 1.0
    phase_resolution: float = 1.0
    cnv_evaluation: float = 1.0
    diplotype_determinism: float = 1.0
    cpic_applicability: float = 1.0

    # Penalty log — human-readable audit trail
    penalties_applied: List[str] = field(default_factory=list)

    # --------------- derived ------------------

    def all_scores(self) -> List[float]:
        """Return all component scores (excludes metadata fields)."""
        return [
            self.variant_quality,
            self.allele_coverage,
            self.phase_resolution,
            self.cnv_evaluation,
            self.diplotype_determinism,
            self.cpic_applicability,
        ]

    @property
    def final(self) -> float:
        """
        Final confidence = min of all component scores.
        Never exceeds the weakest link.
        """
        return max(0.0, min(self.all_scores()))

    def to_dict(self) -> Dict[str, object]:
        return {
            "variant_quality": round(self.variant_quality, 4),
            "allele_coverage": round(self.allele_coverage, 4),
            "phase_resolution": round(self.phase_resolution, 4),
            "cnv_evaluation": round(self.cnv_evaluation, 4),
            "diplotype_determinism": round(self.diplotype_determinism, 4),
            "cpic_applicability": round(self.cpic_applicability, 4),
            "final": round(self.final, 4),
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

# Genes requiring CNV evaluation
CNV_REQUIRED_GENES = frozenset({"CYP2D6"})


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
        """Penalise phase_resolution for unphased compound heterozygotes."""
        if is_compound_het and not has_phasing:
            bd.phase_resolution = max(0.0, 1.0 - PENALTY_UNPHASED_COMPOUND_HET)
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
    ) -> None:
        """Penalise diplotype_determinism based on resolution quality."""
        deduction = 0.0

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

    def calculate_confidence(
        self,
        gene: str,
        phenotype: str,
        diplotype: str,
        variant_count: int,
        has_quality_data: bool = True,
    ) -> float:
        """
        Orchestrate full confidence calculation.
        """
        bd = ConfidenceBreakdown()
        
        # 1. Variant quality
        if not has_quality_data:
            bd.variant_quality = 0.85
            
        # 2. Allele coverage
        self.apply_allele_coverage_penalties(
            bd, gene, [], [], has_coverage_data=(variant_count > 0)
        )

        # 3. Phase: Compound het without phasing gets 0.90
        # If >1 allele and not wildtype, and no phasing data
        is_compound_het = "/" in diplotype and len(set(diplotype.split('/'))) > 1 and "*1" not in diplotype
        self.apply_phase_penalties(
            bd, 
            is_compound_het=is_compound_het, 
            has_phasing=False
        )

        # 4. CNV
        self.apply_cnv_penalties(bd, gene, cnv_evaluated=False)

        # 5. Diplotype
        self.apply_diplotype_penalties(bd, diplotype)

        # 6. CPIC
        is_indeterminate = phenotype.lower() in ["indeterminate", "unknown"]
        self.apply_cpic_penalties(
            bd, 
            has_cpic_rule=(not is_indeterminate),
            phenotype_is_indeterminate=is_indeterminate
        )

        return bd.final


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
                "phase_resolution": 0.10,
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
