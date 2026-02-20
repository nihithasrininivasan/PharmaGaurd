"""
Risk Scoring Module - Continuous risk score calculation (0-100 scale).

Implements weighted composite scoring with:
- Base severity mapping
- Phenotype-specific modifiers
- Confidence penalties
- Population rarity adjustments
- Feedback learning integration
"""

from typing import Dict, Optional
import math


# ============================================================================
# Base Severity Scores
# ============================================================================

SEVERITY_BASE_SCORES: Dict[str, float] = {
    "critical": 95.0,
    "high": 80.0,
    "moderate": 55.0,
    "low": 30.0,
    "none": 10.0,
}


# ============================================================================
# Phenotype Impact Modifiers
# ============================================================================

PHENOTYPE_MODIFIERS: Dict[str, float] = {
    # CYP metabolizer phenotypes
    "PM": +10.0,                    # Poor Metabolizer → highest risk
    "IM": +5.0,                     # Intermediate Metabolizer
    "RM": +3.0,                     # Rapid Metabolizer (toxicity risk)
    "UM": +8.0,                     # Ultra-rapid Metabolizer
    "NM": 0.0,                      # Normal Metabolizer

    # Alternative nomenclature
    "Poor Metabolizer": +10.0,
    "Intermediate Metabolizer": +5.0,
    "Rapid Metabolizer": +3.0,
    "Ultra-rapid Metabolizer": +8.0,
    "Normal Metabolizer": 0.0,
    "Extensive Metabolizer": 0.0,

    # TPMT activity phenotypes
    "Low Activity": +8.0,
    "Intermediate Activity": +5.0,
    "Normal Activity": 0.0,
    "High Activity": +2.0,

    # DPYD function phenotypes
    "Poor Function": +10.0,
    "Intermediate Function": +5.0,
    "Normal Function": 0.0,

    # Generic/fallback
    "Indeterminate": -5.0,          # Uncertainty → reduce risk score
    "Unknown": -5.0,
}


# ============================================================================
# Risk Score Calculator
# ============================================================================

class RiskScoreCalculator:
    """
    Calculate continuous risk scores (0-100) from multiple weighted factors.

    Components:
        1. Base severity score (primary driver)
        2. Phenotype-specific modifier (±10 points)
        3. Confidence penalty (low confidence → reduced score)
        4. Population rarity bonus (rare variant → slight increase)
        5. Feedback learning boost (historical corrections)

    Formula:
        risk_score = (base + phenotype_mod) × confidence_factor
                     + rarity_bonus × feedback_boost

    All scores bounded to [0, 100].
    """

    def __init__(
        self,
        severity_scores: Optional[Dict[str, float]] = None,
        phenotype_modifiers: Optional[Dict[str, float]] = None,
    ):
        self.severity_scores = severity_scores or SEVERITY_BASE_SCORES
        self.phenotype_modifiers = phenotype_modifiers or PHENOTYPE_MODIFIERS

    def calculate_risk_score(
        self,
        severity: str,
        phenotype: str,
        confidence: float,
        population_frequency: float = 0.5,
        feedback_boost: float = 1.0,
    ) -> float:
        """
        Calculate final risk score from all components.

        Args:
            severity: Categorical severity (critical/high/moderate/low/none)
            phenotype: Metabolizer/function phenotype
            confidence: Confidence score [0, 1]
            population_frequency: Diplotype frequency in population [0, 1]
            feedback_boost: Learning prior from clinician feedback [0.8, 1.5]

        Returns:
            Risk score in range [0, 100]
        """
        # 1. Base severity score
        base_score = self.severity_scores.get(severity, 55.0)

        # 2. Phenotype modifier
        phenotype_mod = self.phenotype_modifiers.get(phenotype, 0.0)

        # 3. Confidence penalty
        # Low confidence → reduce score to avoid false alarms
        # confidence=1.0 → factor=1.0 (no penalty)
        # confidence=0.5 → factor=0.85 (15% reduction)
        # confidence=0.0 → factor=0.70 (30% reduction)
        confidence_factor = 0.70 + (0.30 * confidence)

        # 4. Population rarity bonus
        # Rare variants → slight increase (novel risk potential)
        if population_frequency < 0.001:        # <0.1%
            rarity_bonus = 8.0
        elif population_frequency < 0.01:       # <1%
            rarity_bonus = 5.0
        elif population_frequency < 0.05:       # <5%
            rarity_bonus = 2.0
        else:
            rarity_bonus = 0.0

        # 5. Feedback learning integration
        # feedback_boost > 1.0 → increase score (clinicians flagged this)
        # feedback_boost < 1.0 → decrease score (overcalled in past)
        feedback_factor = feedback_boost

        # Calculate final score
        # Apply confidence penalty to base score, then add rarity bonus,
        # then scale by feedback
        risk_score = (
            (base_score + phenotype_mod) * confidence_factor + rarity_bonus
        ) * feedback_factor

        # Clamp to valid range [0, 100]
        return max(0.0, min(100.0, risk_score))

    def get_risk_level(self, risk_score: float) -> str:
        """
        Map continuous risk score to categorical risk level.

        Bands:
            90-100: critical
            70-89:  high
            40-69:  moderate
            20-39:  low
            0-19:   none
        """
        if risk_score >= 90:
            return "critical"
        elif risk_score >= 70:
            return "high"
        elif risk_score >= 40:
            return "moderate"
        elif risk_score >= 20:
            return "low"
        else:
            return "none"

    def calculate_with_explanation(
        self,
        severity: str,
        phenotype: str,
        confidence: float,
        population_frequency: float = 0.5,
        feedback_boost: float = 1.0,
    ) -> Dict:
        """
        Calculate risk score with detailed breakdown.

        Returns:
            {
                "risk_score": 78.5,
                "risk_level": "high",
                "components": {
                    "base_severity_score": 80.0,
                    "phenotype_modifier": +5.0,
                    "confidence_factor": 0.85,
                    "rarity_bonus": 2.0,
                    "feedback_boost": 1.05
                },
                "explanation": "..."
            }
        """
        base_score = self.severity_scores.get(severity, 55.0)
        phenotype_mod = self.phenotype_modifiers.get(phenotype, 0.0)
        confidence_factor = 0.70 + (0.30 * confidence)

        if population_frequency < 0.001:
            rarity_bonus = 8.0
        elif population_frequency < 0.01:
            rarity_bonus = 5.0
        elif population_frequency < 0.05:
            rarity_bonus = 2.0
        else:
            rarity_bonus = 0.0

        risk_score = (
            (base_score + phenotype_mod) * confidence_factor + rarity_bonus
        ) * feedback_boost
        risk_score = max(0.0, min(100.0, risk_score))

        risk_level = self.get_risk_level(risk_score)

        explanation = (
            f"Base severity '{severity}' ({base_score:.1f}) "
            f"+ phenotype '{phenotype}' ({phenotype_mod:+.1f}) "
            f"× confidence penalty ({confidence_factor:.2f}) "
            f"+ rarity bonus ({rarity_bonus:.1f}) "
            f"× feedback learning ({feedback_boost:.2f}) "
            f"= {risk_score:.1f}"
        )

        return {
            "risk_score": round(risk_score, 2),
            "risk_level": risk_level,
            "components": {
                "base_severity_score": base_score,
                "phenotype_modifier": phenotype_mod,
                "confidence_factor": round(confidence_factor, 3),
                "rarity_bonus": rarity_bonus,
                "feedback_boost": round(feedback_boost, 3),
            },
            "explanation": explanation,
        }


# ============================================================================
# Convenience Functions
# ============================================================================

def calculate_risk_score(
    severity: str,
    phenotype: str,
    confidence: float,
    population_frequency: float = 0.5,
    feedback_boost: float = 1.0,
) -> float:
    """
    Convenience function for calculating risk score.

    See RiskScoreCalculator.calculate_risk_score for details.
    """
    calculator = RiskScoreCalculator()
    return calculator.calculate_risk_score(
        severity=severity,
        phenotype=phenotype,
        confidence=confidence,
        population_frequency=population_frequency,
        feedback_boost=feedback_boost,
    )


def get_risk_level(risk_score: float) -> str:
    """Map risk score to categorical level."""
    calculator = RiskScoreCalculator()
    return calculator.get_risk_level(risk_score)
