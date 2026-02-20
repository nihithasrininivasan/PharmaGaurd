"""
Risk Engine - Evaluates pharmacogenomic risk based on drug-gene-phenotype combinations.
Enhanced with adaptive learning, numeric risk scoring, and confidence calibration.

Features:
- Continuous risk scoring (0-100 scale)
- Bayesian feedback learning from clinician corrections
- Calibrated confidence scores
- Structured recommendation generation
- Model versioning and drift detection
"""

from typing import Tuple, Optional, Dict, List
from pathlib import Path
from .models import (
    RiskAssessment,
    ClinicalRecommendation,
    DrugAssessment,
    DiplotypeResult,
    PatientProfile
)
from .cpic_loader import get_cpic_loader
from .confidence import ConfidenceBreakdown, ConfidenceCalculator
from .risk_scoring import RiskScoreCalculator
from .feedback_learning import (
    load_learning_priors,
    LearningPriors,
    get_diplotype_boost
)
from .recommendation_engine import RecommendationEngine
from .model_calibration import ConfidenceCalibrator
from .population_data import PopulationDataLoader


# ---------------------------------------------------------------------------
# Deterministic severity table
# ---------------------------------------------------------------------------

# Canonical risk label → severity mapping.
# These are the ONLY allowed risk labels in the system.
RISK_SEVERITY_TABLE: Dict[str, str] = {
    "Toxic":                       "critical",
    "Ineffective":                 "high",
    "Avoid":                       "critical",
    "Use Alternative":             "high",
    "Adjust Dosage":               "moderate",
    "Standard dosing recommended": "none",
    "Safe":                        "none",
    "Unknown":                     "moderate",
    "No specific CPIC recommendation": "moderate",
}


# Map CPIC recommendation keywords → canonical risk labels
# This avoids NLP interpretation — pattern matches on structured CPIC text
def _classify_risk_from_cpic_text(
    risk_text: str,
    implication_text: str = "",
    severity: str = "",
) -> str:
    """
    Map CPIC recommendation text to a canonical risk label.
    Uses deterministic keyword matching on structured fields only.
    Examines both short risk summary AND full implication text.
    """
    # Combine both fields for a broader match
    combined = (risk_text + " " + implication_text).lower()
    sev_lower = severity.lower() if severity else ""

    # ---------- Action keywords (always checked) ----------

    if any(kw in combined for kw in ("avoid", "contraindicated", "do not use")):
        return "Avoid"

    if any(kw in combined for kw in (
        "increased risk of toxicity", "life-threatening", "fatal",
        "severe toxicity",
    )):
        return "Toxic"

    if any(kw in combined for kw in ("lack of efficacy", "ineffective", "no therapeutic effect")):
        return "Ineffective"

    if any(kw in combined for kw in (
        "alternative antiplatelet", "alternative therapy", "consider alternative",
        "use an alternative",
    )):
        return "Use Alternative"

    if any(kw in combined for kw in (
        "reduce dose", "lower dose", "decreased dose",
        "dose reduction", "reduced starting dose",
        "20-50%", "25-50%",  # common CPIC dosing fractions
    )):
        return "Adjust Dosage"

    # ---------- Standard dosing (only when severity is NOT high/critical) ----------
    # CPIC texts for high-severity entries often mention "standard starting dose"
    # as a calculation reference (e.g., "Initiate at 20-50% of standard starting dose")
    # rather than as a recommendation to use standard dosing.
    if sev_lower not in ("high", "critical"):
        if any(kw in combined for kw in (
            "standard starting dose", "standard dose",
            "label recommended", "no change", "use standard",
            "no clinical intervention",
        )):
            return "Standard dosing recommended"

    # Fallback: derive risk label from stored severity
    severity_to_label = {
        "critical": "Use Alternative",
        "high": "Adjust Dosage",
        "moderate": "Adjust Dosage",
        "low": "Standard dosing recommended",
        "none": "Standard dosing recommended",
    }
    if sev_lower in severity_to_label:
        return severity_to_label[sev_lower]

    return "Unknown"


class RiskEngine:
    """
    Evaluates pharmacogenomic risk for drug-gene-phenotype combinations.

    Enhanced with:
    - Numeric risk scoring (0-100)
    - Feedback learning from clinician corrections
    - Calibrated confidence scores
    - Structured recommendations
    """

    # Drugs without their own CPIC alert files can reuse a related drug's
    # recommendations from the same gene pathway.
    _DRUG_ALIASES: Dict[str, str] = {
        "warfarin":       "warfarin",      # CYP2C9 — no alert file yet, but keep mapping
        "azathioprine":   "thioguanine",   # Both TPMT pathway
        "5-fluorouracil": "fluorouracil",  # Synonyms
    }

    # Model version
    MODEL_VERSION = "2.0.0"

    def __init__(
        self,
        learning_priors_path: Path = Path("data/learning_priors.json"),
        enable_feedback_learning: bool = True,
        enable_calibration: bool = True,
    ):
        # Core components
        self.loader = get_cpic_loader()
        self.confidence_calc = ConfidenceCalculator()

        # New enhanced components
        self.risk_scorer = RiskScoreCalculator()
        self.recommendation_engine = RecommendationEngine()

        # Feedback learning
        self.enable_feedback_learning = enable_feedback_learning
        if enable_feedback_learning:
            self.learning_priors = load_learning_priors(learning_priors_path)
        else:
            self.learning_priors = LearningPriors(genes={}, metadata={})

        # Confidence calibration
        self.enable_calibration = enable_calibration
        if enable_calibration:
            self.calibrator = ConfidenceCalibrator()
        else:
            self.calibrator = None

        # Population data for Bayesian priors
        try:
            self.population_loader = PopulationDataLoader()
        except Exception:
            self.population_loader = None

    def _resolve_drug(self, drug: str) -> str:
        """Resolve drug aliases and normalize casing."""
        d = drug.lower()
        return self._DRUG_ALIASES.get(d, d)

    def evaluate_risk(
        self,
        drug: str,
        gene: str,
        phenotype: str,
        diplotype: str,
        diplotype_confidence: float = 1.0,
        diplotype_confidence_breakdown: Optional[Dict] = None,
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """
        Evaluate risk for a specific drug-gene-phenotype combination.

        Returns:
            Tuple of (RiskAssessment, ClinicalRecommendation)
        """
        # Resolve drug aliases (e.g., azathioprine → thioguanine)
        resolved_drug = self._resolve_drug(drug)

        # Validate drug is supported
        if not self.loader.is_drug_supported(resolved_drug):
            return self._create_unsupported_drug_response(drug)

        # Validate gene matches drug
        expected_gene = self.loader.get_drug_gene(resolved_drug)
        if expected_gene and expected_gene != gene:
            return self._create_gene_mismatch_response(drug, gene, expected_gene)

        # Get recommendation for this drug-phenotype combination
        recommendation_data = self.loader.get_drug_recommendation_for_phenotype(resolved_drug, phenotype)

        if not recommendation_data:
            return self._create_no_recommendation_response(
                drug=drug,
                gene=gene,
                diplotype=diplotype,
                phenotype=phenotype,
                diplotype_confidence=diplotype_confidence,
                diplotype_confidence_breakdown=diplotype_confidence_breakdown,
            )

        # Map recommendation to risk assessment
        risk_assessment = self._create_risk_assessment(
            recommendation_data=recommendation_data,
            diplotype_confidence=diplotype_confidence,
            diplotype_confidence_breakdown=diplotype_confidence_breakdown,
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
        )

        clinical_recommendation = self._create_clinical_recommendation(
            recommendation_data=recommendation_data,
            drug=drug,
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            risk_assessment=risk_assessment,
        )

        return risk_assessment, clinical_recommendation

    def evaluate_drug_for_patient(
        self,
        drug: str,
        patient_profile: PatientProfile
    ) -> Optional[DrugAssessment]:
        """Evaluate a drug against a patient's complete pharmacogenomic profile."""
        gene = self.loader.get_drug_gene(drug)
        if not gene:
            return None

        diplotype_result = patient_profile.diplotypes.get(gene)
        if not diplotype_result:
            return None

        risk, recommendation = self.evaluate_risk(
            drug=drug,
            gene=gene,
            phenotype=diplotype_result.phenotype,
            diplotype=diplotype_result.diplotype,
            diplotype_confidence=diplotype_result.confidence,
            diplotype_confidence_breakdown=diplotype_result.confidence_breakdown,
        )

        return DrugAssessment(
            drug=drug,
            gene=gene,
            diplotype=diplotype_result.diplotype,
            phenotype=diplotype_result.phenotype,
            risk=risk,
            recommendation=recommendation
        )

    def evaluate_multiple_drugs(
        self,
        drugs: List[str],
        patient_profile: PatientProfile
    ) -> List[DrugAssessment]:
        """Evaluate multiple drugs for a patient."""
        assessments = []
        for drug in drugs:
            assessment = self.evaluate_drug_for_patient(drug, patient_profile)
            if assessment:
                assessments.append(assessment)
        return assessments

    # ===== Helper Methods =====

    def _create_risk_assessment(
        self,
        recommendation_data: Dict,
        diplotype_confidence: float,
        diplotype_confidence_breakdown: Optional[Dict] = None,
        gene: str = "",
        diplotype: str = "",
        phenotype: str = "",
    ) -> RiskAssessment:
        """
        Create RiskAssessment from CPIC recommendation data.

        Enhanced with:
        - Numeric risk score (0-100)
        - Feedback learning integration
        - Population frequency adjustment
        """
        # Get severity directly from CPIC cache (deterministic)
        severity = recommendation_data.get('severity', 'moderate')

        # Classify risk using both short text, full implication, and stored severity
        risk_text = recommendation_data.get('risk', '')
        implication_text = recommendation_data.get('implication', '')
        risk_label = _classify_risk_from_cpic_text(risk_text, implication_text, severity)

        # Use CPIC severity if valid, otherwise derive from risk label
        if severity in ('none', 'low', 'moderate', 'high', 'critical'):
            final_severity = severity
        else:
            final_severity = RISK_SEVERITY_TABLE.get(risk_label, "moderate")

        # Build confidence breakdown for the risk assessment
        bd = ConfidenceBreakdown()

        # Start from diplotype confidence breakdown if available
        if diplotype_confidence_breakdown:
            bd.variant_quality = diplotype_confidence_breakdown.get('variant_quality', 1.0)
            bd.allele_coverage = diplotype_confidence_breakdown.get('allele_coverage', 1.0)
            bd.phase_resolution = diplotype_confidence_breakdown.get('phase_resolution', 1.0)
            bd.cnv_evaluation = diplotype_confidence_breakdown.get('cnv_evaluation', 1.0)
            bd.diplotype_determinism = diplotype_confidence_breakdown.get('diplotype_determinism', 1.0)

        # CPIC applicability is always 1.0 here (we found a rule)
        bd.cpic_applicability = 1.0

        # Final confidence = min(all components) — never exceeds weakest link
        # Cap by diplotype_confidence so we never report higher confidence than the genotype call
        final_confidence = min(bd.final, diplotype_confidence)

        # Apply calibration if enabled
        if self.enable_calibration and self.calibrator:
            calibration_factor = self.calibrator.get_calibration_factor(final_confidence)
            final_confidence = final_confidence * calibration_factor

        # Get feedback learning boost
        feedback_boost = 1.0
        if self.enable_feedback_learning and gene and diplotype:
            feedback_boost = get_diplotype_boost(
                self.learning_priors,
                gene,
                diplotype
            )

        # Get population frequency for Bayesian adjustment
        population_freq = 0.5  # Default
        if self.population_loader and gene and diplotype:
            try:
                pop_data = self.population_loader.get_diplotype_frequency(
                    gene,
                    diplotype
                )
                if pop_data:
                    population_freq = pop_data.get('frequency', 0.5)
            except Exception:
                pass

        # Calculate numeric risk score (0-100)
        risk_score = self.risk_scorer.calculate_risk_score(
            severity=final_severity,
            phenotype=phenotype,
            confidence=final_confidence,
            population_frequency=population_freq,
            feedback_boost=feedback_boost,
        )

        # Get risk level from score
        risk_level = self.risk_scorer.get_risk_level(risk_score)

        return RiskAssessment(
            risk_label=risk_label,
            confidence_score=max(0.0, min(1.0, final_confidence)),
            severity=final_severity,
            confidence_breakdown=bd.to_dict(),
            risk_score=risk_score,
            risk_level=risk_level,
        )

    def _create_clinical_recommendation(
        self,
        recommendation_data: Dict,
        drug: str = "",
        gene: str = "",
        diplotype: str = "",
        phenotype: str = "",
        risk_assessment: Optional[RiskAssessment] = None,
    ) -> ClinicalRecommendation:
        """
        Create ClinicalRecommendation from CPIC recommendation data.

        Enhanced with structured recommendation generation.
        """
        # Extract CPIC data
        cpic_text = recommendation_data.get('risk', 'No specific recommendation')
        cpic_implication = recommendation_data.get('implication', 'Standard considerations apply')
        cpic_url = recommendation_data.get('url')

        # If we have risk_assessment with risk_score, generate structured recommendation
        if risk_assessment and risk_assessment.risk_score is not None:
            try:
                structured_rec = self.recommendation_engine.generate_recommendation(
                    risk_score=risk_assessment.risk_score,
                    risk_level=risk_assessment.risk_level or risk_assessment.severity,
                    drug=drug,
                    gene=gene,
                    diplotype=diplotype,
                    phenotype=phenotype,
                    confidence_score=risk_assessment.confidence_score,
                    cpic_text=cpic_text,
                    cpic_implication=cpic_implication,
                    cpic_url=cpic_url,
                )

                # Convert structured recommendation to ClinicalRecommendation
                # Combine primary action with monitoring
                full_text = (
                    f"{structured_rec.clinical_recommendation.primary_action}\n\n"
                    f"Monitoring: {structured_rec.clinical_recommendation.monitoring}"
                )

                # Add dosing guidance if available
                if structured_rec.clinical_recommendation.dosing_guidance:
                    full_text = (
                        f"{full_text}\n\n"
                        f"Dosing: {structured_rec.clinical_recommendation.dosing_guidance}"
                    )

                return ClinicalRecommendation(
                    text=full_text,
                    implication=cpic_implication,
                    recommendation_url=cpic_url,
                )
            except Exception as e:
                # Fallback to simple CPIC text if recommendation generation fails
                pass

        # Fallback: simple CPIC text
        return ClinicalRecommendation(
            text=cpic_text,
            implication=cpic_implication,
            recommendation_url=cpic_url,
        )

    def _create_unsupported_drug_response(
        self, drug: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """Create response for unsupported drug."""
        risk_score = 0.0
        risk_level = "none"
        risk = RiskAssessment(
            risk_label=f"Drug '{drug}' not in CPIC database",
            confidence_score=0.0,
            severity="none",
            confidence_breakdown=None,
            risk_score=risk_score,
            risk_level=risk_level,
        )

        recommendation = ClinicalRecommendation(
            text=f"No pharmacogenomic guidance available for {drug}",
            implication="Drug not currently supported by CPIC guidelines",
            recommendation_url=None
        )

        return risk, recommendation

    def _create_gene_mismatch_response(
        self, drug: str, provided_gene: str, expected_gene: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """Create response for gene-drug mismatch."""
        risk_score = 0.0
        risk_level = "none"
        risk = RiskAssessment(
            risk_label=f"Gene mismatch: expected {expected_gene}, got {provided_gene}",
            confidence_score=0.0,
            severity="none",
            confidence_breakdown=None,
            risk_score=risk_score,
            risk_level=risk_level,
        )

        recommendation = ClinicalRecommendation(
            text=f"Incorrect gene provided for {drug}. Expected {expected_gene}.",
            implication=f"Drug {drug} is primarily metabolized by {expected_gene}",
            recommendation_url=None
        )

        return risk, recommendation

    def _create_no_recommendation_response(
        self,
        drug: str,
        gene: str,
        diplotype: str,
        phenotype: str,
        diplotype_confidence: float,
        diplotype_confidence_breakdown: Optional[Dict] = None,
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """
        Create response when no specific CPIC recommendation exists for this phenotype.
        """
        is_normal = phenotype in (
            "NM", "Normal Metabolizer",
            "Normal Function", "Increased Function",
            "Extensive Metabolizer",
            "Normal Activity",
        )

        # Build confidence breakdown
        bd = ConfidenceBreakdown()
        if diplotype_confidence_breakdown:
            bd.variant_quality = diplotype_confidence_breakdown.get('variant_quality', 1.0)
            bd.allele_coverage = diplotype_confidence_breakdown.get('allele_coverage', 1.0)
            bd.phase_resolution = diplotype_confidence_breakdown.get('phase_resolution', 1.0)
            bd.cnv_evaluation = diplotype_confidence_breakdown.get('cnv_evaluation', 1.0)
            bd.diplotype_determinism = diplotype_confidence_breakdown.get('diplotype_determinism', 1.0)

        if is_normal:
            # Normal metabolizer — standard dosing
            bd.cpic_applicability = 1.0
            raw_final = bd.final
            final_conf = min(raw_final, diplotype_confidence)
            severity = "none"
            risk_score = self.risk_scorer.calculate_risk_score(
                severity=severity,
                phenotype=phenotype,
                confidence=final_conf,
                population_frequency=0.5,
                feedback_boost=get_diplotype_boost(self.learning_priors, gene, diplotype) if self.enable_feedback_learning else 1.0,
            )
            risk_level = self.risk_scorer.get_risk_level(risk_score)
            risk = RiskAssessment(
                risk_label="Standard dosing recommended",
                confidence_score=max(0.0, min(1.0, final_conf)),
                severity=severity,
                confidence_breakdown=bd.to_dict(),
                risk_score=risk_score,
                risk_level=risk_level,
            )
            recommendation = ClinicalRecommendation(
                text=f"Use standard {drug} dosing guidelines",
                implication="Normal drug metabolism expected",
                recommendation_url=None
            )
        else:
            # Phenotype exists but no CPIC recommendation — penalise
            self.confidence_calc.apply_cpic_penalties(
                bd,
                has_cpic_rule=False,
                phenotype_is_indeterminate=(phenotype in ("Indeterminate", "Unknown")),
            )
            raw_final = bd.final
            final_conf = min(raw_final, diplotype_confidence)
            severity = "moderate"
            risk_score = self.risk_scorer.calculate_risk_score(
                severity=severity,
                phenotype=phenotype,
                confidence=final_conf,
                population_frequency=0.5,
                feedback_boost=get_diplotype_boost(self.learning_priors, gene, diplotype) if self.enable_feedback_learning else 1.0,
            )
            risk_level = self.risk_scorer.get_risk_level(risk_score)
            risk = RiskAssessment(
                risk_label="No specific CPIC recommendation",
                confidence_score=max(0.0, min(1.0, final_conf)),
                severity=severity,
                confidence_breakdown=bd.to_dict(),
                risk_score=risk_score,
                risk_level=risk_level,
            )
            recommendation = ClinicalRecommendation(
                text=f"No CPIC recommendation for {drug} in {phenotype} individuals. Exercise clinical judgment.",
                implication=f"Limited guidance for {phenotype} phenotype with {drug}",
                recommendation_url=None
            )

        return risk, recommendation


def create_risk_engine() -> RiskEngine:
    """Factory function to create a RiskEngine instance."""
    return RiskEngine()
