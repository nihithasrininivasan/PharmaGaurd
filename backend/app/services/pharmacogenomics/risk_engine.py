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
from .recommendation_engine import RecommendationEngine
from .population_data import PopulationDataLoader
from .feedback_learning import load_learning_priors, LearningPriors, get_diplotype_boost
from .model_calibration import ConfidenceCalibrator
from .pharmgkb_loader import get_pharmgkb_loader as get_pharmgkb, harmonize_annotation_associations


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

    # Drug name normalization — map common synonyms to canonical names.
    # NOTE: Do NOT map to a DIFFERENT drug (e.g. azathioprine→thioguanine)
    # as this corrupts gene_drug_confirmation.drug integrity.
    _DRUG_ALIASES: Dict[str, str] = {
        "warfarin":       "warfarin",      # CYP2C9 — identity mapping
        "azathioprine":   "azathioprine",  # TPMT — must NOT alias to thioguanine
        "5-fluorouracil": "fluorouracil",  # Synonyms
    }

    # Canonical drug→gene mapping (works even without CPIC alert Excel files).
    # This is the authoritative mapping used by evaluate_drug_for_patient()
    # when the CPIC cache doesn't contain the drug.
    _GENE_DRUG_MAP: Dict[str, str] = {
        "codeine":        "CYP2D6",
        "clopidogrel":    "CYP2C19",
        "warfarin":       "CYP2C9",
        "simvastatin":    "SLCO1B1",
        "azathioprine":   "TPMT",
        "fluorouracil":   "DPYD",
        "thioguanine":    "TPMT",
        "5-fluorouracil": "DPYD",
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

        try:
            self.population_loader = PopulationDataLoader()
        except Exception:
            self.population_loader = None

        # PharmGKB dataset loader for variant annotation, evidence scoring,
        # gene-drug confirmation, and clinical annotation linking
        try:
            self.pharmgkb = get_pharmgkb()
        except Exception:
            self.pharmgkb = None

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
        # Resolve drug aliases (normalize synonyms, e.g. 5-fluorouracil → fluorouracil)
        resolved_drug = self._resolve_drug(drug)

        # Validate drug is supported (check PharmGKB for Level 1A/1B evidence)
        # This aggregates across all genes: warfarin is supported if it has
        # 1A/1B evidence for ANY of CYP2C9, VKORC1, CYP4F2, etc.
        drug_is_supported = False
        if self.pharmgkb:
            drug_is_supported = self.pharmgkb.is_drug_supported(resolved_drug)
        else:
            # Fallback to CPIC loader if PharmGKB not available
            drug_is_supported = self.loader.is_drug_supported(resolved_drug)

        if not drug_is_supported:
            return self._create_unsupported_drug_response(drug)

        # Validate gene matches drug (only for CPIC guideline files)
        expected_gene = self.loader.get_drug_gene(resolved_drug)
        if expected_gene and expected_gene != gene:
            return self._create_gene_mismatch_response(drug, gene, expected_gene)

        # ---- PharmGKB Gene-Drug Confirmation ----
        gene_drug_confirmed = True
        gene_drug_confirmation = None
        evidence_info = None
        clinical_annotations = None
        knowledge_confidence = 1.0  # Default: assume well-studied

        if self.pharmgkb:
            # Gene-drug integrity check: confirm_gene_drug_pair returns
            # the drug name as stored in PharmGKB. If it doesn't match
            # the input drug (case-insensitive), reject with hard error.
            _confirmation_preview = self.pharmgkb.confirm_gene_drug_pair(gene, resolved_drug)
            _confirmed_drug = _confirmation_preview.get("drug", "").strip().lower()
            if _confirmation_preview.get("confirmed") and _confirmed_drug != resolved_drug.lower():
                return self._create_gene_drug_integrity_error(drug, gene, _confirmed_drug)
            # Gate 4: Confirm gene-drug pair from relationships.tsv
            gene_drug_confirmation = self.pharmgkb.confirm_gene_drug_pair(gene, resolved_drug)
            if not gene_drug_confirmation.get("confirmed", False):
                gene_drug_confirmed = False
                # Build breakdown showing WHY automation is blocked
                bd = ConfidenceBreakdown()
                bd.gene_drug_confirmed = False
                bd.knowledge_confidence = 0.0
                bd.diplotype_determinism = 0.0
                bd.allele_coverage = 0.3
                bd.cnv_evaluation = 0.5
                bd.penalties_applied.append(
                    f"Gene-drug pair ({gene}, {drug}) not found in PharmGKB relationships"
                )
                auto_status = bd.get_automation_status()
                risk = RiskAssessment(
                    risk_label="Unsupported in current knowledge base",
                    confidence_score=bd.final,
                    severity="undetermined",
                    confidence_breakdown=bd.to_dict(),
                    gene_drug_confirmation=gene_drug_confirmation,
                    automation_status=auto_status,
                )
                recommendation = ClinicalRecommendation(
                    text=f"No PharmGKB gene-drug relationship found for {drug} and {gene}. "
                         f"Automated recommendation blocked. "
                         f"Blocked gate: Gene-drug pair not confirmed.",
                    implication=f"Gene-drug pair ({gene}, {drug}) absent from PharmGKB relationships dataset",
                    recommendation_url=None,
                )
                return risk, recommendation

            # Evidence level → knowledge_confidence (Layer 1)
            evidence_info = self.pharmgkb.get_evidence_level(gene, resolved_drug)
            knowledge_confidence = evidence_info.get("confidence_weight", 1.0)
            # Mark role explicitly
            evidence_info["role"] = "knowledge_confidence_only"

            # Clinical annotations (deduplicated)
            clinical_annotations = self.pharmgkb.get_clinical_annotations(gene, resolved_drug)

            # Harmonize clinical annotation associations with top-level classification
            # This ensures hierarchical consistency: if top-level = "established",
            # individual annotations should not contradict with "ambiguous"
            if clinical_annotations and gene_drug_confirmation:
                top_level_association = gene_drug_confirmation.get("association", "")
                clinical_annotations = harmonize_annotation_associations(
                    clinical_annotations,
                    top_level_association
                )

        # Block CPIC mapping when evidence is insufficient
        # Unresolved diplotype / Indeterminate phenotype → no CPIC lookup
        if diplotype in ("Unresolved", "Indeterminate", "Unknown") or \
           phenotype in ("Indeterminate", "Unknown"):
            return self._create_no_recommendation_response(
                drug=drug,
                gene=gene,
                diplotype=diplotype,
                phenotype=phenotype,
                diplotype_confidence=diplotype_confidence,
                diplotype_confidence_breakdown=diplotype_confidence_breakdown,
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_info=evidence_info,
                clinical_annotations=clinical_annotations,
                knowledge_confidence=knowledge_confidence,
                gene_drug_confirmed=gene_drug_confirmed,
            )

        # Calculate activity score if possible (for drugs like codeine)
        activity_score = None
        if gene and diplotype:
            try:
                activity_score = self.loader.calculate_total_activity_score(gene, diplotype)
            except Exception:
                pass

        # Get recommendation for this drug-phenotype combination
        recommendation_data = self.loader.get_drug_recommendation_for_phenotype(
            resolved_drug, phenotype, activity_score=activity_score
        )

        if not recommendation_data:
            return self._create_no_recommendation_response(
                drug=drug,
                gene=gene,
                diplotype=diplotype,
                phenotype=phenotype,
                diplotype_confidence=diplotype_confidence,
                diplotype_confidence_breakdown=diplotype_confidence_breakdown,
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_info=evidence_info,
                clinical_annotations=clinical_annotations,
                knowledge_confidence=knowledge_confidence,
                gene_drug_confirmed=gene_drug_confirmed,
            )

        # Map recommendation to risk assessment
        risk_assessment = self._create_risk_assessment(
            recommendation_data=recommendation_data,
            diplotype_confidence=diplotype_confidence,
            diplotype_confidence_breakdown=diplotype_confidence_breakdown,
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            knowledge_confidence=knowledge_confidence,
            gene_drug_confirmed=gene_drug_confirmed,
        )

        # Attach PharmGKB data layers to the risk assessment
        if gene_drug_confirmation:
            risk_assessment.gene_drug_confirmation = gene_drug_confirmation
        if evidence_info:
            risk_assessment.evidence_level = evidence_info
        if clinical_annotations:
            risk_assessment.clinical_annotations = clinical_annotations

        # ---- 4-Gate Automation Check ----
        # Read from model field (single location, no duplication)
        auto_status = risk_assessment.automation_status or {}
        if not auto_status.get("allowed", True):
            blocked_reasons = auto_status.get("blocked_reasons", [])
            blocked_text = "; ".join(blocked_reasons) if blocked_reasons else "Unknown gate failure"
            recommendation = ClinicalRecommendation(
                text=(
                    f"Automated CPIC-based dosing guidance for {drug} ({gene}) is blocked. "
                    f"Blocked gates: {blocked_text}. "
                    f"Exercise clinical judgment and consider comprehensive pharmacogenomic testing."
                ),
                implication=f"Automation blocked: {blocked_text}",
                recommendation_url=None,
            )
            return risk_assessment, recommendation

        clinical_recommendation = self._create_clinical_recommendation(
            recommendation_data=recommendation_data,
            drug=drug,
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            risk_assessment=risk_assessment,
        )
        
        # Runtime Invariant Check
        if not self._validate_invariants(risk_assessment, phenotype, drug, gene):
             return self._create_no_recommendation_response(
                 drug=drug,
                 gene=gene,
                 diplotype=diplotype,
                 phenotype=phenotype,
                 diplotype_confidence=diplotype_confidence,
                 diplotype_confidence_breakdown=diplotype_confidence_breakdown,
                 gene_drug_confirmation=gene_drug_confirmation,
                 evidence_info=evidence_info,
                 clinical_annotations=clinical_annotations,
                 knowledge_confidence=knowledge_confidence,
                 gene_drug_confirmed=gene_drug_confirmed,
             )

        return risk_assessment, clinical_recommendation
        
    def _validate_invariants(
        self,
        risk: RiskAssessment,
        phenotype: str,
        drug: str,
        gene: str
    ) -> bool:
        """
        Enforce critical safety invariants.
        Returns False if invariant violated (block output).
        """
        # Invariant 1: Normal Metabolizer != Critical Risk
        # (Unless explicitly Ultrarapid, which is distinctly named)
        is_normal = phenotype in (
            "NM", "Normal Metabolizer",
            "Normal Function", "Extended Metabolizer" # EM is sometimes Normal
        )
        if is_normal and risk.severity == "critical":
            print(f"CRITICAL INVARIANT FAILURE: {drug} ({gene}) - Phenotype '{phenotype}' but Severity 'critical'. blocking.")
            return False
            
        return True

    def _get_drug_gene(self, drug: str) -> Optional[str]:
        """Get primary gene for a drug, checking CPIC cache then canonical map."""
        resolved = self._resolve_drug(drug)
        gene = self.loader.get_drug_gene(resolved)
        if gene:
            return gene
        return self._GENE_DRUG_MAP.get(resolved)

    def evaluate_drug_for_patient(
        self,
        drug: str,
        patient_profile: PatientProfile
    ) -> Optional[DrugAssessment]:
        """Evaluate a drug against a patient's complete pharmacogenomic profile."""
        gene = self._get_drug_gene(drug)
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

    def calculate_confidence_score(
        self, base: float, coverage: float, ambiguity: float
    ) -> float:
        """Calculate confidence score from base, coverage, and ambiguity factors."""
        return max(0.0, min(1.0, base * coverage * ambiguity))

    @staticmethod
    def map_severity_level(implication_text: str) -> str:
        """Map CPIC implication text to a severity level."""
        text = implication_text.lower()
        if any(kw in text for kw in (
            "contraindicated", "life-threatening", "avoid use", "fatal",
        )):
            return "critical"
        if any(kw in text for kw in (
            "high risk", "therapeutic failure", "ineffective",
        )):
            return "high"
        if any(kw in text for kw in (
            "dose adjustment", "alternative therapy", "consider",
            "reduce dose", "decreased",
        )):
            return "moderate"
        if any(kw in text for kw in (
            "informative", "minor", "low",
        )):
            return "low"
        if any(kw in text for kw in (
            "normal", "standard", "no change",
        )):
            return "none"
        return "moderate"

    # ===== Helper Methods =====

    def _create_risk_assessment(
        self,
        recommendation_data: Dict,
        diplotype_confidence: float,
        diplotype_confidence_breakdown: Optional[Dict] = None,
        gene: str = "",
        diplotype: str = "",
        phenotype: str = "",
        knowledge_confidence: float = 1.0,
        gene_drug_confirmed: bool = True,
    ) -> RiskAssessment:
        """
        Create RiskAssessment from CPIC recommendation data.

        Uses two-axis confidence:
          - phenotype_confidence: genotype_quality × diplotype_determinism
          - classification_confidence: how confident we are in the LABEL
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

        # Build confidence breakdown
        bd = ConfidenceBreakdown()

        # Knowledge confidence (external evidence)
        bd.knowledge_confidence = knowledge_confidence

        # Gene-drug confirmation
        bd.gene_drug_confirmed = gene_drug_confirmed

        # Genotype components (from diplotype resolution)
        if diplotype_confidence_breakdown:
            bd.variant_quality = diplotype_confidence_breakdown.get('variant_quality', 1.0)
            bd.allele_coverage = diplotype_confidence_breakdown.get('allele_coverage', 1.0)
            bd.cnv_evaluation = diplotype_confidence_breakdown.get('cnv_evaluation', 1.0)
            bd.genome_build_validity = diplotype_confidence_breakdown.get('genome_build_validity', 1.0)
            bd.diplotype_determinism = diplotype_confidence_breakdown.get('diplotype_determinism', 1.0)

        # CPIC applicability is 1.0 here (we found a rule)
        bd.cpic_applicability = 1.0

        # Classification confidence (label correctness)
        final_confidence = bd.final

        # Automation status (single location: on model only)
        auto_status = bd.get_automation_status()

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
            automation_status=auto_status,
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
        """
        Create response for unsupported drug.

        This is ONLY called when the drug has no Level 1A or 1B evidence
        in the PharmGKB dataset across ALL genes.
        """
        risk = RiskAssessment(
            risk_label="Drug not currently supported by CPIC guidelines",
            confidence_score=0.0,
            severity="none",
            confidence_breakdown=None,
        )

        recommendation = ClinicalRecommendation(
            text="Drug not currently supported by CPIC guidelines",
            implication=f"No annotation with evidence_level 1A or 1B exists for {drug}, "
                        f"and no gene_drug_confirmation confirms support",
            recommendation_url=None
        )

        return risk, recommendation

    def _create_gene_mismatch_response(
        self, drug: str, provided_gene: str, expected_gene: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """Create response for gene-drug mismatch."""
        risk = RiskAssessment(
            risk_label=f"Gene mismatch: expected {expected_gene}, got {provided_gene}",
            confidence_score=0.0,
            severity="none",
            confidence_breakdown=None,
        )

        recommendation = ClinicalRecommendation(
            text=f"Incorrect gene provided for {drug}. Expected {expected_gene}.",
            implication=f"Drug {drug} is primarily metabolized by {expected_gene}",
            recommendation_url=None
        )

        return risk, recommendation

    def _create_gene_drug_integrity_error(
        self, drug: str, gene: str, confirmed_drug: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """
        Hard validation error when gene_drug_confirmation.drug doesn't match input drug.

        This prevents analysis from proceeding with corrupted drug identity
        (e.g. azathioprine being confirmed as thioguanine).
        """
        risk = RiskAssessment(
            risk_label=f"Gene-drug integrity error: expected {drug.lower()}, got {confirmed_drug}",
            confidence_score=0.0,
            severity="none",
            confidence_breakdown=None,
        )

        recommendation = ClinicalRecommendation(
            text=(
                f"Gene-drug integrity check failed for {drug}/{gene}. "
                f"Confirmed drug '{confirmed_drug}' does not match input drug '{drug.lower()}'. "
                f"Analysis halted."
            ),
            implication="Gene-drug identity mismatch — cannot proceed with analysis",
            recommendation_url=None,
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
        gene_drug_confirmation: Optional[Dict] = None,
        evidence_info: Optional[Dict] = None,
        clinical_annotations: Optional[List] = None,
        knowledge_confidence: float = 1.0,
        gene_drug_confirmed: bool = True,
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """
        Create response when no specific CPIC recommendation exists.

        Uses two-axis confidence: phenotype_confidence and
        classification_confidence. Genotype components are capped
        when diplotype is unresolved.
        """
        is_normal = phenotype in (
            "NM", "Normal Metabolizer",
            "Normal Function", "Increased Function",
            "Extensive Metabolizer",
            "Normal Activity",
        )

        # Build confidence breakdown
        bd = ConfidenceBreakdown()
        bd.knowledge_confidence = knowledge_confidence
        bd.gene_drug_confirmed = gene_drug_confirmed

        if diplotype_confidence_breakdown:
            bd.variant_quality = diplotype_confidence_breakdown.get('variant_quality', 1.0)
            bd.allele_coverage = diplotype_confidence_breakdown.get('allele_coverage', 1.0)
            bd.cnv_evaluation = diplotype_confidence_breakdown.get('cnv_evaluation', 1.0)
            bd.genome_build_validity = diplotype_confidence_breakdown.get('genome_build_validity', 1.0)
            bd.diplotype_determinism = diplotype_confidence_breakdown.get('diplotype_determinism', 1.0)

        # Handle Unresolved diplotype — insufficient evidence for ANY phenotype
        is_unresolved = diplotype in ("Unresolved",) or (
            phenotype in ("Indeterminate", "Unknown")
            and diplotype not in ("Indeterminate", "Unknown")
        )

        if is_unresolved or diplotype == "Unresolved":
            # Phenotype unresolved → diplotype_determinism = 0
            bd.diplotype_determinism = 0.0
            # Cap genotype components to reflect unresolved state
            bd.allele_coverage = min(bd.allele_coverage, 0.3)
            bd.cnv_evaluation = min(bd.cnv_evaluation, 0.5)
            PENALTY_GUIDELINE_BLOCKED = 0.30
            bd.cpic_applicability = max(0.0, 1.0 - PENALTY_GUIDELINE_BLOCKED)
            bd.penalties_applied.append(
                f"Guideline mapping blocked due to unresolved phenotype (−{PENALTY_GUIDELINE_BLOCKED:.2f})"
            )
            auto_status = bd.get_automation_status()
            risk = RiskAssessment(
                risk_label="Supported drug — insufficient genotype resolution for recommendation",
                confidence_score=bd.final,
                severity="undetermined",
                confidence_breakdown=bd.to_dict(),
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_level=evidence_info,
                clinical_annotations=clinical_annotations,
                automation_status=auto_status,
            )
            recommendation = ClinicalRecommendation(
                text=(
                    f"Phenotype could not be determined due to insufficient {gene} "
                    f"genetic coverage. Automated CPIC-based dosing guidance is therefore "
                    f"blocked. Recommend comprehensive {gene} testing including CNV assessment."
                ),
                implication="Genetic data insufficient for phenotype determination — automation blocked",
                recommendation_url=None,
            )
            return risk, recommendation

        if is_normal:
            # Normal metabolizer — standard dosing
            bd.cpic_applicability = 0.90
            bd.penalties_applied.append(
                "No specific CPIC rule for known Normal phenotype (−0.10)"
            )
            auto_status = bd.get_automation_status()
            risk = RiskAssessment(
                risk_label="Standard dosing recommended",
                confidence_score=bd.final,
                severity="none",
                confidence_breakdown=bd.to_dict(),
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_level=evidence_info,
                clinical_annotations=clinical_annotations,
                automation_status=auto_status,
            )
            recommendation = ClinicalRecommendation(
                text=f"Use standard {drug} dosing guidelines",
                implication="Normal drug metabolism expected",
                recommendation_url=None
            )
        elif phenotype in ("Indeterminate", "Unknown"):
            # Phenotype unresolvable → diplotype_determinism = 0
            bd.diplotype_determinism = 0.0
            bd.allele_coverage = min(bd.allele_coverage, 0.3)
            bd.cnv_evaluation = min(bd.cnv_evaluation, 0.5)
            PENALTY_GUIDELINE_BLOCKED = 0.30
            bd.cpic_applicability = max(0.0, 1.0 - PENALTY_GUIDELINE_BLOCKED)
            bd.penalties_applied.append(
                f"Guideline mapping blocked due to unresolved phenotype (−{PENALTY_GUIDELINE_BLOCKED:.2f})"
            )
            auto_status = bd.get_automation_status()
            risk = RiskAssessment(
                risk_label="Indeterminate Phenotype",
                confidence_score=bd.final,
                severity="undetermined",
                confidence_breakdown=bd.to_dict(),
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_level=evidence_info,
                clinical_annotations=clinical_annotations,
                automation_status=auto_status,
            )
            recommendation = ClinicalRecommendation(
                text=(
                    f"Phenotype for {gene} could not be determined. "
                    f"Automated CPIC-based dosing guidance is blocked. "
                    f"Consult a pharmacogenomics specialist."
                ),
                implication="Genetic data insufficient to predict drug response — automation blocked",
                recommendation_url=None
            )
        else:
            # Phenotype exists but no CPIC alert Excel file for this drug.
            # Generate clinically correct risk label based on phenotype.
            is_poor = phenotype in ("PM", "Poor Metabolizer")
            is_intermediate = phenotype in ("IM", "Intermediate Metabolizer")
            is_ultrarapid = phenotype in ("UM", "Ultrarapid Metabolizer")
            is_rapid = phenotype in ("RM", "Rapid Metabolizer")

            if is_poor or is_ultrarapid:
                # High-risk phenotype → clinically significant
                risk_label, severity = self._phenotype_risk_for_drug(
                    drug, gene, phenotype, "high"
                )
                bd.cpic_applicability = 0.85
                bd.penalties_applied.append(
                    "No CPIC alert file — using phenotype-driven risk classification (−0.15)"
                )
            elif is_intermediate:
                # Moderate-risk phenotype
                risk_label, severity = self._phenotype_risk_for_drug(
                    drug, gene, phenotype, "moderate"
                )
                bd.cpic_applicability = 0.85
                bd.penalties_applied.append(
                    "No CPIC alert file — using phenotype-driven risk classification (−0.15)"
                )
            else:
                # Other phenotype (Rapid, etc.) — moderate by default
                risk_label = f"No specific CPIC recommendation"
                severity = "moderate"
                self.confidence_calc.apply_cpic_penalties(
                    bd,
                    has_cpic_rule=False,
                    phenotype_is_indeterminate=False,
                )

            auto_status = bd.get_automation_status()
            risk = RiskAssessment(
                risk_label=risk_label,
                confidence_score=bd.final,
                severity=severity,
                confidence_breakdown=bd.to_dict(),
                gene_drug_confirmation=gene_drug_confirmation,
                evidence_level=evidence_info,
                clinical_annotations=clinical_annotations,
                automation_status=auto_status,
            )
            recommendation = self._phenotype_recommendation_for_drug(
                drug, gene, phenotype
            )

        return risk, recommendation

    def _phenotype_risk_for_drug(
        self, drug: str, gene: str, phenotype: str, base_severity: str
    ) -> Tuple[str, str]:
        """
        Determine risk label and severity for a resolved phenotype
        when no CPIC alert Excel file exists.

        Uses deterministic phenotype → risk mapping based on
        pharmacogenomic knowledge.
        """
        is_poor = phenotype in ("PM", "Poor Metabolizer")
        is_intermediate = phenotype in ("IM", "Intermediate Metabolizer")
        is_ultrarapid = phenotype in ("UM", "Ultrarapid Metabolizer")

        drug_lower = drug.lower()

        # Warfarin: CYP2C9 PM/IM → increased bleeding risk, reduce dose
        if drug_lower == "warfarin":
            if is_poor:
                return "Adjust Dosage", "high"
            elif is_intermediate:
                return "Adjust Dosage", "moderate"

        # Clopidogrel: CYP2C19 PM → reduced activation, use alternative
        if drug_lower == "clopidogrel":
            if is_poor:
                return "Use Alternative", "high"
            elif is_intermediate:
                return "Adjust Dosage", "moderate"

        # Generic phenotype-based classification
        if is_poor or is_ultrarapid:
            return "Adjust Dosage", base_severity
        elif is_intermediate:
            return "Adjust Dosage", base_severity

        return "Unknown", base_severity

    def _phenotype_recommendation_for_drug(
        self, drug: str, gene: str, phenotype: str
    ) -> ClinicalRecommendation:
        """
        Generate clinically correct recommendation for a resolved phenotype
        when no CPIC alert Excel file exists.
        """
        is_poor = phenotype in ("PM", "Poor Metabolizer")
        is_intermediate = phenotype in ("IM", "Intermediate Metabolizer")
        is_ultrarapid = phenotype in ("UM", "Ultrarapid Metabolizer")

        drug_lower = drug.lower()

        # ── Warfarin ──
        if drug_lower == "warfarin":
            if is_poor:
                return ClinicalRecommendation(
                    text=(
                        f"Reduce {drug} dose significantly. {gene} Poor Metabolizer status "
                        f"leads to reduced metabolism and increased drug exposure. "
                        f"Consider alternative anticoagulant or reduce dose by ≥50%."
                    ),
                    implication=(
                        f"Increased risk of bleeding due to reduced {gene} metabolism of warfarin. "
                        f"Consider pharmacogenomic-guided dosing."
                    ),
                    recommendation_url="https://cpicpgx.org/guidelines/",
                )
            elif is_intermediate:
                return ClinicalRecommendation(
                    text=(
                        f"Consider reducing {drug} dose. {gene} Intermediate Metabolizer "
                        f"status may lead to moderately reduced metabolism."
                    ),
                    implication=(
                        f"Moderately increased risk of bleeding due to reduced {gene} "
                        f"metabolism of warfarin. Monitor INR closely."
                    ),
                    recommendation_url="https://cpicpgx.org/guidelines/",
                )

        # ── Clopidogrel ──
        if drug_lower == "clopidogrel":
            if is_poor:
                return ClinicalRecommendation(
                    text=(
                        f"Use alternative antiplatelet therapy (e.g., prasugrel, ticagrelor). "
                        f"{gene} Poor Metabolizer status results in significantly reduced "
                        f"clopidogrel activation."
                    ),
                    implication=(
                        f"Reduced platelet inhibition due to decreased {gene}-mediated "
                        f"activation of clopidogrel. High risk of adverse cardiovascular events."
                    ),
                    recommendation_url="https://cpicpgx.org/guidelines/",
                )
            elif is_intermediate:
                return ClinicalRecommendation(
                    text=(
                        f"Consider alternative antiplatelet therapy or monitor closely. "
                        f"{gene} Intermediate Metabolizer status may reduce clopidogrel activation."
                    ),
                    implication=(
                        f"Moderately reduced platelet inhibition due to decreased {gene}-mediated "
                        f"activation of clopidogrel."
                    ),
                    recommendation_url="https://cpicpgx.org/guidelines/",
                )

        # ── Generic fallback ──
        phenotype_display = phenotype
        if is_poor:
            return ClinicalRecommendation(
                text=(
                    f"Exercise caution with {drug}. {gene} {phenotype_display} status may "
                    f"significantly alter drug metabolism. Consider dose adjustment or alternative."
                ),
                implication=f"Altered {gene} metabolism may affect {drug} response",
                recommendation_url=None,
            )
        elif is_intermediate:
            return ClinicalRecommendation(
                text=(
                    f"Monitor closely with {drug}. {gene} {phenotype_display} status may "
                    f"moderately alter drug metabolism. Consider dose adjustment."
                ),
                implication=f"Moderately altered {gene} metabolism may affect {drug} response",
                recommendation_url=None,
            )
        else:
            return ClinicalRecommendation(
                text=f"Exercise clinical judgment for {drug} dosing. {phenotype_display} phenotype.",
                implication=f"{phenotype_display} phenotype for {gene} with {drug}",
                recommendation_url=None,
            )


def create_risk_engine() -> RiskEngine:
    """Factory function to create a RiskEngine instance."""
    return RiskEngine()
