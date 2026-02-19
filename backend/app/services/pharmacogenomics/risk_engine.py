"""
Risk Engine - Evaluates pharmacogenomic risk based on drug-gene-phenotype combinations.
Deterministic, rule-based system following CPIC guidelines.
"""

from typing import Tuple, Optional, Dict, List
from .models import (
    RiskAssessment,
    ClinicalRecommendation,
    DrugAssessment,
    DiplotypeResult,
    PatientProfile
)
from .cpic_loader import get_cpic_loader


class RiskEngine:
    """
    Evaluates pharmacogenomic risk for drug-gene-phenotype combinations.
    All logic is deterministic and based on CPIC guidelines.
    """

    def __init__(self):
        self.loader = get_cpic_loader()

    def evaluate_risk(
        self,
        drug: str,
        gene: str,
        phenotype: str,
        diplotype: str,
        diplotype_confidence: float = 1.0
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """
        Evaluate risk for a specific drug-gene-phenotype combination.

        Args:
            drug: Drug name (e.g., "codeine")
            gene: Gene symbol (e.g., "CYP2D6")
            phenotype: Phenotype (e.g., "PM", "IM", "NM", "RM", "UM")
            diplotype: Diplotype (e.g., "*1/*4")
            diplotype_confidence: Confidence score from diplotype calling (0-1)

        Returns:
            Tuple of (RiskAssessment, ClinicalRecommendation)
        """
        # Validate drug is supported
        if not self.loader.is_drug_supported(drug):
            return self._create_unsupported_drug_response(drug)

        # Validate gene matches drug
        expected_gene = self.loader.get_drug_gene(drug)
        if expected_gene and expected_gene != gene:
            return self._create_gene_mismatch_response(drug, gene, expected_gene)

        # Get recommendation for this drug-phenotype combination
        recommendation_data = self.loader.get_drug_recommendation_for_phenotype(drug, phenotype)

        if not recommendation_data:
            # No specific recommendation for this phenotype
            # This typically means normal metabolizer with no special considerations
            return self._create_standard_dosing_response(drug, phenotype)

        # Map recommendation to risk assessment
        risk_assessment = self._create_risk_assessment(
            recommendation_data, diplotype_confidence
        )

        clinical_recommendation = self._create_clinical_recommendation(
            recommendation_data
        )

        return risk_assessment, clinical_recommendation

    def evaluate_drug_for_patient(
        self,
        drug: str,
        patient_profile: PatientProfile
    ) -> Optional[DrugAssessment]:
        """
        Evaluate a drug against a patient's complete pharmacogenomic profile.

        Args:
            drug: Drug name
            patient_profile: Patient's complete diplotype profile

        Returns:
            DrugAssessment if gene is in profile, None otherwise
        """
        # Get gene associated with this drug
        gene = self.loader.get_drug_gene(drug)
        if not gene:
            return None

        # Get patient's diplotype for this gene
        diplotype_result = patient_profile.diplotypes.get(gene)
        if not diplotype_result:
            return None

        # Evaluate risk
        risk, recommendation = self.evaluate_risk(
            drug=drug,
            gene=gene,
            phenotype=diplotype_result.phenotype,
            diplotype=diplotype_result.diplotype,
            diplotype_confidence=diplotype_result.confidence
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
        """
        Evaluate multiple drugs for a patient.

        Args:
            drugs: List of drug names
            patient_profile: Patient's pharmacogenomic profile

        Returns:
            List of DrugAssessments
        """
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
        diplotype_confidence: float
    ) -> RiskAssessment:
        """Create RiskAssessment from CPIC recommendation data."""
        # Extract severity from recommendation data
        severity = recommendation_data.get('severity', 'none')

        # Get risk label
        risk_label = recommendation_data.get('risk', 'Unknown risk')

        # Calculate final confidence score
        # Start with diplotype confidence and adjust based on recommendation certainty
        confidence_score = diplotype_confidence

        # Reduce confidence for indeterminate or uncertain cases
        if 'indeterminate' in risk_label.lower() or 'unknown' in risk_label.lower():
            confidence_score *= 0.5

        return RiskAssessment(
            risk_label=risk_label,
            confidence_score=min(1.0, max(0.0, confidence_score)),
            severity=severity
        )

    def _create_clinical_recommendation(
        self,
        recommendation_data: Dict
    ) -> ClinicalRecommendation:
        """Create ClinicalRecommendation from CPIC recommendation data."""
        return ClinicalRecommendation(
            text=recommendation_data.get('risk', 'No specific recommendation'),
            implication=recommendation_data.get('implication', 'Standard considerations apply'),
            recommendation_url=recommendation_data.get('url')
        )

    def _create_unsupported_drug_response(
        self, drug: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """Create response for unsupported drug."""
        risk = RiskAssessment(
            risk_label=f"Drug '{drug}' not in CPIC database",
            confidence_score=0.0,
            severity="none"
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
        risk = RiskAssessment(
            risk_label=f"Gene mismatch: expected {expected_gene}, got {provided_gene}",
            confidence_score=0.0,
            severity="none"
        )

        recommendation = ClinicalRecommendation(
            text=f"Incorrect gene provided for {drug}. Expected {expected_gene}.",
            implication=f"Drug {drug} is primarily metabolized by {expected_gene}",
            recommendation_url=None
        )

        return risk, recommendation

    def _create_standard_dosing_response(
        self, drug: str, phenotype: str
    ) -> Tuple[RiskAssessment, ClinicalRecommendation]:
        """Create response for standard/normal metabolizer with no special considerations."""
        # Determine if this is truly a normal metabolizer
        is_normal = phenotype in ["NM", "Normal Metabolizer"]

        if is_normal:
            risk = RiskAssessment(
                risk_label="Standard dosing recommended",
                confidence_score=0.95,
                severity="none"
            )

            recommendation = ClinicalRecommendation(
                text=f"Use standard {drug} dosing guidelines",
                implication="Normal drug metabolism expected",
                recommendation_url=None
            )
        else:
            # Phenotype exists but no specific CPIC recommendation
            risk = RiskAssessment(
                risk_label="No specific CPIC recommendation for this phenotype",
                confidence_score=0.6,
                severity="low"
            )

            recommendation = ClinicalRecommendation(
                text=f"Consider clinical judgment for {drug} dosing in {phenotype} individuals",
                implication=f"Limited guidance for {phenotype} phenotype with {drug}",
                recommendation_url=None
            )

        return risk, recommendation

    def map_severity_level(self, cpic_implication: str) -> str:
        """
        Map CPIC implication text to severity level.

        Severity levels: "none", "low", "moderate", "high", "critical"
        """
        implication_lower = cpic_implication.lower()

        # Critical severity indicators
        if any(term in implication_lower for term in [
            'contraindicated', 'avoid', 'life-threatening', 'severe adverse',
            'black box', 'fatal'
        ]):
            return "critical"

        # High severity indicators
        if any(term in implication_lower for term in [
            'high risk', 'significantly increased', 'major toxicity',
            'therapeutic failure', 'ineffective'
        ]):
            return "high"

        # Moderate severity indicators
        if any(term in implication_lower for term in [
            'moderate', 'adjust dose', 'consider alternative',
            'increased risk', 'reduced efficacy', 'monitor closely'
        ]):
            return "moderate"

        # Low severity indicators
        if any(term in implication_lower for term in [
            'informative', 'minor', 'slight', 'possible'
        ]):
            return "low"

        # Default to none for normal/standard cases
        if any(term in implication_lower for term in [
            'normal', 'standard', 'no change', 'typical'
        ]):
            return "none"

        # Default to moderate if unclear
        return "moderate"

    def calculate_confidence_score(
        self,
        base_confidence: float,
        coverage_quality: float,
        ambiguity_penalty: float
    ) -> float:
        """
        Calculate final confidence score using the formula:
        Confidence = Base * Coverage * Ambiguity

        Args:
            base_confidence: Base confidence from diplotype calling (0-1)
            coverage_quality: Coverage quality factor (0-1)
            ambiguity_penalty: Penalty for ambiguous calls (0-1)

        Returns:
            Final confidence score (0-1)
        """
        confidence = base_confidence * coverage_quality * ambiguity_penalty
        return min(1.0, max(0.0, confidence))


def create_risk_engine() -> RiskEngine:
    """Factory function to create a RiskEngine instance."""
    return RiskEngine()
