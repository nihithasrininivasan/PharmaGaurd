"""
Unit tests for risk engine.
Tests drug-specific risk assessment and clinical recommendations.
"""

import pytest
from app.services.pharmacogenomics.risk_engine import RiskEngine, create_risk_engine
from app.services.pharmacogenomics.models import (
    PatientProfile,
    DiplotypeResult,
    DrugAssessment
)


class TestRiskEngine:
    """Test RiskEngine risk assessment logic."""

    @pytest.fixture
    def engine(self):
        """Create a RiskEngine instance."""
        return create_risk_engine()

    # ===== Codeine - CYP2D6 Tests =====

    def test_codeine_poor_metabolizer(self, engine):
        """Test codeine risk for CYP2D6 Poor Metabolizer -> Avoid use"""
        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            diplotype="*4/*4",
            diplotype_confidence=0.95
        )

        assert "avoid" in risk.risk_label.lower() or "avoid" in recommendation.text.lower()
        assert risk.severity in ["high", "critical"]
        assert risk.confidence_score > 0.5
        assert "morphine" in recommendation.implication.lower()

    def test_codeine_intermediate_metabolizer(self, engine):
        """Test codeine risk for CYP2D6 Intermediate Metabolizer"""
        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="IM",
            diplotype="*1/*4",
            diplotype_confidence=0.90
        )

        assert risk.severity in ["moderate", "high"]
        assert risk.confidence_score > 0.5
        assert "alternative" in recommendation.text.lower() or "monitor" in recommendation.text.lower()

    def test_codeine_normal_metabolizer(self, engine):
        """Test codeine for Normal Metabolizer -> Standard dosing"""
        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=1.0
        )

        assert risk.severity == "none"
        assert "standard" in risk.risk_label.lower() or "normal" in recommendation.implication.lower()
        assert risk.confidence_score > 0.8

    def test_codeine_ultrarapid_metabolizer(self, engine):
        """Test codeine for Ultra-rapid Metabolizer -> Avoid (toxicity risk)"""
        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="UM",
            diplotype="*1/*2x2",  # Gene duplication
            diplotype_confidence=0.85
        )

        assert risk.severity in ["high", "critical"]
        assert "avoid" in risk.risk_label.lower() or "toxicity" in recommendation.implication.lower()

    # ===== Warfarin - CYP2C9 Tests =====

    def test_warfarin_poor_metabolizer(self, engine):
        """Test warfarin for CYP2C9 Poor Metabolizer -> Reduced dose"""
        risk, recommendation = engine.evaluate_risk(
            drug="warfarin",
            gene="CYP2C9",
            phenotype="PM",
            diplotype="*3/*3",
            diplotype_confidence=0.92
        )

        assert risk.severity in ["moderate", "high"]
        assert "dose" in risk.risk_label.lower() or "reduce" in risk.risk_label.lower()
        assert "bleeding" in recommendation.implication.lower()

    def test_warfarin_intermediate_metabolizer(self, engine):
        """Test warfarin for Intermediate Metabolizer"""
        risk, recommendation = engine.evaluate_risk(
            drug="warfarin",
            gene="CYP2C9",
            phenotype="IM",
            diplotype="*1/*3",
            diplotype_confidence=0.88
        )

        assert risk.severity in ["moderate", "high"]
        assert "dose" in risk.risk_label.lower() or "reduce" in recommendation.text.lower()

    def test_warfarin_normal_metabolizer(self, engine):
        """Test warfarin for Normal Metabolizer -> Standard dosing"""
        risk, recommendation = engine.evaluate_risk(
            drug="warfarin",
            gene="CYP2C9",
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=1.0
        )

        assert risk.severity == "none"
        assert "standard" in risk.risk_label.lower()

    # ===== Clopidogrel - CYP2C19 Tests =====

    def test_clopidogrel_poor_metabolizer(self, engine):
        """Test clopidogrel for CYP2C19 PM -> Use alternative"""
        risk, recommendation = engine.evaluate_risk(
            drug="clopidogrel",
            gene="CYP2C19",
            phenotype="PM",
            diplotype="*2/*2",
            diplotype_confidence=0.93
        )

        assert risk.severity in ["moderate", "high"]
        assert ("alternative" in risk.risk_label.lower() or
                "alternative" in recommendation.text.lower())
        assert "platelet" in recommendation.implication.lower()

    def test_clopidogrel_intermediate_metabolizer(self, engine):
        """Test clopidogrel for Intermediate Metabolizer"""
        risk, recommendation = engine.evaluate_risk(
            drug="clopidogrel",
            gene="CYP2C19",
            phenotype="IM",
            diplotype="*1/*2",
            diplotype_confidence=0.90
        )

        assert risk.severity in ["moderate", "high"]
        assert risk.confidence_score > 0.5

    def test_clopidogrel_normal_metabolizer(self, engine):
        """Test clopidogrel for Normal Metabolizer"""
        risk, recommendation = engine.evaluate_risk(
            drug="clopidogrel",
            gene="CYP2C19",
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=1.0
        )

        assert risk.severity == "none"
        assert "standard" in risk.risk_label.lower()

    # ===== Error Handling Tests =====

    def test_unsupported_drug(self, engine):
        """Test handling of unsupported drug"""
        risk, recommendation = engine.evaluate_risk(
            drug="unknown_drug_xyz",
            gene="CYP2D6",
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=1.0
        )

        assert risk.confidence_score == 0.0
        assert risk.severity == "none"
        assert "not in" in risk.risk_label.lower() or "not" in recommendation.text.lower()

    def test_gene_drug_mismatch(self, engine):
        """Test handling of incorrect gene for drug"""
        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2C19",  # Wrong gene (should be CYP2D6)
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=1.0
        )

        assert risk.confidence_score == 0.0
        assert "mismatch" in risk.risk_label.lower() or "incorrect" in recommendation.text.lower()


class TestPatientProfileEvaluation:
    """Test evaluation of drugs against patient profiles."""

    @pytest.fixture
    def engine(self):
        return create_risk_engine()

    @pytest.fixture
    def patient_profile(self):
        """Create a sample patient profile."""
        return PatientProfile(
            sample_id="PATIENT001",
            diplotypes={
                "CYP2D6": DiplotypeResult(
                    gene="CYP2D6",
                    diplotype="*1/*4",
                    phenotype="IM",
                    confidence=0.90,
                    is_indeterminate=False,
                    notes="Heterozygous"
                ),
                "CYP2C19": DiplotypeResult(
                    gene="CYP2C19",
                    diplotype="*1/*2",
                    phenotype="IM",
                    confidence=0.88,
                    is_indeterminate=False,
                    notes="Heterozygous"
                ),
                "CYP2C9": DiplotypeResult(
                    gene="CYP2C9",
                    diplotype="*1/*1",
                    phenotype="NM",
                    confidence=1.0,
                    is_indeterminate=False,
                    notes="Wildtype"
                )
            }
        )

    def test_evaluate_single_drug_for_patient(self, engine, patient_profile):
        """Test evaluating a single drug for a patient"""
        assessment = engine.evaluate_drug_for_patient("codeine", patient_profile)

        assert assessment is not None
        assert assessment.drug == "codeine"
        assert assessment.gene == "CYP2D6"
        assert assessment.phenotype == "IM"
        assert assessment.diplotype == "*1/*4"
        assert assessment.risk is not None
        assert assessment.recommendation is not None

    def test_evaluate_multiple_drugs_for_patient(self, engine, patient_profile):
        """Test evaluating multiple drugs for a patient"""
        drugs = ["codeine", "warfarin", "clopidogrel"]
        assessments = engine.evaluate_multiple_drugs(drugs, patient_profile)

        assert len(assessments) == 3
        assert all(isinstance(a, DrugAssessment) for a in assessments)

        drug_names = [a.drug for a in assessments]
        assert "codeine" in drug_names
        assert "warfarin" in drug_names
        assert "clopidogrel" in drug_names

    def test_evaluate_drug_not_in_profile(self, engine):
        """Test evaluating a drug when gene is not in patient profile"""
        limited_profile = PatientProfile(
            sample_id="PATIENT002",
            diplotypes={
                "CYP2D6": DiplotypeResult(
                    gene="CYP2D6",
                    diplotype="*1/*1",
                    phenotype="NM",
                    confidence=1.0,
                    is_indeterminate=False,
                    notes="Wildtype"
                )
            }
        )

        # Try to evaluate warfarin (needs CYP2C9, which is not in profile)
        assessment = engine.evaluate_drug_for_patient("warfarin", limited_profile)

        assert assessment is None  # Should return None when gene not in profile


class TestConfidenceScoring:
    """Test confidence score calculations."""

    @pytest.fixture
    def engine(self):
        return create_risk_engine()

    def test_confidence_formula(self, engine):
        """Test the confidence score calculation formula"""
        base = 0.9
        coverage = 0.95
        ambiguity = 1.0

        confidence = engine.calculate_confidence_score(base, coverage, ambiguity)

        expected = 0.9 * 0.95 * 1.0
        assert abs(confidence - expected) < 0.001

    def test_confidence_bounds(self, engine):
        """Test that confidence is always between 0 and 1"""
        # Test lower bound
        confidence = engine.calculate_confidence_score(0.0, 0.5, 0.5)
        assert 0.0 <= confidence <= 1.0

        # Test upper bound (should cap at 1.0)
        confidence = engine.calculate_confidence_score(1.0, 1.0, 1.0)
        assert confidence == 1.0

    def test_low_diplotype_confidence_affects_risk(self, engine):
        """Test that low diplotype confidence reduces risk assessment confidence"""
        # High confidence diplotype
        risk1, _ = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            diplotype="*4/*4",
            diplotype_confidence=0.95
        )

        # Low confidence diplotype
        risk2, _ = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            diplotype="*4/*4",
            diplotype_confidence=0.50
        )

        # Risk with higher diplotype confidence should have higher risk confidence
        assert risk1.confidence_score > risk2.confidence_score


class TestSeverityMapping:
    """Test severity level mapping from CPIC implications."""

    @pytest.fixture
    def engine(self):
        return create_risk_engine()

    def test_critical_severity(self, engine):
        """Test critical severity mapping"""
        severity = engine.map_severity_level("Contraindicated - avoid use due to life-threatening risk")
        assert severity == "critical"

    def test_high_severity(self, engine):
        """Test high severity mapping"""
        severity = engine.map_severity_level("High risk of therapeutic failure")
        assert severity == "high"

    def test_moderate_severity(self, engine):
        """Test moderate severity mapping"""
        severity = engine.map_severity_level("Consider dose adjustment or alternative therapy")
        assert severity == "moderate"

    def test_low_severity(self, engine):
        """Test low severity mapping"""
        severity = engine.map_severity_level("Informative - minor clinical impact")
        assert severity == "low"

    def test_none_severity(self, engine):
        """Test none severity mapping"""
        severity = engine.map_severity_level("Normal metabolism - use standard dosing")
        assert severity == "none"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
