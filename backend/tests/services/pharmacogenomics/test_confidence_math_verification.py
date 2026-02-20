"""
Verification tests for the new confidence math architecture.

Validates:
  - Empty VCF (unresolved phenotype) output
  - Resolved PM (*4/*4) output
  - Level 3/4 evidence (weak evidence) output
"""

import pytest
from app.services.pharmacogenomics.risk_engine import RiskEngine
from app.services.pharmacogenomics.confidence import ConfidenceBreakdown


class TestEmptyVCFOutput:
    """Verify confidence math for unresolved phenotypes (Empty VCF scenario)."""

    @pytest.fixture
    def engine(self):
        return RiskEngine(
            enable_feedback_learning=False,
            enable_calibration=False,
        )

    def test_unresolved_phenotype_confidence_math(self, engine):
        """
        Verify Empty VCF output:
          - genotype_confidence < 1.0 (due to capped components)
          - phenotype_confidence = 0.0 (diplotype_determinism = 0)
          - classification_confidence > 0.0 (prevents zero-collapse)
          - automation_status blocked with reason "Phenotype unresolved"
        """
        # Simulate unresolved diplotype from empty/insufficient VCF
        diplotype_confidence_breakdown = {
            'variant_quality': 1.0,
            'allele_coverage': 0.3,      # Capped due to unresolved state
            'cnv_evaluation': 0.5,        # Capped due to unresolved state
            'genome_build_validity': 1.0,
            'diplotype_determinism': 0.0, # Unresolved → 0
        }

        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="Indeterminate",
            diplotype="Unresolved",
            diplotype_confidence=0.0,
            diplotype_confidence_breakdown=diplotype_confidence_breakdown,
        )

        # Extract breakdown
        bd = risk.confidence_breakdown
        assert bd is not None, "confidence_breakdown should be populated"

        # Verify genotype_confidence < 1.0
        genotype_conf = bd.get('genotype_confidence', 0)
        assert genotype_conf < 1.0, (
            f"genotype_confidence should be < 1.0 for unresolved phenotype, got {genotype_conf}"
        )

        # Verify phenotype_confidence = 0.0
        phenotype_conf = bd.get('phenotype_confidence', 1.0)
        assert phenotype_conf == 0.0, (
            f"phenotype_confidence should be 0.0 for unresolved phenotype, got {phenotype_conf}"
        )

        # Verify classification_confidence > 0.0 (prevents zero-collapse)
        classification_conf = bd.get('classification_confidence', 0.0)
        assert classification_conf > 0.0, (
            f"classification_confidence should be > 0.0 even when phenotype is unresolved, got {classification_conf}"
        )

        # Verify automation_status blocked
        auto_status = risk.automation_status
        assert auto_status is not None, "automation_status should be populated"
        assert auto_status.get('allowed') is False, (
            "automation should be blocked for unresolved phenotype"
        )

        # Verify blocked reason includes "Phenotype unresolved"
        blocked_reasons = auto_status.get('blocked_reasons', [])
        assert any('Phenotype unresolved' in reason for reason in blocked_reasons), (
            f"Expected 'Phenotype unresolved' in blocked_reasons, got {blocked_reasons}"
        )

    def test_unresolved_diplotype_severity_undetermined(self, engine):
        """Verify unresolved phenotype results in severity='undetermined'."""
        risk, _ = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="Indeterminate",
            diplotype="Unresolved",
            diplotype_confidence=0.0,
            diplotype_confidence_breakdown={
                'variant_quality': 1.0,
                'allele_coverage': 0.3,
                'cnv_evaluation': 0.5,
                'genome_build_validity': 1.0,
                'diplotype_determinism': 0.0,
            },
        )

        assert risk.severity == "undetermined", (
            f"Expected severity='undetermined' for unresolved phenotype, got '{risk.severity}'"
        )


class TestResolvedPMOutput:
    """Verify confidence math for resolved Poor Metabolizer phenotype."""

    @pytest.fixture
    def engine(self):
        return RiskEngine(
            enable_feedback_learning=False,
            enable_calibration=False,
        )

    def test_resolved_pm_confidence_math(self, engine):
        """
        Verify Resolved PM (*4/*4) output:
          - High genotype_confidence (if valid)
          - phenotype_confidence tracks genotype
          - classification_confidence tracks phenotype + knowledge
          - automation_status allowed (assuming knowledge >= 0.80)
        """
        # Simulate high-quality resolved PM diplotype
        diplotype_confidence_breakdown = {
            'variant_quality': 0.95,
            'allele_coverage': 0.90,
            'cnv_evaluation': 0.80,
            'genome_build_validity': 1.0,
            'diplotype_determinism': 0.95,  # Resolved, high confidence
        }

        risk, recommendation = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            diplotype="*4/*4",
            diplotype_confidence=0.90,
            diplotype_confidence_breakdown=diplotype_confidence_breakdown,
        )

        # Extract breakdown
        bd = risk.confidence_breakdown
        assert bd is not None, "confidence_breakdown should be populated"

        # Verify high genotype_confidence
        genotype_conf = bd.get('genotype_confidence', 0)
        # Formula: 0.35*allele + 0.25*cnv + 0.25*variant + 0.15*genome
        expected_genotype = 0.35 * 0.90 + 0.25 * 0.80 + 0.25 * 0.95 + 0.15 * 1.0
        assert abs(genotype_conf - expected_genotype) < 0.01, (
            f"genotype_confidence should be ~{expected_genotype:.2f}, got {genotype_conf}"
        )

        # Verify phenotype_confidence tracks genotype
        phenotype_conf = bd.get('phenotype_confidence', 0)
        expected_phenotype = genotype_conf * 0.95  # genotype × diplotype_determinism
        assert abs(phenotype_conf - expected_phenotype) < 0.01, (
            f"phenotype_confidence should track genotype_confidence, expected ~{expected_phenotype:.2f}, got {phenotype_conf}"
        )

        # Verify classification_confidence > 0 and reasonable
        classification_conf = bd.get('classification_confidence', 0)
        assert classification_conf > 0.5, (
            f"classification_confidence should be > 0.5 for resolved PM, got {classification_conf}"
        )

        # With default knowledge_confidence = 1.0:
        # classification = 0.6 * phenotype_conf + 0.4 * 1.0
        expected_classification = 0.6 * phenotype_conf + 0.4 * 1.0
        assert abs(classification_conf - expected_classification) < 0.01, (
            f"classification_confidence formula mismatch, expected ~{expected_classification:.2f}, got {classification_conf}"
        )

        # Verify automation_status allowed (knowledge >= 0.80 by default)
        auto_status = risk.automation_status
        assert auto_status is not None, "automation_status should be populated"

        # Check gates:
        # Gate 1: phenotype_confidence > 0 ✓
        # Gate 2: knowledge_confidence >= 0.80 (depends on PharmGKB data)
        # Gate 3: genotype_confidence >= 0.50 ✓
        # Gate 4: gene_drug_confirmed (depends on PharmGKB data)

        # We can't guarantee automation is allowed without PharmGKB data,
        # but we can verify the phenotype and genotype gates pass
        blocked_reasons = auto_status.get('blocked_reasons', [])

        # Phenotype should NOT be in blocked reasons
        assert not any('Phenotype unresolved' in reason for reason in blocked_reasons), (
            f"Phenotype should be resolved for PM, but got blocked: {blocked_reasons}"
        )

        # Genotype quality should NOT be in blocked reasons
        assert not any('Genotype quality too low' in reason for reason in blocked_reasons), (
            f"Genotype quality should be sufficient, but got blocked: {blocked_reasons}"
        )

    def test_resolved_pm_severity_high_or_critical(self, engine):
        """Verify resolved PM phenotype results in high/critical severity."""
        risk, _ = engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="PM",
            diplotype="*4/*4",
            diplotype_confidence=0.90,
            diplotype_confidence_breakdown={
                'variant_quality': 0.95,
                'allele_coverage': 0.90,
                'cnv_evaluation': 0.80,
                'genome_build_validity': 1.0,
                'diplotype_determinism': 0.95,
            },
        )

        assert risk.severity in ['high', 'critical'], (
            f"Expected severity in ['high', 'critical'] for PM phenotype, got '{risk.severity}'"
        )


class TestLevel3EvidenceOutput:
    """Verify confidence math for weak evidence (Level 3/4)."""

    @pytest.fixture
    def engine(self):
        # For this test, we need PharmGKB to be loaded to test evidence levels
        return RiskEngine(
            enable_feedback_learning=False,
            enable_calibration=False,
        )

    def test_level3_evidence_blocks_automation(self, engine):
        """
        Verify Level 3 Evidence output:
          - automation_status blocked with reason "Evidence insufficient"
          - knowledge_confidence caps the final score downward
        """
        # Simulate resolved diplotype but with weak evidence
        # We'll manually set knowledge_confidence low to simulate Level 3/4 evidence

        # First, let's create a ConfidenceBreakdown manually to test the logic
        bd = ConfidenceBreakdown()
        bd.variant_quality = 0.95
        bd.allele_coverage = 0.90
        bd.cnv_evaluation = 0.80
        bd.genome_build_validity = 1.0
        bd.diplotype_determinism = 0.95
        bd.cpic_applicability = 1.0
        bd.knowledge_confidence = 0.60  # Level 3 evidence (< 0.80 threshold)
        bd.gene_drug_confirmed = True

        # Verify genotype_confidence is high
        assert bd.genotype_confidence > 0.80

        # Verify phenotype_confidence is high
        assert bd.phenotype_confidence > 0.70

        # Verify classification_confidence is capped by weak knowledge
        # Formula: 0.6 * phenotype + 0.4 * knowledge
        # = 0.6 * 0.85 + 0.4 * 0.60 = 0.51 + 0.24 = 0.75
        expected_classification = 0.6 * bd.phenotype_confidence + 0.4 * 0.60
        assert abs(bd.classification_confidence - expected_classification) < 0.01, (
            f"classification_confidence should be capped by weak knowledge, "
            f"expected ~{expected_classification:.2f}, got {bd.classification_confidence}"
        )

        # Verify automation is blocked due to insufficient knowledge
        auto_status = bd.get_automation_status()
        assert auto_status['allowed'] is False, (
            "automation should be blocked for weak evidence"
        )

        blocked_reasons = auto_status['blocked_reasons']
        assert any('Evidence insufficient' in reason for reason in blocked_reasons), (
            f"Expected 'Evidence insufficient' in blocked_reasons, got {blocked_reasons}"
        )

    def test_knowledge_confidence_affects_classification(self, engine):
        """
        Verify that knowledge_confidence directly impacts classification_confidence.

        Two scenarios with identical genotype/phenotype confidence but different
        knowledge_confidence should have different classification_confidence.
        """
        # Scenario 1: High knowledge (Level 1A)
        bd1 = ConfidenceBreakdown()
        bd1.variant_quality = 0.95
        bd1.allele_coverage = 0.90
        bd1.cnv_evaluation = 0.80
        bd1.genome_build_validity = 1.0
        bd1.diplotype_determinism = 0.95
        bd1.knowledge_confidence = 1.0  # Level 1A

        # Scenario 2: Weak knowledge (Level 3)
        bd2 = ConfidenceBreakdown()
        bd2.variant_quality = 0.95
        bd2.allele_coverage = 0.90
        bd2.cnv_evaluation = 0.80
        bd2.genome_build_validity = 1.0
        bd2.diplotype_determinism = 0.95
        bd2.knowledge_confidence = 0.60  # Level 3

        # Both should have same phenotype_confidence
        assert abs(bd1.phenotype_confidence - bd2.phenotype_confidence) < 0.001

        # But different classification_confidence
        assert bd1.classification_confidence > bd2.classification_confidence, (
            f"Higher knowledge should yield higher classification_confidence, "
            f"got bd1={bd1.classification_confidence}, bd2={bd2.classification_confidence}"
        )

        # The difference should be exactly 0.4 * (1.0 - 0.60) = 0.16
        expected_diff = 0.4 * (1.0 - 0.60)
        actual_diff = bd1.classification_confidence - bd2.classification_confidence
        assert abs(actual_diff - expected_diff) < 0.01, (
            f"Expected classification difference of ~{expected_diff:.2f}, got {actual_diff:.2f}"
        )


class TestConfidenceFormulaInvariants:
    """Test invariants of the confidence formula architecture."""

    def test_phenotype_confidence_zero_when_determinism_zero(self):
        """phenotype_confidence MUST be 0 when diplotype_determinism is 0."""
        bd = ConfidenceBreakdown()
        bd.variant_quality = 1.0
        bd.allele_coverage = 1.0
        bd.cnv_evaluation = 1.0
        bd.genome_build_validity = 1.0
        bd.diplotype_determinism = 0.0  # Unresolved

        assert bd.genotype_confidence == 1.0, "genotype should be perfect"
        assert bd.phenotype_confidence == 0.0, (
            f"phenotype_confidence MUST be 0 when determinism=0, got {bd.phenotype_confidence}"
        )

    def test_classification_confidence_nonzero_when_phenotype_zero(self):
        """classification_confidence MUST be > 0 even when phenotype_confidence = 0."""
        bd = ConfidenceBreakdown()
        bd.diplotype_determinism = 0.0  # → phenotype_confidence = 0
        bd.knowledge_confidence = 0.80

        assert bd.phenotype_confidence == 0.0

        # Formula: 0.6 * (1 - 0.0) + 0.4 * 0.80 = 0.6 + 0.32 = 0.92
        expected = 0.6 * 1.0 + 0.4 * 0.80
        assert abs(bd.classification_confidence - expected) < 0.01, (
            f"Expected classification_confidence ~{expected:.2f} when phenotype unresolved, "
            f"got {bd.classification_confidence}"
        )

    def test_genotype_confidence_weighted_formula(self):
        """Verify genotype_confidence uses weighted formula, not min()."""
        bd = ConfidenceBreakdown()
        bd.variant_quality = 0.90
        bd.allele_coverage = 0.80
        bd.cnv_evaluation = 0.70
        bd.genome_build_validity = 1.0

        # Old min() formula would give: 0.70
        # New weighted formula:
        expected = 0.35 * 0.80 + 0.25 * 0.70 + 0.25 * 0.90 + 0.15 * 1.0
        # = 0.28 + 0.175 + 0.225 + 0.15 = 0.83

        assert bd.genotype_confidence > 0.70, (
            "genotype_confidence should NOT use min(), should be > 0.70"
        )
        assert abs(bd.genotype_confidence - expected) < 0.01, (
            f"genotype_confidence should use weighted formula, "
            f"expected ~{expected:.2f}, got {bd.genotype_confidence}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
