"""
Integration tests for the complete pharmacogenomics pipeline.
Tests end-to-end workflow from genotype data to risk assessment.
"""

import pytest
from app.services.pharmacogenomics.models import (
    GenotypeData,
    VariantCall,
    PatientProfile
)
from app.services.pharmacogenomics.phenotype_mapper import PhenotypeMapper
from app.services.pharmacogenomics.risk_engine import RiskEngine


class TestEndToEndWorkflow:
    """Test complete workflow from genotype to risk assessment."""

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    @pytest.fixture
    def risk_engine(self):
        return RiskEngine()

    def test_complete_cyp2d6_workflow(self, mapper, risk_engine):
        """
        Test complete workflow:
        1. Start with genotype data (variants)
        2. Resolve diplotype and phenotype
        3. Generate risk assessment for codeine
        """
        # Step 1: Create genotype data for CYP2D6 *1/*4 (Intermediate Metabolizer)
        genotype = GenotypeData(
            sample_id="INTEGRATION_TEST_001",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    rsid="rs1135840",
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=55.0,
            covered_positions=[42126611]
        )

        # Step 2: Resolve diplotype and phenotype
        diplotype_result = mapper.process_genotype(genotype)

        assert diplotype_result.gene == "CYP2D6"
        assert diplotype_result.diplotype is not None
        assert diplotype_result.phenotype in ["PM", "IM", "NM", "RM", "UM", "Indeterminate"]
        assert 0.0 <= diplotype_result.confidence <= 1.0

        # Step 3: Generate risk assessment for codeine
        risk, recommendation = risk_engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype=diplotype_result.phenotype,
            diplotype=diplotype_result.diplotype,
            diplotype_confidence=diplotype_result.confidence
        )

        assert risk is not None
        assert recommendation is not None
        assert risk.severity in ["none", "low", "moderate", "high", "critical"]
        assert 0.0 <= risk.confidence_score <= 1.0
        assert len(recommendation.text) > 0

    def test_multi_gene_patient_profile(self, mapper, risk_engine):
        """
        Test processing multiple genes for a patient and evaluating multiple drugs.
        """
        # Create genotype data for multiple genes
        genotypes = [
            GenotypeData(
                sample_id="PATIENT_MULTI_001",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=50.0,
                covered_positions=[]
            ),
            GenotypeData(
                sample_id="PATIENT_MULTI_001",
                gene_symbol="CYP2C19",
                variants=[],
                coverage_mean=48.0,
                covered_positions=[]
            ),
            GenotypeData(
                sample_id="PATIENT_MULTI_001",
                gene_symbol="CYP2C9",
                variants=[],
                coverage_mean=52.0,
                covered_positions=[]
            )
        ]

        # Process all genotypes
        diplotype_results = mapper.process_multiple_genes(genotypes)

        assert len(diplotype_results) == 3
        assert "CYP2D6" in diplotype_results
        assert "CYP2C19" in diplotype_results
        assert "CYP2C9" in diplotype_results

        # Create patient profile
        patient_profile = PatientProfile(
            sample_id="PATIENT_MULTI_001",
            diplotypes=diplotype_results
        )

        # Evaluate multiple drugs
        drugs_to_test = ["codeine", "warfarin", "clopidogrel"]
        assessments = risk_engine.evaluate_multiple_drugs(drugs_to_test, patient_profile)

        assert len(assessments) == 3
        for assessment in assessments:
            assert assessment.drug in drugs_to_test
            assert assessment.risk is not None
            assert assessment.recommendation is not None

    def test_poor_metabolizer_high_risk_workflow(self, mapper, risk_engine):
        """
        Test workflow for a poor metabolizer with high-risk scenario.
        CYP2D6 *4/*4 (PM) + Codeine = High Risk
        """
        # Create genotype for CYP2D6 *4/*4 (homozygous variant)
        genotype = GenotypeData(
            sample_id="HIGH_RISK_001",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    rsid="rs1135840",
                    ref="C",
                    alt="G",
                    zygosity="HOM_ALT",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=60.0,
            covered_positions=[42126611]
        )

        # Resolve diplotype
        diplotype_result = mapper.process_genotype(genotype)

        # Should be poor metabolizer or similar
        assert diplotype_result.confidence > 0.5

        # Evaluate codeine risk
        risk, recommendation = risk_engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype=diplotype_result.phenotype,
            diplotype=diplotype_result.diplotype,
            diplotype_confidence=diplotype_result.confidence
        )

        # For PM phenotype with codeine, should be high severity
        if diplotype_result.phenotype == "PM":
            assert risk.severity in ["high", "critical"]
            assert "avoid" in recommendation.text.lower() or "alternative" in recommendation.text.lower()

    def test_normal_metabolizer_low_risk_workflow(self, mapper, risk_engine):
        """
        Test workflow for normal metabolizer (low risk scenario).
        CYP2D6 *1/*1 (NM) + Codeine = Low/No Risk
        """
        # Wildtype genotype
        genotype = GenotypeData(
            sample_id="LOW_RISK_001",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=55.0,
            covered_positions=list(range(42126500, 42127000))
        )

        # Resolve diplotype
        diplotype_result = mapper.process_genotype(genotype)

        assert diplotype_result.diplotype == "*1/*1"
        assert diplotype_result.phenotype == "NM"

        # Evaluate codeine risk
        risk, recommendation = risk_engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="NM",
            diplotype="*1/*1",
            diplotype_confidence=diplotype_result.confidence
        )

        # Should be low/no risk
        assert risk.severity in ["none", "low"]
        assert risk.confidence_score > 0.8

    def test_indeterminate_handling(self, mapper, risk_engine):
        """Test handling of indeterminate diplotype calls."""
        # Create genotype with unknown/ambiguous variants
        genotype = GenotypeData(
            sample_id="INDETERMINATE_001",
            gene_symbol="FAKE_GENE_XYZ",
            variants=[],
            coverage_mean=10.0,  # Low coverage
            covered_positions=[]
        )

        # Try to resolve
        diplotype_result = mapper.process_genotype(genotype)

        # Should handle gracefully
        assert diplotype_result is not None
        # Likely indeterminate or unknown
        assert diplotype_result.diplotype in ["Unknown", "Indeterminate"] or diplotype_result.is_indeterminate


class TestDataValidation:
    """Test data validation and error handling."""

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    def test_empty_variant_list(self, mapper):
        """Test handling of empty variant list"""
        genotype = GenotypeData(
            sample_id="EMPTY_001",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=50.0,
            covered_positions=[]
        )

        result = mapper.process_genotype(genotype)

        # Should assume wildtype
        assert result.diplotype == "*1/*1"
        assert result.phenotype == "NM"

    def test_variant_key_generation(self):
        """Test variant key generation"""
        variant = VariantCall(
            chrom="chr22",
            pos=42126611,
            rsid="rs1135840",
            ref="C",
            alt="G",
            zygosity="HET",
            quality=99.0,
            filter="PASS"
        )

        key = variant.variant_key()
        assert key == "42126611:C:G"

    def test_multiple_variants_same_position(self, mapper):
        """Test handling of multiple variants at the same position"""
        genotype = GenotypeData(
            sample_id="MULTI_VAR_001",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    rsid="rs1135840",
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                ),
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    rsid="rs1135840",
                    ref="C",
                    alt="T",  # Different alt allele
                    zygosity="HET",
                    quality=95.0,
                    filter="PASS"
                )
            ],
            coverage_mean=50.0,
            covered_positions=[42126611]
        )

        # Should handle without crashing
        result = mapper.process_genotype(genotype)
        assert result is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    @pytest.fixture
    def risk_engine(self):
        return RiskEngine()

    def test_very_low_quality_variant(self, mapper):
        """Test handling of low quality variants"""
        genotype = GenotypeData(
            sample_id="LOW_QUAL_001",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    rsid=None,
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=5.0,  # Very low quality
                    filter="LOW_QUAL"
                )
            ],
            coverage_mean=15.0,
            covered_positions=[42126611]
        )

        # Should still process but may have lower confidence
        result = mapper.process_genotype(genotype)
        assert result is not None

    def test_zero_confidence_diplotype(self, risk_engine):
        """Test risk assessment with zero confidence diplotype"""
        risk, recommendation = risk_engine.evaluate_risk(
            drug="codeine",
            gene="CYP2D6",
            phenotype="Indeterminate",
            diplotype="Unknown",
            diplotype_confidence=0.0
        )

        assert risk.confidence_score <= 0.5  # Should be low confidence
        assert risk is not None
        assert recommendation is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
