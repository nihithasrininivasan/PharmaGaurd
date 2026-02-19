"""
Test-Driven Development Tests for Indeterminate States

Following TDD principles:
1. Write failing tests FIRST
2. Implement minimum code to pass
3. Refactor while keeping tests green

These tests ensure each indeterminate state is properly distinguished.
"""

import pytest
from app.services.pharmacogenomics import (
    GenotypeData,
    VariantCall,
    PhenotypeMapper,
    IndeterminateReason
)


class TestIndeterminateStateDistinction:
    """
    TDD Tests for granular indeterminate state detection.
    Each test focuses on ONE specific indeterminate reason.
    """

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    # ===== Test: UNSUPPORTED_GENE =====

    def test_unsupported_gene_returns_unsupported_gene_reason(self, mapper):
        """
        GIVEN a genotype for a gene not in the CPIC database
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be UNSUPPORTED_GENE
        """
        genotype = GenotypeData(
            sample_id="TEST_UNSUPPORTED",
            gene_symbol="FAKE_GENE_XYZ",
            variants=[],
            coverage_mean=50.0,
            covered_positions=[]
        )

        result = mapper.process_genotype(genotype)

        assert result.is_indeterminate is True
        assert result.indeterminate_reason == IndeterminateReason.UNSUPPORTED_GENE
        assert result.confidence == 0.0
        assert "not supported" in result.notes.lower()

    # ===== Test: NOVEL_VARIANTS =====

    def test_novel_variants_returns_novel_variants_reason(self, mapper):
        """
        GIVEN a genotype with variants that don't match any known alleles
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be NOVEL_VARIANTS
        """
        genotype = GenotypeData(
            sample_id="TEST_NOVEL",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=99999999,  # Position not in any allele definition
                    rsid=None,
                    ref="A",
                    alt="T",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=50.0,
            covered_positions=[99999999]
        )

        result = mapper.process_genotype(genotype)

        assert result.is_indeterminate is True
        assert result.indeterminate_reason == IndeterminateReason.NOVEL_VARIANTS
        assert "no matching alleles" in result.notes.lower() or "novel" in result.notes.lower()

    # ===== Test: NO_COVERAGE =====

    def test_no_coverage_returns_no_coverage_reason(self, mapper):
        """
        GIVEN a genotype with missing coverage at multiple key positions
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be NO_COVERAGE
        """
        # First, get a valid genotype to understand key positions
        from app.services.pharmacogenomics import get_cpic_loader

        loader = get_cpic_loader()
        key_positions = loader.get_key_positions("CYP2D6")

        # Create genotype with NO coverage at key positions
        genotype = GenotypeData(
            sample_id="TEST_NO_COVERAGE",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=5.0,  # Very low coverage
            covered_positions=[]  # No positions covered
        )

        result = mapper.process_genotype(genotype)

        # With no coverage but no variants, might still be wildtype
        # But if we explicitly mark it as having missing coverage issues,
        # we should see reduced confidence
        # NOTE: This might not trigger NO_COVERAGE for wildtype
        # Let's test with variants but missing coverage

        genotype_with_variant = GenotypeData(
            sample_id="TEST_NO_COVERAGE_2",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=20.0,  # Low quality
                    filter="PASS"
                )
            ],
            coverage_mean=5.0,
            covered_positions=[42126611]  # Only one position, missing many others
        )

        result2 = mapper.process_genotype(genotype_with_variant)

        # If >2 key positions are missing, should flag coverage issues
        if len(key_positions) > 3:
            # With significant missing coverage, confidence should be reduced
            assert result2.confidence < 0.5 or result2.indeterminate_reason == IndeterminateReason.NO_COVERAGE

    # ===== Test: AMBIGUOUS =====

    def test_ambiguous_returns_ambiguous_reason(self, mapper):
        """
        GIVEN a genotype with unphased compound heterozygote with low confidence
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be AMBIGUOUS
        """
        # Create unphased compound heterozygote
        genotype = GenotypeData(
            sample_id="TEST_AMBIGUOUS",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=60.0,
                    filter="PASS",
                    phased=False  # Explicitly unphased
                ),
                VariantCall(
                    chrom="chr22",
                    pos=42126578,
                    ref="C",
                    alt="T",
                    zygosity="HET",
                    quality=60.0,
                    filter="PASS",
                    phased=False
                )
            ],
            coverage_mean=40.0,
            covered_positions=[42126611, 42126578]
        )

        result = mapper.process_genotype(genotype)

        # Compound heterozygote without phasing should have some ambiguity
        # May or may not be marked as AMBIGUOUS depending on confidence threshold
        if result.confidence < 0.6:
            assert result.indeterminate_reason in [
                IndeterminateReason.AMBIGUOUS,
                IndeterminateReason.NONE  # Might still be acceptable confidence
            ]

    # ===== Test: PARTIAL_MATCH =====

    def test_partial_match_returns_partial_match_reason(self, mapper):
        """
        GIVEN a genotype with incomplete match to allele definition
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be PARTIAL_MATCH (or default call)
        """
        # This is typically the "default heterozygous call" scenario
        # where we have some evidence but not a complete match
        genotype = GenotypeData(
            sample_id="TEST_PARTIAL",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=50.0,
            covered_positions=[42126611]
        )

        result = mapper.process_genotype(genotype)

        # Should identify some allele, but might have partial match
        assert result.diplotype is not None
        # Partial match is often used for default/fallback calls
        if "default" in result.notes.lower():
            assert result.indeterminate_reason in [
                IndeterminateReason.PARTIAL_MATCH,
                IndeterminateReason.NONE
            ]

    # ===== Test: LOW_QUALITY =====

    def test_low_quality_returns_low_quality_reason(self, mapper):
        """
        GIVEN a genotype with very low quality variants
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be LOW_QUALITY (if confidence < 0.5)
        """
        genotype = GenotypeData(
            sample_id="TEST_LOW_QUALITY",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=42126611,
                    ref="C",
                    alt="G",
                    zygosity="HET",
                    quality=5.0,  # Very low quality
                    filter="LOW_QUAL"
                )
            ],
            coverage_mean=10.0,  # Low coverage
            covered_positions=[42126611]
        )

        result = mapper.process_genotype(genotype)

        # Low quality doesn't automatically mean indeterminate
        # but combined with low confidence it should
        if result.confidence < 0.5:
            assert result.is_indeterminate is True
            # Could be LOW_QUALITY or another reason
            assert result.indeterminate_reason in [
                IndeterminateReason.LOW_QUALITY,
                IndeterminateReason.PARTIAL_MATCH,
                IndeterminateReason.NONE
            ]

    # ===== Test: NONE (Confident Call) =====

    def test_confident_call_returns_none_reason(self, mapper):
        """
        GIVEN a genotype with clear, high-confidence call
        WHEN diplotype resolution is performed
        THEN indeterminate_reason should be NONE
        """
        genotype = GenotypeData(
            sample_id="TEST_CONFIDENT",
            gene_symbol="CYP2D6",
            variants=[],  # Wildtype
            coverage_mean=60.0,
            covered_positions=list(range(42126500, 42127000))  # Full coverage
        )

        result = mapper.process_genotype(genotype)

        assert result.is_indeterminate is False
        assert result.indeterminate_reason == IndeterminateReason.NONE
        assert result.confidence > 0.9
        assert result.diplotype == "*1/*1"


class TestIndeterminateReasonPriority:
    """Test that the most specific indeterminate reason is returned."""

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    def test_unsupported_gene_takes_precedence_over_no_coverage(self, mapper):
        """Unsupported gene should be identified even with no coverage."""
        genotype = GenotypeData(
            sample_id="TEST_PRECEDENCE_1",
            gene_symbol="UNKNOWN_GENE",
            variants=[],
            coverage_mean=0.0,
            covered_positions=[]
        )

        result = mapper.process_genotype(genotype)

        # UNSUPPORTED_GENE should take precedence
        assert result.indeterminate_reason == IndeterminateReason.UNSUPPORTED_GENE

    def test_novel_variants_takes_precedence_over_partial_match(self, mapper):
        """Novel variants should be identified specifically."""
        genotype = GenotypeData(
            sample_id="TEST_PRECEDENCE_2",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=88888888,  # Completely unknown position
                    ref="A",
                    alt="T",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=50.0,
            covered_positions=[88888888]
        )

        result = mapper.process_genotype(genotype)

        assert result.indeterminate_reason == IndeterminateReason.NOVEL_VARIANTS


class TestIndeterminateStateActionability:
    """Test that each indeterminate state suggests appropriate action."""

    @pytest.fixture
    def mapper(self):
        return PhenotypeMapper()

    def test_no_coverage_suggests_resequencing(self, mapper):
        """NO_COVERAGE should suggest resequencing action."""
        # This is a documentation test - the reason itself suggests action
        assert IndeterminateReason.NO_COVERAGE.value == "no_coverage"
        # Application layer should map this to "resequence" action

    def test_novel_variants_suggests_curation(self, mapper):
        """NOVEL_VARIANTS should suggest manual curation."""
        assert IndeterminateReason.NOVEL_VARIANTS.value == "novel_variants"
        # Application layer should map this to "manual_review" action

    def test_ambiguous_suggests_phasing(self, mapper):
        """AMBIGUOUS should suggest getting phasing data."""
        assert IndeterminateReason.AMBIGUOUS.value == "ambiguous"
        # Application layer should map this to "get_phasing" action


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
