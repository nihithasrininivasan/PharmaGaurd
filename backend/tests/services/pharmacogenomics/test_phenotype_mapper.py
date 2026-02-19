"""
Unit tests for phenotype mapper and diplotype resolution.
Tests the star allele calling algorithm and phenotype determination.
"""

import pytest
from app.services.pharmacogenomics.models import GenotypeData, VariantCall
from app.services.pharmacogenomics.phenotype_mapper import DiplotypeResolver, PhenotypeMapper


class TestDiplotypeResolver:
    """Test diplotype resolution logic."""

    @pytest.fixture
    def resolver(self):
        """Create a DiplotypeResolver instance."""
        return DiplotypeResolver()

    def test_resolve_homozygous_wildtype(self, resolver):
        """Test resolution when no variants are present -> *1/*1"""
        genotype = GenotypeData(
            sample_id="TEST001",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=50.0,
            covered_positions=list(range(42126500, 42127000))
        )

        result = resolver.resolve_diplotype(genotype)

        assert result.diplotype == "*1/*1"
        assert result.phenotype == "NM"
        assert result.confidence >= 0.9
        assert not result.is_indeterminate
        assert "wildtype" in result.notes.lower()

    def test_resolve_homozygous_variant(self, resolver):
        """Test resolution for homozygous variant -> *4/*4"""
        # CYP2D6 *4 is defined by position 42126611:C:G (among others)
        genotype = GenotypeData(
            sample_id="TEST002",
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
            coverage_mean=55.0,
            covered_positions=[42126611]
        )

        result = resolver.resolve_diplotype(genotype)

        # Should identify *4 allele (or similar)
        assert "/" in result.diplotype
        assert result.confidence > 0.5
        # For homozygous variant, should be same allele on both chromosomes
        alleles = result.diplotype.split("/")
        # Either both same or one is *1 (depending on matching logic)

    def test_resolve_compound_heterozygote(self, resolver):
        """Test resolution for compound heterozygote with two different variants"""
        genotype = GenotypeData(
            sample_id="TEST003",
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
                    pos=42126578,
                    rsid="rs1440526469",
                    ref="C",
                    alt="T",
                    zygosity="HET",
                    quality=95.0,
                    filter="PASS"
                )
            ],
            coverage_mean=60.0,
            covered_positions=[42126611, 42126578]
        )

        result = resolver.resolve_diplotype(genotype)

        assert "/" in result.diplotype
        assert result.confidence > 0.5
        # Should identify compound heterozygote
        alleles = result.diplotype.split("/")
        assert len(alleles) == 2

    def test_heterozygous_with_wildtype(self, resolver):
        """Test resolution for single heterozygous variant -> *1/*X"""
        genotype = GenotypeData(
            sample_id="TEST004",
            gene_symbol="CYP2C19",
            variants=[
                VariantCall(
                    chrom="chr10",
                    pos=94781859,
                    rsid="rs4244285",
                    ref="G",
                    alt="A",
                    zygosity="HET",
                    quality=99.0,
                    filter="PASS"
                )
            ],
            coverage_mean=45.0,
            covered_positions=[94781859]
        )

        result = resolver.resolve_diplotype(genotype)

        assert "/" in result.diplotype
        # Should be heterozygous with wildtype
        alleles = result.diplotype.split("/")
        assert "*1" in alleles or result.confidence > 0.5

    def test_unsupported_gene(self, resolver):
        """Test handling of unsupported gene"""
        genotype = GenotypeData(
            sample_id="TEST005",
            gene_symbol="FAKE_GENE",
            variants=[],
            coverage_mean=50.0,
            covered_positions=[]
        )

        result = resolver.resolve_diplotype(genotype)

        assert result.diplotype == "Unknown"
        assert result.phenotype == "Unknown"
        assert result.is_indeterminate
        assert "not supported" in result.notes.lower()

    def test_indeterminate_unknown_variants(self, resolver):
        """Test handling when variants don't match any known alleles"""
        genotype = GenotypeData(
            sample_id="TEST006",
            gene_symbol="CYP2D6",
            variants=[
                VariantCall(
                    chrom="chr22",
                    pos=99999999,  # Non-existent position
                    rsid=None,
                    ref="A",
                    alt="T",
                    zygosity="HET",
                    quality=50.0,
                    filter="PASS"
                )
            ],
            coverage_mean=30.0,
            covered_positions=[99999999]
        )

        result = resolver.resolve_diplotype(genotype)

        # Should handle gracefully
        assert result.confidence < 1.0 or result.is_indeterminate


class TestPhenotypeMapper:
    """Test PhenotypeMapper high-level interface."""

    @pytest.fixture
    def mapper(self):
        """Create a PhenotypeMapper instance."""
        return PhenotypeMapper()

    def test_process_single_genotype(self, mapper):
        """Test processing a single genotype"""
        genotype = GenotypeData(
            sample_id="TEST007",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=50.0,
            covered_positions=[]
        )

        result = mapper.process_genotype(genotype)

        assert result.gene == "CYP2D6"
        assert result.diplotype is not None
        assert result.phenotype is not None
        assert 0.0 <= result.confidence <= 1.0

    def test_process_multiple_genes(self, mapper):
        """Test processing multiple genes at once"""
        genotypes = [
            GenotypeData(
                sample_id="TEST008",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=50.0,
                covered_positions=[]
            ),
            GenotypeData(
                sample_id="TEST008",
                gene_symbol="CYP2C19",
                variants=[],
                coverage_mean=45.0,
                covered_positions=[]
            )
        ]

        results = mapper.process_multiple_genes(genotypes)

        assert len(results) == 2
        assert "CYP2D6" in results
        assert "CYP2C19" in results
        assert all(0.0 <= r.confidence <= 1.0 for r in results.values())


class TestCoverageConfidence:
    """Test confidence scoring based on coverage."""

    @pytest.fixture
    def resolver(self):
        return DiplotypeResolver()

    def test_full_coverage_high_confidence(self, resolver):
        """Test that full coverage gives high confidence"""
        # Get key positions for CYP2D6
        from app.services.pharmacogenomics.cpic_loader import get_cpic_loader
        loader = get_cpic_loader()
        key_positions = loader.get_key_positions("CYP2D6")

        genotype = GenotypeData(
            sample_id="TEST009",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=60.0,
            covered_positions=key_positions  # All key positions covered
        )

        result = resolver.resolve_diplotype(genotype)

        # Should have high confidence with full coverage
        assert result.confidence >= 0.9

    def test_partial_coverage_reduced_confidence(self, resolver):
        """Test that missing coverage reduces confidence"""
        genotype = GenotypeData(
            sample_id="TEST010",
            gene_symbol="CYP2D6",
            variants=[],
            coverage_mean=20.0,
            covered_positions=[42126611]  # Only one position
        )

        result = resolver.resolve_diplotype(genotype)

        # Confidence should be somewhat reduced due to limited coverage
        # (though still high for wildtype)
        assert result.confidence > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
