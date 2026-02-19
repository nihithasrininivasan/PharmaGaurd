#!/usr/bin/env python3
"""
Pharmacogenomics Service - Improvements Demo

Demonstrates all new features:
1. Configuration system
2. Granular indeterminate states
3. VCF phasing support
4. Activity scores from config
5. Auto-ETL capability
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.services.pharmacogenomics import (
    VariantCall,
    GenotypeData,
    PhenotypeMapper,
    RiskEngine,
    IndeterminateReason,
    get_config,
    update_config,
    get_confidence_penalties,
    get_diplotype_config,
    get_activity_scores
)


def demo_configuration_system():
    """Demo 1: Configuration System"""
    print("=" * 80)
    print("DEMO 1: Configuration System")
    print("=" * 80)

    # View current configuration
    config = get_config()
    penalties = get_confidence_penalties()
    diplotype_config = get_diplotype_config()

    print("\n📋 Current Configuration:")
    print(f"  Auto-run ETL: {config.auto_run_etl}")
    print(f"  Verbose Logging: {config.verbose_logging}")

    print("\n⚙️  Confidence Penalties:")
    print(f"  Missing Key Position: {penalties.missing_key_position} (reduces confidence by {int((1-penalties.missing_key_position)*100)}%)")
    print(f"  Unphased Heterozygote: {penalties.unphased_heterozygote} (reduces confidence by {int((1-penalties.unphased_heterozygote)*100)}%)")
    print(f"  Partial Match: {penalties.partial_allele_match} (reduces confidence by {int((1-penalties.partial_allele_match)*100)}%)")

    print("\n🧬 Diplotype Resolution Config:")
    print(f"  Homozygous Score Threshold: {diplotype_config.homozygous_score_threshold}")
    print(f"  Default Phase Assumption: {diplotype_config.default_phase_assumption}")

    # Demonstrate updating configuration
    print("\n🔧 Updating Configuration:")
    original_penalty = penalties.missing_key_position
    update_config(**{"confidence_penalties.missing_key_position": 0.75})

    updated_penalties = get_confidence_penalties()
    print(f"  Changed missing_key_position: {original_penalty} → {updated_penalties.missing_key_position}")

    # Restore original
    update_config(**{"confidence_penalties.missing_key_position": original_penalty})
    print(f"  Restored to: {original_penalty}")
    print()


def demo_granular_indeterminate_states():
    """Demo 2: Granular Indeterminate States"""
    print("=" * 80)
    print("DEMO 2: Granular Indeterminate States")
    print("=" * 80)

    mapper = PhenotypeMapper()

    scenarios = [
        {
            "name": "Unsupported Gene",
            "genotype": GenotypeData(
                sample_id="INDET_001",
                gene_symbol="UNKNOWN_GENE",
                variants=[],
                coverage_mean=50.0
            ),
            "expected": IndeterminateReason.UNSUPPORTED_GENE
        },
        {
            "name": "Novel Variants",
            "genotype": GenotypeData(
                sample_id="INDET_002",
                gene_symbol="CYP2D6",
                variants=[
                    VariantCall(
                        chrom="chr22", pos=99999999, ref="A", alt="T",
                        zygosity="HET", quality=50.0, filter="PASS"
                    )
                ],
                coverage_mean=40.0
            ),
            "expected": IndeterminateReason.NOVEL_VARIANTS
        },
        {
            "name": "No Coverage",
            "genotype": GenotypeData(
                sample_id="INDET_003",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=5.0,
                covered_positions=[]  # No coverage data
            ),
            "expected": None  # Will be wildtype but with low confidence
        }
    ]

    for scenario in scenarios:
        print(f"\n📊 Scenario: {scenario['name']}")
        result = mapper.process_genotype(scenario["genotype"])

        print(f"  Diplotype: {result.diplotype}")
        print(f"  Is Indeterminate: {result.is_indeterminate}")
        print(f"  Reason: {result.indeterminate_reason.value}")
        print(f"  Confidence: {result.confidence:.2%}")
        print(f"  Notes: {result.notes}")

        if scenario["expected"] and result.is_indeterminate:
            assert result.indeterminate_reason == scenario["expected"], \
                f"Expected {scenario['expected']}, got {result.indeterminate_reason}"
            print(f"  ✓ Correctly identified as {scenario['expected'].value}")
    print()


def demo_phasing_support():
    """Demo 3: VCF Phasing Support"""
    print("=" * 80)
    print("DEMO 3: VCF Phasing Support")
    print("=" * 80)

    mapper = PhenotypeMapper()

    # Scenario 1: Unphased compound heterozygote
    print("\n🔬 Scenario 1: Unphased Compound Heterozygote")
    unphased_genotype = GenotypeData(
        sample_id="PHASE_001",
        gene_symbol="CYP2D6",
        variants=[
            VariantCall(
                chrom="chr22", pos=42126611, ref="C", alt="G",
                zygosity="HET", quality=99.0, filter="PASS",
                phased=False
            ),
            VariantCall(
                chrom="chr22", pos=42126578, ref="C", alt="T",
                zygosity="HET", quality=95.0, filter="PASS",
                phased=False
            )
        ],
        coverage_mean=60.0,
        covered_positions=[42126611, 42126578]
    )

    unphased_result = mapper.process_genotype(unphased_genotype)
    print(f"  Diplotype: {unphased_result.diplotype}")
    print(f"  Phased: {unphased_result.phased}")
    print(f"  Confidence: {unphased_result.confidence:.2%}")
    print(f"  Notes: {unphased_result.notes}")

    # Scenario 2: Phased compound heterozygote
    print("\n🔬 Scenario 2: Phased Compound Heterozygote (Same Variants)")
    phased_genotype = GenotypeData(
        sample_id="PHASE_002",
        gene_symbol="CYP2D6",
        variants=[
            VariantCall(
                chrom="chr22", pos=42126611, ref="C", alt="G",
                zygosity="HET", quality=99.0, filter="PASS",
                phased=True, phase_set="12345", haplotype=1
            ),
            VariantCall(
                chrom="chr22", pos=42126578, ref="C", alt="T",
                zygosity="HET", quality=95.0, filter="PASS",
                phased=True, phase_set="12345", haplotype=0
            )
        ],
        coverage_mean=60.0,
        covered_positions=[42126611, 42126578]
    )

    phased_result = mapper.process_genotype(phased_genotype)
    print(f"  Diplotype: {phased_result.diplotype}")
    print(f"  Phased: {phased_result.phased}")
    print(f"  Confidence: {phased_result.confidence:.2%}")
    print(f"  Notes: {phased_result.notes}")

    # Compare
    print(f"\n📈 Impact of Phasing:")
    confidence_increase = phased_result.confidence - unphased_result.confidence
    print(f"  Confidence Increase: {confidence_increase:.2%}")
    print(f"  Phasing eliminates ambiguity penalty!")
    print()


def demo_activity_scores():
    """Demo 4: Activity Scores from Configuration"""
    print("=" * 80)
    print("DEMO 4: Activity Scores from Configuration")
    print("=" * 80)

    from app.services.pharmacogenomics import get_cpic_loader

    loader = get_cpic_loader()
    activity_config = get_activity_scores()

    print("\n📊 CYP2D6 Activity Scores:")

    alleles = ["*1", "*2", "*3", "*4", "*5", "*9", "*10", "*41"]
    for allele in alleles:
        score = loader.get_activity_score("CYP2D6", allele)
        function = "Normal" if score >= 1.0 else "Decreased" if score > 0 else "Non-functional"
        print(f"  {allele:5} → {score:.1f} ({function})")

    print("\n🎯 Phenotype Thresholds:")
    print(f"  Poor Metabolizer (PM):     0.0 - {activity_config.poor_metabolizer_max}")
    print(f"  Intermediate Metabolizer:  {activity_config.poor_metabolizer_max} - {activity_config.intermediate_metabolizer_max}")
    print(f"  Normal Metabolizer:        {activity_config.intermediate_metabolizer_max} - {activity_config.normal_metabolizer_max}")
    print(f"  Ultra-rapid Metabolizer:   > {activity_config.normal_metabolizer_max}")

    print("\n💡 Example Diplotype Scores:")
    diplotypes = ["*1/*1", "*1/*4", "*4/*4", "*1/*10"]
    for diplotype in diplotypes:
        total_score = loader.calculate_total_activity_score("CYP2D6", diplotype)
        alleles = diplotype.split("/")
        scores = [loader.get_activity_score("CYP2D6", a) for a in alleles]
        print(f"  {diplotype:8} = {scores[0]:.1f} + {scores[1]:.1f} = {total_score:.1f}")
    print()


def demo_auto_etl():
    """Demo 5: Auto-ETL Capability"""
    print("=" * 80)
    print("DEMO 5: Auto-ETL Capability")
    print("=" * 80)

    config = get_config()

    print("\n⚙️  Auto-ETL Configuration:")
    print(f"  Enabled: {config.auto_run_etl}")
    print(f"  CPIC Data Directory: {config.cpic_data_dir}")
    print(f"  Cache Path: {config.cpic_cache_path}")

    print("\n📝 How Auto-ETL Works:")
    print("  1. On first import, CPICDataLoader checks if cache exists")
    print("  2. If cache missing and auto_run_etl=True:")
    print("     → Looks for CPIC data files in cpic_data_dir")
    print("     → Automatically runs ETL process")
    print("     → Generates cpic_cache.json")
    print("     → Continues loading normally")
    print("  3. If cache exists:")
    print("     → Loads cache directly (fast)")

    print("\n✅ Benefits:")
    print("  • Zero-config setup for development")
    print("  • Automatic cache refresh when data updates")
    print("  • CI/CD friendly (no manual ETL step)")
    print("  • Clear error messages if data missing")
    print()


def main():
    """Run all demos"""
    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 15 + "PHARMACOGENOMICS SERVICE - IMPROVEMENTS DEMO" + " " * 19 + "║")
    print("╚" + "=" * 78 + "╝")
    print()

    try:
        demo_configuration_system()
        demo_granular_indeterminate_states()
        demo_phasing_support()
        demo_activity_scores()
        demo_auto_etl()

        print("=" * 80)
        print("✅ All demos completed successfully!")
        print("=" * 80)
        print("\n📚 See IMPROVEMENTS_IMPLEMENTED.md for full documentation")
        print()

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
