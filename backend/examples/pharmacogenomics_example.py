#!/usr/bin/env python3
"""
Pharmacogenomics Service Example

Demonstrates the complete workflow from variant data to risk assessment.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.services.pharmacogenomics import (
    VariantCall,
    GenotypeData,
    PhenotypeMapper,
    RiskEngine,
    PatientProfile
)


def example_1_simple_wildtype():
    """Example 1: Simple wildtype patient (no variants)"""
    print("=" * 80)
    print("Example 1: Wildtype Patient (Normal Metabolizer)")
    print("=" * 80)

    mapper = PhenotypeMapper()
    engine = RiskEngine()

    # Patient with no variants in CYP2D6 (wildtype)
    genotype = GenotypeData(
        sample_id="PATIENT_WT_001",
        gene_symbol="CYP2D6",
        variants=[],
        coverage_mean=55.0,
        covered_positions=list(range(42126500, 42127000))
    )

    # Resolve diplotype
    result = mapper.process_genotype(genotype)
    print(f"\nGenotype Analysis:")
    print(f"  Gene: {result.gene}")
    print(f"  Diplotype: {result.diplotype}")
    print(f"  Phenotype: {result.phenotype}")
    print(f"  Confidence: {result.confidence:.2%}")
    print(f"  Notes: {result.notes}")

    # Evaluate codeine risk
    risk, recommendation = engine.evaluate_risk(
        drug="codeine",
        gene="CYP2D6",
        phenotype=result.phenotype,
        diplotype=result.diplotype,
        diplotype_confidence=result.confidence
    )

    print(f"\nCodeine Risk Assessment:")
    print(f"  Risk Label: {risk.risk_label}")
    print(f"  Severity: {risk.severity}")
    print(f"  Confidence: {risk.confidence_score:.2%}")
    print(f"  Recommendation: {recommendation.text}")
    print(f"  Implication: {recommendation.implication}")
    print()


def example_2_poor_metabolizer():
    """Example 2: CYP2D6 Poor Metabolizer with high risk"""
    print("=" * 80)
    print("Example 2: Poor Metabolizer (High Risk)")
    print("=" * 80)

    mapper = PhenotypeMapper()
    engine = RiskEngine()

    # Patient with CYP2D6 *4 variant (homozygous)
    genotype = GenotypeData(
        sample_id="PATIENT_PM_001",
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
                filter="PASS",
                ad=[0, 45]
            )
        ],
        coverage_mean=60.0,
        covered_positions=[42126611, 42126605, 42126578]
    )

    # Resolve diplotype
    result = mapper.process_genotype(genotype)
    print(f"\nGenotype Analysis:")
    print(f"  Gene: {result.gene}")
    print(f"  Diplotype: {result.diplotype}")
    print(f"  Phenotype: {result.phenotype}")
    print(f"  Confidence: {result.confidence:.2%}")

    # Evaluate codeine risk
    risk, recommendation = engine.evaluate_risk(
        drug="codeine",
        gene="CYP2D6",
        phenotype=result.phenotype,
        diplotype=result.diplotype,
        diplotype_confidence=result.confidence
    )

    print(f"\nCodeine Risk Assessment:")
    print(f"  ⚠️  Risk Label: {risk.risk_label}")
    print(f"  ⚠️  Severity: {risk.severity.upper()}")
    print(f"  Confidence: {risk.confidence_score:.2%}")
    print(f"  Recommendation: {recommendation.text}")
    print(f"  Implication: {recommendation.implication}")
    if recommendation.recommendation_url:
        print(f"  Guidelines: {recommendation.recommendation_url}")
    print()


def example_3_compound_heterozygote():
    """Example 3: Compound heterozygote (two different variants)"""
    print("=" * 80)
    print("Example 3: Compound Heterozygote")
    print("=" * 80)

    mapper = PhenotypeMapper()
    engine = RiskEngine()

    # Patient with two different heterozygous variants
    genotype = GenotypeData(
        sample_id="PATIENT_CH_001",
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
        coverage_mean=58.0,
        covered_positions=[42126611, 42126578]
    )

    result = mapper.process_genotype(genotype)
    print(f"\nGenotype Analysis:")
    print(f"  Gene: {result.gene}")
    print(f"  Diplotype: {result.diplotype}")
    print(f"  Phenotype: {result.phenotype}")
    print(f"  Confidence: {result.confidence:.2%}")
    print(f"  Notes: {result.notes}")
    print()


def example_4_multi_gene_patient():
    """Example 4: Multi-gene patient profile with multiple drug assessments"""
    print("=" * 80)
    print("Example 4: Multi-Gene Patient Profile")
    print("=" * 80)

    mapper = PhenotypeMapper()
    engine = RiskEngine()

    # Create genotypes for multiple genes
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

    # Process all genes
    diplotype_results = mapper.process_multiple_genes(genotypes)

    print(f"\nPatient Pharmacogenomic Profile:")
    print(f"  Sample ID: PATIENT_MULTI_001")
    print(f"\n  Diplotypes:")
    for gene, result in diplotype_results.items():
        print(f"    {gene:12} {result.diplotype:10} → {result.phenotype:3} (confidence: {result.confidence:.2%})")

    # Create patient profile
    patient = PatientProfile(
        sample_id="PATIENT_MULTI_001",
        diplotypes=diplotype_results
    )

    # Evaluate multiple drugs
    drugs = ["codeine", "warfarin", "clopidogrel"]
    print(f"\n  Drug Risk Assessments:")

    assessments = engine.evaluate_multiple_drugs(drugs, patient)
    for assessment in assessments:
        severity_icon = {
            "critical": "🔴",
            "high": "🟠",
            "moderate": "🟡",
            "low": "🟢",
            "none": "✅"
        }.get(assessment.risk.severity, "❓")

        print(f"\n    {severity_icon} {assessment.drug.upper()} (via {assessment.gene})")
        print(f"       Phenotype: {assessment.phenotype}")
        print(f"       Risk: {assessment.risk.risk_label}")
        print(f"       Severity: {assessment.risk.severity}")
        print(f"       Recommendation: {assessment.recommendation.text}")
    print()


def example_5_confidence_scenarios():
    """Example 5: Different confidence scenarios"""
    print("=" * 80)
    print("Example 5: Confidence Scoring Scenarios")
    print("=" * 80)

    mapper = PhenotypeMapper()

    scenarios = [
        {
            "name": "High Confidence (Full Coverage)",
            "genotype": GenotypeData(
                sample_id="CONF_HIGH",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=60.0,
                covered_positions=list(range(42126500, 42127000))
            )
        },
        {
            "name": "Moderate Confidence (Partial Coverage)",
            "genotype": GenotypeData(
                sample_id="CONF_MED",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=30.0,
                covered_positions=[42126611]
            )
        },
        {
            "name": "Lower Confidence (Low Coverage)",
            "genotype": GenotypeData(
                sample_id="CONF_LOW",
                gene_symbol="CYP2D6",
                variants=[],
                coverage_mean=10.0,
                covered_positions=[]
            )
        }
    ]

    for scenario in scenarios:
        result = mapper.process_genotype(scenario["genotype"])
        print(f"\n  {scenario['name']}")
        print(f"    Diplotype: {result.diplotype}")
        print(f"    Confidence: {result.confidence:.2%}")
        print(f"    Coverage: {scenario['genotype'].coverage_mean}x")
    print()


def main():
    """Run all examples"""
    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 20 + "PHARMACOGENOMICS SERVICE EXAMPLES" + " " * 25 + "║")
    print("╚" + "=" * 78 + "╝")
    print()

    try:
        example_1_simple_wildtype()
        example_2_poor_metabolizer()
        example_3_compound_heterozygote()
        example_4_multi_gene_patient()
        example_5_confidence_scenarios()

        print("=" * 80)
        print("✅ All examples completed successfully!")
        print("=" * 80)
        print()

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
