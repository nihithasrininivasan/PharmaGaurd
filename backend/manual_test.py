import sys
import os

# Add current directory to path to allow imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.services.pharmacogenomics.models import GenotypeData, VariantCall
from app.services.pharmacogenomics.phenotype_mapper import PhenotypeMapper
from app.services.pharmacogenomics.risk_engine import RiskEngine
from app.services.pharmacogenomics.cpic_loader import get_cpic_loader

def get_user_input(prompt):
    return input(prompt).strip()

def main():
    print("=== CPIC Pharmacogenomic Engine - Manual Test ===\n")

    # Initialize loader to check cache
    try:
        loader = get_cpic_loader()
    except Exception as e:
        print(f"Error initializing loader: {e}")
        print("Tip: Run 'python3 -m app.utils.cpic_etl' to generate cache first.")
        return

    # 1. Get Drug and Gene
    drug = get_user_input("Enter Drug Name (e.g., codeine): ") or "codeine"
    gene = get_user_input("Enter Gene Symbol (e.g., CYP2D6): ") or "CYP2D6"

    print(f"\n--- Testing for {drug} ({gene}) ---\n")

    # 2. Get Variants
    variants = []
    print("Enter observed variants (press Enter when done):")
    
    while True:
        print(f"\nVariant #{len(variants)+1}")
        pos_str = get_user_input("  Position (e.g., 42126611): ")
        if not pos_str:
            break
            
        ref = get_user_input("  Ref Allele (e.g., C): ").upper()
        alt = get_user_input("  Alt Allele (e.g., G): ").upper()
        zygosity = get_user_input("  Zygosity (HET/HOM_ALT): ").upper()
        
        if zygosity not in ["HET", "HOM_ALT"]:
            print("  Invalid zygosity. Must be HET or HOM_ALT.")
            continue
            
        try:
            pos = int(pos_str)
            variants.append(VariantCall(
                chrom="chr22", # Placeholder
                pos=pos,
                ref=ref,
                alt=alt,
                zygosity=zygosity,
                quality=100.0,
                filter="PASS"
            ))
            print("  Variant added.")
        except ValueError:
            print("  Invalid position number.")

    # 3. Process
    genotype_data = GenotypeData(
        sample_id="MANUAL_TEST",
        gene_symbol=gene,
        variants=variants,
        covered_positions=[], # Assume full coverage for test
        coverage_mean=30.0
    )

    print("\nProcessing Genotype...")
    mapper = PhenotypeMapper()
    diplotype_result = mapper.process_genotype(genotype_data)

    print("\n=== Diplotype Result ===")
    print(f"Diplotype: {diplotype_result.diplotype}")
    print(f"Phenotype: {diplotype_result.phenotype}")
    print(f"Confidence: {diplotype_result.confidence:.2f}")
    if diplotype_result.notes:
        print(f"Notes: {diplotype_result.notes}")

    print("\n=== Risk Assessment ===")
    engine = RiskEngine()
    risk, recommendation = engine.evaluate_risk(
        drug=drug,
        gene=gene,
        phenotype=diplotype_result.phenotype,
        diplotype=diplotype_result.diplotype,
        diplotype_confidence=diplotype_result.confidence
    )

    print(f"Risk: {risk.risk_label}")
    print(f"Severity: {risk.severity}")
    print(f"Confidence: {risk.confidence_score:.2f}")
    print(f"Recommendation: {recommendation.text}")
    print(f"Implication: {recommendation.implication}")

if __name__ == "__main__":
    main()
