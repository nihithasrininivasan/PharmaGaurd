#!/usr/bin/env python3
"""
Validation script for drug support detection fix.

Tests that warfarin (and other drugs with Level 1A/1B PharmGKB evidence)
are correctly detected as "supported" even without dedicated CPIC guideline files.
"""

import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent / "backend"
sys.path.insert(0, str(backend_path))

from app.services.pharmacogenomics.pharmgkb_loader import get_pharmgkb_loader


def test_warfarin_support_detection():
    """
    Test Case 1: WARFARIN Support Detection

    Problem: System incorrectly reports "Drug 'WARFARIN' not in CPIC database"
    Root Cause: Drug support checked against CPIC guideline files, not PharmGKB evidence
    Expected Fix: Check PharmGKB for Level 1A/1B evidence across all genes
    """
    print("=" * 80)
    print("TEST 1: WARFARIN Support Detection")
    print("=" * 80)

    try:
        pharmgkb = get_pharmgkb_loader()
    except Exception as e:
        print(f"✗ FAIL: Could not load PharmGKB data: {e}")
        return False

    # Test case-insensitive matching
    test_cases = [
        "warfarin",    # lowercase
        "WARFARIN",    # uppercase
        "Warfarin",    # title case
        " warfarin ",  # with whitespace
    ]

    print(f"\nInput: Drug name variations for warfarin")
    print(f"Expected: is_drug_supported() returns True for all variations")
    print()

    all_passed = True
    for drug_name in test_cases:
        is_supported = pharmgkb.is_drug_supported(drug_name)
        status = "✓ PASS" if is_supported else "✗ FAIL"
        print(f"{status}: '{drug_name}' → is_supported = {is_supported}")
        if not is_supported:
            all_passed = False

    # Explain why warfarin should be supported
    print("\n" + "-" * 80)
    print("RATIONALE:")
    print("Warfarin MUST be detected as supported because:")
    print("  1. PharmGKB contains Level 1A annotations for warfarin")
    print("  2. Multiple genes have 1A evidence: CYP2C9, VKORC1, CYP4F2")
    print("  3. CPIC has published warfarin-pharmacogenomic guidelines")
    print("  4. Detection should aggregate across ALL genes, not require single file")

    print("\n" + "-" * 80)
    if all_passed:
        print("✓✓✓ TEST 1 PASSED: WARFARIN correctly detected as supported ✓✓✓")
    else:
        print("✗✗✗ TEST 1 FAILED: WARFARIN incorrectly marked as unsupported ✗✗✗")

    return all_passed


def test_case_normalization():
    """
    Test Case 2: Case Normalization

    Verify that drug support detection is case-insensitive.
    """
    print("\n\n" + "=" * 80)
    print("TEST 2: Case Normalization for Common Drugs")
    print("=" * 80)

    try:
        pharmgkb = get_pharmgkb_loader()
    except Exception as e:
        print(f"✗ FAIL: Could not load PharmGKB data: {e}")
        return False

    # Test drugs that should have Level 1A/1B evidence
    known_supported_drugs = [
        ("clopidogrel", "CLOPIDOGREL"),
        ("codeine", "CODEINE"),
        ("simvastatin", "SIMVASTATIN"),
    ]

    print(f"\nTesting case-insensitive matching for known supported drugs:")
    print()

    all_passed = True
    for lowercase, uppercase in known_supported_drugs:
        lower_result = pharmgkb.is_drug_supported(lowercase)
        upper_result = pharmgkb.is_drug_supported(uppercase)

        if lower_result == upper_result and lower_result:
            print(f"✓ PASS: '{lowercase}' and '{uppercase}' both → {lower_result}")
        else:
            print(f"✗ FAIL: '{lowercase}' → {lower_result}, '{uppercase}' → {upper_result}")
            all_passed = False

    print("\n" + "-" * 80)
    if all_passed:
        print("✓✓✓ TEST 2 PASSED: Case normalization working correctly ✓✓✓")
    else:
        print("✗✗✗ TEST 2 FAILED: Inconsistent case handling ✗✗✗")

    return all_passed


def test_unsupported_drug():
    """
    Test Case 3: Truly Unsupported Drug

    Verify that drugs with no Level 1A/1B evidence are correctly marked as unsupported.
    """
    print("\n\n" + "=" * 80)
    print("TEST 3: Truly Unsupported Drug Detection")
    print("=" * 80)

    try:
        pharmgkb = get_pharmgkb_loader()
    except Exception as e:
        print(f"✗ FAIL: Could not load PharmGKB data: {e}")
        return False

    # Test a drug that should NOT have 1A/1B evidence
    # Using a made-up drug name
    unsupported_drugs = [
        "notarealdrug123",
        "fakemedicine",
        "testdrug999",
    ]

    print(f"\nTesting that truly unsupported drugs are correctly identified:")
    print()

    all_passed = True
    for drug_name in unsupported_drugs:
        is_supported = pharmgkb.is_drug_supported(drug_name)
        if not is_supported:
            print(f"✓ PASS: '{drug_name}' → is_supported = False (correct)")
        else:
            print(f"✗ FAIL: '{drug_name}' → is_supported = True (should be False)")
            all_passed = False

    print("\n" + "-" * 80)
    if all_passed:
        print("✓✓✓ TEST 3 PASSED: Unsupported drugs correctly detected ✓✓✓")
    else:
        print("✗✗✗ TEST 3 FAILED: False positives for unsupported drugs ✗✗✗")

    return all_passed


def test_aggregation_across_genes():
    """
    Test Case 4: Aggregation Across Multiple Genes

    Verify that a drug is marked as supported if it has 1A/1B evidence
    for ANY gene, not just a specific gene.
    """
    print("\n\n" + "=" * 80)
    print("TEST 4: Aggregation Across Genes")
    print("=" * 80)

    try:
        pharmgkb = get_pharmgkb_loader()
    except Exception as e:
        print(f"✗ FAIL: Could not load PharmGKB data: {e}")
        return False

    print(f"\nVerifying that warfarin support aggregates across multiple genes:")
    print()

    # Warfarin has evidence for multiple genes
    # Check that it's supported at the drug level
    is_supported = pharmgkb.is_drug_supported("warfarin")

    if is_supported:
        print(f"✓ PASS: warfarin is_supported = True")
        print(f"         (Aggregates evidence from CYP2C9, VKORC1, CYP4F2, etc.)")
    else:
        print(f"✗ FAIL: warfarin is_supported = False")
        print(f"         (Should aggregate evidence across all genes)")
        return False

    # Check evidence level for specific gene-drug pairs
    print(f"\nChecking evidence levels for warfarin-gene pairs:")

    genes_to_check = ["CYP2C9", "VKORC1", "CYP4F2"]
    evidence_found = []

    for gene in genes_to_check:
        evidence_info = pharmgkb.get_evidence_level(gene, "warfarin")
        level = evidence_info.get("level", "none")
        if level in ("1A", "1B"):
            print(f"  ✓ {gene} + warfarin: evidence_level = {level}")
            evidence_found.append(gene)
        else:
            print(f"    {gene} + warfarin: evidence_level = {level} (not 1A/1B)")

    print("\n" + "-" * 80)
    if evidence_found:
        print(f"✓✓✓ TEST 4 PASSED: Found 1A/1B evidence for {len(evidence_found)} gene(s) ✓✓✓")
        print(f"    Genes with 1A/1B evidence: {', '.join(evidence_found)}")
        return True
    else:
        print(f"✗✗✗ TEST 4 FAILED: No 1A/1B evidence found for warfarin ✗✗✗")
        return False


def test_decision_tree_logic():
    """
    Test Case 5: Decision Tree Logic

    Verify the complete decision tree:
    1. Normalize drug name (lowercase, strip)
    2. Check clinical variants for 1A/1B annotations
    3. Check gene-drug evidence index for 1A/1B
    4. Return True if ANY match found
    """
    print("\n\n" + "=" * 80)
    print("TEST 5: Decision Tree Logic Validation")
    print("=" * 80)

    try:
        pharmgkb = get_pharmgkb_loader()
    except Exception as e:
        print(f"✗ FAIL: Could not load PharmGKB data: {e}")
        return False

    print(f"\nDecision Tree:")
    print(f"  1. Normalize drug name → lowercase, strip whitespace")
    print(f"  2. Check clinical_variants for 1A/1B annotations")
    print(f"  3. Check gene_drug_evidence_index for 1A/1B")
    print(f"  4. Return True if ANY match found")
    print()

    # Test the decision tree with warfarin
    test_drug = "  WARFARIN  "  # uppercase with whitespace
    print(f"Input: '{test_drug}'")
    print(f"Step 1: Normalize → '{test_drug.strip().lower()}'")

    is_supported = pharmgkb.is_drug_supported(test_drug)
    print(f"Step 2-3: Check PharmGKB datasets...")
    print(f"Step 4: Result → is_supported = {is_supported}")

    print("\n" + "-" * 80)
    if is_supported:
        print(f"✓✓✓ TEST 5 PASSED: Decision tree executed correctly ✓✓✓")
        return True
    else:
        print(f"✗✗✗ TEST 5 FAILED: Decision tree logic error ✗✗✗")
        return False


def main():
    print("\n" + "=" * 80)
    print("DRUG SUPPORT DETECTION VALIDATION")
    print("=" * 80)
    print("\nValidating fix for:")
    print("  BUG: 'Drug WARFARIN not in CPIC database' (false negative)")
    print("  ROOT CAUSE: Checking CPIC guideline files instead of PharmGKB evidence")
    print("  FIX: Check PharmGKB for Level 1A/1B evidence across all genes")
    print()

    try:
        test1_pass = test_warfarin_support_detection()
        test2_pass = test_case_normalization()
        test3_pass = test_unsupported_drug()
        test4_pass = test_aggregation_across_genes()
        test5_pass = test_decision_tree_logic()

        print("\n\n" + "=" * 80)
        print("FINAL RESULTS")
        print("=" * 80)
        print(f"Test 1 (WARFARIN Support Detection):  {'PASSED ✓' if test1_pass else 'FAILED ✗'}")
        print(f"Test 2 (Case Normalization):          {'PASSED ✓' if test2_pass else 'FAILED ✗'}")
        print(f"Test 3 (Unsupported Drug Detection):  {'PASSED ✓' if test3_pass else 'FAILED ✗'}")
        print(f"Test 4 (Aggregation Across Genes):    {'PASSED ✓' if test4_pass else 'FAILED ✗'}")
        print(f"Test 5 (Decision Tree Logic):         {'PASSED ✓' if test5_pass else 'FAILED ✗'}")

        if test1_pass and test2_pass and test3_pass and test4_pass and test5_pass:
            print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
            print("\nDrug support detection fix successfully implemented:")
            print("  1. WARFARIN now correctly detected as supported")
            print("  2. Case-insensitive matching working")
            print("  3. Aggregation across genes functioning")
            print("  4. False positives prevented")
            return 0
        else:
            print("\n✗✗✗ SOME TESTS FAILED ✗✗✗")
            return 1

    except Exception as e:
        print(f"\n✗✗✗ ERROR: {e} ✗✗✗")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
