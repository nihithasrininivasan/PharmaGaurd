#!/usr/bin/env python3
"""
Validation script for confidence scoring and association classification fixes.

Tests the two logical inconsistencies that were fixed:
1. Confidence score = 1 with phenotype_confidence = 0 and automation blocked
2. Association = "ambiguous" with evidence_level = 1A and confirmed = true
"""

import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent / "backend"
sys.path.insert(0, str(backend_path))

from app.services.pharmacogenomics.confidence import ConfidenceBreakdown
from app.services.pharmacogenomics.pharmgkb_loader import _classify_association


def test_confidence_scoring_fix():
    """
    Test Case 1: Confidence Scoring Fix

    Problem: confidence_score = 1 while phenotype_confidence = 0 and automation blocked
    Expected: confidence_score should be capped at 0.50 (phenotype cap)
    """
    print("=" * 80)
    print("TEST 1: Confidence Scoring Fix")
    print("=" * 80)

    # Create confidence breakdown with problematic case
    bd = ConfidenceBreakdown()

    # High quality genotype components
    bd.variant_quality = 0.9
    bd.allele_coverage = 0.9
    bd.cnv_evaluation = 0.9
    bd.genome_build_validity = 1.0

    # BUT phenotype unresolved
    bd.diplotype_determinism = 0.0  # This makes phenotype_confidence = 0

    # High knowledge confidence
    bd.knowledge_confidence = 1.0

    # Gene-drug confirmed
    bd.gene_drug_confirmed = True

    print(f"\nInput:")
    print(f"  variant_quality: {bd.variant_quality}")
    print(f"  allele_coverage: {bd.allele_coverage}")
    print(f"  diplotype_determinism: {bd.diplotype_determinism}")
    print(f"  knowledge_confidence: {bd.knowledge_confidence}")

    print(f"\nDerived values:")
    print(f"  genotype_confidence: {bd.genotype_confidence}")
    print(f"  phenotype_confidence: {bd.phenotype_confidence}")
    print(f"  classification_confidence: {bd.classification_confidence}")

    automation_status = bd.get_automation_status()
    print(f"\nAutomation status:")
    print(f"  allowed: {automation_status['allowed']}")
    if not automation_status['allowed']:
        print(f"  blocked_reasons: {automation_status['blocked_reasons']}")

    final_confidence = bd.final
    print(f"\nFinal confidence score: {final_confidence}")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    if bd.phenotype_confidence == 0:
        print(f"✓ phenotype_confidence = {bd.phenotype_confidence}")
        if final_confidence <= 0.50:
            print(f"✓ PASS: confidence_score = {final_confidence} (capped at ≤0.50)")
            test1_pass = True
        else:
            print(f"✗ FAIL: confidence_score = {final_confidence} (should be ≤0.50)")
            test1_pass = False
    else:
        print(f"✗ Test setup error: phenotype_confidence should be 0, got {bd.phenotype_confidence}")
        test1_pass = False

    if not automation_status['allowed']:
        print(f"✓ automation_status.allowed = {automation_status['allowed']}")
        if final_confidence <= 0.70:
            print(f"✓ PASS: confidence_score = {final_confidence} (capped at ≤0.70)")
        else:
            print(f"✗ FAIL: confidence_score = {final_confidence} (should be ≤0.70)")
            test1_pass = False

    if test1_pass:
        print("\n✓✓✓ TEST 1 PASSED: Confidence scoring is mathematically consistent ✓✓✓")
    else:
        print("\n✗✗✗ TEST 1 FAILED: Confidence scoring inconsistency remains ✗✗✗")

    return test1_pass


def test_association_classification_fix():
    """
    Test Case 2: Association Classification Fix

    Problem: association = "ambiguous" with evidence_level = 1A, confirmed = true, and guideline
    Expected: association should be "established"
    """
    print("\n\n" + "=" * 80)
    print("TEST 2: Association Classification Fix")
    print("=" * 80)

    # Test case from the problem description
    confirmed = True
    evidence_level = "1A"
    evidence_types = ["ClinicalAnnotation", "GuidelineAnnotation"]
    raw_associations = {"associated"}  # Raw data says associated

    print(f"\nInput:")
    print(f"  confirmed: {confirmed}")
    print(f"  evidence_level: {evidence_level}")
    print(f"  evidence_types: {evidence_types}")
    print(f"  raw_associations: {raw_associations}")

    association = _classify_association(
        confirmed=confirmed,
        evidence_level=evidence_level,
        evidence_types=evidence_types,
        raw_associations=raw_associations,
    )

    print(f"\nClassified association: {association}")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    has_guideline = any("Guideline" in et for et in evidence_types)
    print(f"✓ confirmed = {confirmed}")
    print(f"✓ evidence_level = {evidence_level} (high quality)")
    print(f"✓ has_guideline = {has_guideline}")

    if association == "established":
        print(f"✓ PASS: association = '{association}' (correct for 1A + guideline)")
        test2_pass = True
    else:
        print(f"✗ FAIL: association = '{association}' (should be 'established')")
        test2_pass = False

    if test2_pass:
        print("\n✓✓✓ TEST 2 PASSED: Association classification is deterministic ✓✓✓")
    else:
        print("\n✗✗✗ TEST 2 FAILED: Association classification inconsistency remains ✗✗✗")

    return test2_pass


def test_edge_cases():
    """Test edge cases for association classification."""
    print("\n\n" + "=" * 80)
    print("TEST 3: Association Classification Edge Cases")
    print("=" * 80)

    test_cases = [
        {
            "name": "Conflicting evidence",
            "confirmed": True,
            "evidence_level": "1A",
            "evidence_types": ["ClinicalAnnotation"],
            "raw_associations": {"associated", "not associated"},
            "expected": "conflicting",
        },
        {
            "name": "Moderate evidence (2A)",
            "confirmed": True,
            "evidence_level": "2A",
            "evidence_types": ["ClinicalAnnotation"],
            "raw_associations": {"associated"},
            "expected": "moderate",
        },
        {
            "name": "Emerging evidence (level 3 with multiple sources)",
            "confirmed": True,
            "evidence_level": "3",
            "evidence_types": ["ClinicalAnnotation", "VariantAnnotation", "OtherAnnotation"],
            "raw_associations": {"associated"},
            "expected": "emerging",
        },
        {
            "name": "Limited evidence (level 3 with few sources)",
            "confirmed": True,
            "evidence_level": "3",
            "evidence_types": ["ClinicalAnnotation"],
            "raw_associations": {"associated"},
            "expected": "limited",
        },
        {
            "name": "Unconfirmed",
            "confirmed": False,
            "evidence_level": "1A",
            "evidence_types": ["GuidelineAnnotation"],
            "raw_associations": set(),
            "expected": "unconfirmed",
        },
    ]

    all_passed = True
    for tc in test_cases:
        result = _classify_association(
            confirmed=tc["confirmed"],
            evidence_level=tc["evidence_level"],
            evidence_types=tc["evidence_types"],
            raw_associations=tc["raw_associations"],
        )

        passed = result == tc["expected"]
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"\n{status}: {tc['name']}")
        print(f"  Expected: {tc['expected']}, Got: {result}")

        if not passed:
            all_passed = False

    print("\n" + "-" * 80)
    if all_passed:
        print("✓✓✓ TEST 3 PASSED: All edge cases handled correctly ✓✓✓")
    else:
        print("✗✗✗ TEST 3 FAILED: Some edge cases failed ✗✗✗")

    return all_passed


def main():
    print("\n" + "=" * 80)
    print("PHARMACOGENOMICS PIPELINE FIX VALIDATION")
    print("=" * 80)
    print("\nValidating fixes for:")
    print("1. Confidence scoring consistency")
    print("2. Association classification logic")
    print()

    try:
        test1_pass = test_confidence_scoring_fix()
        test2_pass = test_association_classification_fix()
        test3_pass = test_edge_cases()

        print("\n\n" + "=" * 80)
        print("FINAL RESULTS")
        print("=" * 80)
        print(f"Test 1 (Confidence Scoring):      {'PASSED ✓' if test1_pass else 'FAILED ✗'}")
        print(f"Test 2 (Association Classification): {'PASSED ✓' if test2_pass else 'FAILED ✗'}")
        print(f"Test 3 (Edge Cases):                {'PASSED ✓' if test3_pass else 'FAILED ✗'}")

        if test1_pass and test2_pass and test3_pass:
            print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
            print("\nBoth logical inconsistencies have been successfully fixed:")
            print("  1. Confidence score now caps at 0.50 when phenotype_confidence = 0")
            print("  2. Association is now 'established' for 1A/1B + guideline evidence")
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
