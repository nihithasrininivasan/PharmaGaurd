#!/usr/bin/env python3
"""
Validation script for association harmonization fix.

Tests that clinical_annotations[].association values are harmonized
with the top-level gene_drug_confirmation.association classification.
"""

import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent / "backend"
sys.path.insert(0, str(backend_path))

from app.services.pharmacogenomics.pharmgkb_loader import harmonize_annotation_associations


def test_harmonization_established():
    """
    Test Case 1: Top-level = "established"

    Problem: Top-level association = "established" but clinical annotations show "ambiguous"
    Expected: All annotations should be normalized to "supporting"
    """
    print("=" * 80)
    print("TEST 1: Harmonization with Top-Level = 'established'")
    print("=" * 80)

    # Input: clinical annotations with mixed raw associations
    clinical_annotations = [
        {
            "annotation_id": "PA1001",
            "gene": "CYP2D6",
            "drug": "codeine",
            "evidence_type": "ClinicalAnnotation",
            "association": "ambiguous",  # ← Should be harmonized
            "pmids": ["12345678"],
        },
        {
            "annotation_id": "PA1002",
            "gene": "CYP2D6",
            "drug": "codeine",
            "evidence_type": "GuidelineAnnotation",
            "association": "associated",  # ← Should be harmonized
            "pmids": ["23456789"],
        },
        {
            "annotation_id": "PA1003",
            "gene": "CYP2D6",
            "drug": "codeine",
            "evidence_type": "VariantAnnotation",
            "association": "not associated",  # ← Should be kept
            "pmids": ["34567890"],
        },
    ]

    top_level_association = "established"

    print(f"\nInput:")
    print(f"  top_level_association: {top_level_association}")
    print(f"  clinical_annotations:")
    for ann in clinical_annotations:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Harmonize
    harmonized = harmonize_annotation_associations(
        clinical_annotations,
        top_level_association
    )

    print(f"\nOutput:")
    print(f"  harmonized_clinical_annotations:")
    for ann in harmonized:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    all_passed = True

    # Check that "ambiguous" and "associated" were mapped to "supporting"
    for i, ann in enumerate(harmonized):
        original_assoc = clinical_annotations[i]["association"]
        new_assoc = ann["association"]

        if original_assoc == "ambiguous":
            if new_assoc == "supporting":
                print(f"✓ PASS: 'ambiguous' → 'supporting' for {ann['evidence_type']}")
            else:
                print(f"✗ FAIL: 'ambiguous' → '{new_assoc}' (expected 'supporting')")
                all_passed = False
        elif original_assoc == "associated":
            if new_assoc == "supporting":
                print(f"✓ PASS: 'associated' → 'supporting' for {ann['evidence_type']}")
            else:
                print(f"✗ FAIL: 'associated' → '{new_assoc}' (expected 'supporting')")
                all_passed = False
        elif original_assoc == "not associated":
            if new_assoc == "not associated":
                print(f"✓ PASS: 'not associated' kept as-is for {ann['evidence_type']}")
            else:
                print(f"✗ FAIL: 'not associated' → '{new_assoc}' (should be kept)")
                all_passed = False

    if all_passed:
        print("\n✓✓✓ TEST 1 PASSED: Associations harmonized correctly ✓✓✓")
    else:
        print("\n✗✗✗ TEST 1 FAILED: Harmonization errors ✗✗✗")

    return all_passed


def test_harmonization_conflicting():
    """
    Test Case 2: Top-level = "conflicting"

    When top-level = "conflicting", raw values should be preserved
    (because conflict is expected in this case)
    """
    print("\n\n" + "=" * 80)
    print("TEST 2: Harmonization with Top-Level = 'conflicting'")
    print("=" * 80)

    clinical_annotations = [
        {
            "annotation_id": "PA2001",
            "gene": "CYP2D6",
            "drug": "tamoxifen",
            "evidence_type": "ClinicalAnnotation",
            "association": "associated",
            "pmids": ["12345678"],
        },
        {
            "annotation_id": "PA2002",
            "gene": "CYP2D6",
            "drug": "tamoxifen",
            "evidence_type": "ClinicalAnnotation",
            "association": "not associated",
            "pmids": ["23456789"],
        },
    ]

    top_level_association = "conflicting"

    print(f"\nInput:")
    print(f"  top_level_association: {top_level_association}")
    print(f"  clinical_annotations:")
    for ann in clinical_annotations:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Harmonize
    harmonized = harmonize_annotation_associations(
        clinical_annotations,
        top_level_association
    )

    print(f"\nOutput:")
    print(f"  harmonized_clinical_annotations:")
    for ann in harmonized:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    all_passed = True

    # Check that raw values are preserved for conflicting
    for i, ann in enumerate(harmonized):
        original_assoc = clinical_annotations[i]["association"]
        new_assoc = ann["association"]

        if original_assoc == new_assoc:
            print(f"✓ PASS: '{original_assoc}' preserved for conflicting top-level")
        else:
            print(f"✗ FAIL: '{original_assoc}' → '{new_assoc}' (should be preserved)")
            all_passed = False

    if all_passed:
        print("\n✓✓✓ TEST 2 PASSED: Conflicting associations preserved ✓✓✓")
    else:
        print("\n✗✗✗ TEST 2 FAILED: Unexpected modification ✗✗✗")

    return all_passed


def test_harmonization_moderate():
    """
    Test Case 3: Top-level = "moderate"

    Similar to "established", should harmonize to "supporting"
    """
    print("\n\n" + "=" * 80)
    print("TEST 3: Harmonization with Top-Level = 'moderate'")
    print("=" * 80)

    clinical_annotations = [
        {
            "annotation_id": "PA3001",
            "gene": "CYP2C19",
            "drug": "clopidogrel",
            "evidence_type": "ClinicalAnnotation",
            "association": "ambiguous",
            "pmids": ["12345678"],
        },
    ]

    top_level_association = "moderate"

    print(f"\nInput:")
    print(f"  top_level_association: {top_level_association}")
    print(f"  clinical_annotations:")
    for ann in clinical_annotations:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Harmonize
    harmonized = harmonize_annotation_associations(
        clinical_annotations,
        top_level_association
    )

    print(f"\nOutput:")
    print(f"  harmonized_clinical_annotations:")
    for ann in harmonized:
        print(f"    - {ann['evidence_type']}: association = '{ann['association']}'")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    all_passed = True

    if harmonized[0]["association"] == "supporting":
        print(f"✓ PASS: 'ambiguous' → 'supporting' for moderate top-level")
    else:
        print(f"✗ FAIL: 'ambiguous' → '{harmonized[0]['association']}' (expected 'supporting')")
        all_passed = False

    if all_passed:
        print("\n✓✓✓ TEST 3 PASSED: Moderate associations harmonized ✓✓✓")
    else:
        print("\n✗✗✗ TEST 3 FAILED: Harmonization error ✗✗✗")

    return all_passed


def test_hierarchical_consistency():
    """
    Test Case 4: Hierarchical Consistency Check

    Verify that after harmonization:
    - No "ambiguous" exists when top-level = "established"
    - All associations are either "supporting" or "not associated"
    """
    print("\n\n" + "=" * 80)
    print("TEST 4: Hierarchical Consistency Validation")
    print("=" * 80)

    # The problematic case from the user's description
    clinical_annotations = [
        {
            "annotation_id": "PA4001",
            "gene": "CYP2D6",
            "drug": "codeine",
            "evidence_type": "ClinicalAnnotation",
            "association": "ambiguous",  # ← PROBLEM
            "pmids": ["12345678"],
        },
        {
            "annotation_id": "PA4002",
            "gene": "CYP2D6",
            "drug": "codeine",
            "evidence_type": "GuidelineAnnotation",
            "association": "ambiguous",  # ← PROBLEM
            "pmids": ["23456789"],
        },
    ]

    top_level_association = "established"

    print(f"\nProblem:")
    print(f"  Top-level: gene_drug_confirmation.association = '{top_level_association}'")
    print(f"  Nested:    clinical_annotations[].association = 'ambiguous'")
    print(f"  → Hierarchical inconsistency!")

    # Harmonize
    harmonized = harmonize_annotation_associations(
        clinical_annotations,
        top_level_association
    )

    print(f"\nAfter harmonization:")
    print(f"  Top-level: gene_drug_confirmation.association = '{top_level_association}'")
    print(f"  Nested:    clinical_annotations[].association:")
    for ann in harmonized:
        print(f"    - {ann['evidence_type']}: '{ann['association']}'")

    # Validation
    print("\n" + "-" * 80)
    print("VALIDATION:")

    all_passed = True

    # Check that NO "ambiguous" exists when top-level = "established"
    ambiguous_found = False
    for ann in harmonized:
        if ann["association"] == "ambiguous":
            ambiguous_found = True
            print(f"✗ FAIL: Found 'ambiguous' in {ann['evidence_type']} (inconsistent with 'established')")
            all_passed = False

    if not ambiguous_found:
        print(f"✓ PASS: No 'ambiguous' associations (consistent with top-level 'established')")

    # Check that all are "supporting"
    all_supporting = all(
        ann["association"] in ("supporting", "not associated")
        for ann in harmonized
    )
    if all_supporting:
        print(f"✓ PASS: All associations are 'supporting' or 'not associated'")
    else:
        print(f"✗ FAIL: Some associations are not 'supporting' or 'not associated'")
        all_passed = False

    if all_passed:
        print("\n✓✓✓ TEST 4 PASSED: Hierarchical consistency enforced ✓✓✓")
    else:
        print("\n✗✗✗ TEST 4 FAILED: Inconsistency remains ✗✗✗")

    return all_passed


def main():
    print("\n" + "=" * 80)
    print("ASSOCIATION HARMONIZATION VALIDATION")
    print("=" * 80)
    print("\nValidating hierarchical consistency between:")
    print("  - Top-level: gene_drug_confirmation.association")
    print("  - Nested:    clinical_annotations[].association")
    print()

    try:
        test1_pass = test_harmonization_established()
        test2_pass = test_harmonization_conflicting()
        test3_pass = test_harmonization_moderate()
        test4_pass = test_hierarchical_consistency()

        print("\n\n" + "=" * 80)
        print("FINAL RESULTS")
        print("=" * 80)
        print(f"Test 1 (Established Harmonization): {'PASSED ✓' if test1_pass else 'FAILED ✗'}")
        print(f"Test 2 (Conflicting Preserved):     {'PASSED ✓' if test2_pass else 'FAILED ✗'}")
        print(f"Test 3 (Moderate Harmonization):    {'PASSED ✓' if test3_pass else 'FAILED ✗'}")
        print(f"Test 4 (Hierarchical Consistency):  {'PASSED ✓' if test4_pass else 'FAILED ✗'}")

        if test1_pass and test2_pass and test3_pass and test4_pass:
            print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
            print("\nAssociation harmonization successfully implemented:")
            print("  1. Top-level 'established' → nested annotations = 'supporting'")
            print("  2. Conflicting top-level → raw values preserved")
            print("  3. Hierarchical consistency enforced")
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
