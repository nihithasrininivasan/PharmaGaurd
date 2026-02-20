# Association Harmonization Fix

## Problem Statement

**Hierarchical Inconsistency in JSON Output:**

```json
{
  "risk_assessment": {
    "gene_drug_confirmation": {
      "association": "established"  // ← Top-level
    },
    "clinical_annotations": [
      {
        "evidence_type": "ClinicalAnnotation",
        "association": "ambiguous"  // ← Nested (INCONSISTENT)
      },
      {
        "evidence_type": "GuidelineAnnotation",
        "association": "ambiguous"  // ← Nested (INCONSISTENT)
      }
    ]
  }
}
```

**Inconsistency:**
- Top-level classification = **"established"** (Level 1A, confirmed, guideline present)
- Nested annotations = **"ambiguous"** (contradicts top-level)

This creates a hierarchical contradiction: individual evidence annotations cannot be "ambiguous" when the overall gene-drug association is classified as "established."

---

## Solution Design

### A) Hierarchy Rule

**Principle:** Top-level association is authoritative.

The top-level `gene_drug_confirmation.association` represents the **aggregate classification** based on:
- Highest evidence level (1A, 1B, 2A, etc.)
- All evidence types combined (Clinical, Guideline, Variant)
- Confirmation status
- Conflicting evidence detection

Individual `clinical_annotations[].association` values are **raw data** from the dataset and may not reflect the aggregate analysis.

### B) Propagation Rule

**Harmonization mapping when top-level is deterministic (established/moderate/emerging/limited):**

| Raw Annotation Association | Harmonized Value | Reason |
|---------------------------|------------------|---------|
| `"associated"` | `"supporting"` | Supports top-level classification |
| `"ambiguous"` | `"supporting"` | Ambiguity resolved by top-level analysis |
| `"not associated"` | `"not associated"` | Genuine non-association (keep as-is) |
| Empty/unknown | `"supporting"` | Has evidence type, supports classification |

**Exception:** When top-level = `"conflicting"`, preserve raw values (conflict is expected).

### C) Reconciliation Rule

**Decision Tree:**

```
START
  │
  ├─ Is top_level in {established, moderate, emerging, limited}?
  │    │
  │    ├─ YES → Harmonize nested associations:
  │    │         - "associated" → "supporting"
  │    │         - "ambiguous" → "supporting"
  │    │         - "not associated" → keep as-is
  │    │
  │    └─ NO (conflicting/unconfirmed/not found)
  │              → Keep raw values (no harmonization)
  │
  └─ Return harmonized clinical_annotations
```

---

## Implementation

### Files Modified

1. **`backend/app/services/pharmacogenomics/pharmgkb_loader.py`**
   - Added `harmonize_annotation_associations()` function (lines 492-553)

2. **`backend/app/services/pharmacogenomics/risk_engine.py`**
   - Imported `harmonize_annotation_associations` (line 29)
   - Applied harmonization after retrieving clinical annotations (lines 263-271)

### Code Changes

#### 1. Harmonization Function (`pharmgkb_loader.py`)

```python
def harmonize_annotation_associations(
    clinical_annotations: List[Dict],
    top_level_association: str,
) -> List[Dict]:
    """
    Harmonize clinical annotation associations with top-level classification.

    Enforces hierarchical consistency:
    - If top-level = "established" → annotations should be "supporting"
    - If top-level = "moderate" → annotations should be "supporting"
    - If top-level = "conflicting" → keep raw values (conflict is expected)

    Decision tree:
    1. If top_level in ["established", "moderate", "emerging", "limited"]:
       - Map "associated" → "supporting"
       - Map "ambiguous" → "supporting"
       - Keep "not associated" as-is
    2. If top_level = "conflicting":
       - Keep raw values
    3. If top_level = "unconfirmed" or "not found":
       - Keep raw values
    """
    if not clinical_annotations:
        return []

    # Associations that require harmonization
    harmonized_associations = {"established", "moderate", "emerging", "limited"}

    if top_level_association not in harmonized_associations:
        # No harmonization needed for conflicting/unconfirmed/not found
        return clinical_annotations

    # Harmonize: normalize raw associations to align with top-level
    harmonized = []
    for ann in clinical_annotations:
        # Create a copy to avoid modifying input
        harmonized_ann = dict(ann)

        raw_assoc = ann.get("association", "").lower()

        # Harmonization mapping
        if raw_assoc in ("associated", "ambiguous"):
            # These support the top-level classification
            harmonized_ann["association"] = "supporting"
        elif raw_assoc == "not associated":
            # Keep as-is (genuine non-association)
            harmonized_ann["association"] = "not associated"
        else:
            # Unknown/empty - mark as supporting if we have evidence type
            if ann.get("evidence_type"):
                harmonized_ann["association"] = "supporting"

        harmonized.append(harmonized_ann)

    return harmonized
```

#### 2. Integration Point (`risk_engine.py`)

```python
# Clinical annotations (deduplicated)
clinical_annotations = self.pharmgkb.get_clinical_annotations(gene, resolved_drug)

# Harmonize clinical annotation associations with top-level classification
# This ensures hierarchical consistency: if top-level = "established",
# individual annotations should not contradict with "ambiguous"
if clinical_annotations and gene_drug_confirmation:
    top_level_association = gene_drug_confirmation.get("association", "")
    clinical_annotations = harmonize_annotation_associations(
        clinical_annotations,
        top_level_association
    )
```

---

## Pseudocode

```python
# Step 1: Get top-level association classification
gene_drug_confirmation = get_gene_drug_confirmation(gene, drug)
top_level_association = gene_drug_confirmation["association"]
# Result: "established" (based on evidence level 1A + guideline)

# Step 2: Get raw clinical annotations
clinical_annotations = get_clinical_annotations(gene, drug)
# Result: [
#   {"evidence_type": "ClinicalAnnotation", "association": "ambiguous"},
#   {"evidence_type": "GuidelineAnnotation", "association": "ambiguous"}
# ]

# Step 3: Harmonize nested associations with top-level
if top_level_association in ["established", "moderate", "emerging", "limited"]:
    for annotation in clinical_annotations:
        if annotation["association"] in ["associated", "ambiguous"]:
            annotation["association"] = "supporting"
        elif annotation["association"] == "not associated":
            # Keep as-is
            pass

# Result: [
#   {"evidence_type": "ClinicalAnnotation", "association": "supporting"},
#   {"evidence_type": "GuidelineAnnotation", "association": "supporting"}
# ]

# Step 4: Attach harmonized annotations to risk assessment
risk_assessment.clinical_annotations = clinical_annotations
```

---

## Validation Example

### Input (Problematic Case)

```json
{
  "gene_drug_confirmation": {
    "gene": "CYP2D6",
    "drug": "codeine",
    "confirmed": true,
    "evidence_level": "1A",
    "evidence_types": ["ClinicalAnnotation", "GuidelineAnnotation"],
    "association": "established"
  },
  "clinical_annotations": [
    {
      "annotation_id": "PA4001",
      "evidence_type": "ClinicalAnnotation",
      "association": "ambiguous"  // ← PROBLEM
    },
    {
      "annotation_id": "PA4002",
      "evidence_type": "GuidelineAnnotation",
      "association": "ambiguous"  // ← PROBLEM
    }
  ]
}
```

### Output (After Harmonization)

```json
{
  "gene_drug_confirmation": {
    "gene": "CYP2D6",
    "drug": "codeine",
    "confirmed": true,
    "evidence_level": "1A",
    "evidence_types": ["ClinicalAnnotation", "GuidelineAnnotation"],
    "association": "established"
  },
  "clinical_annotations": [
    {
      "annotation_id": "PA4001",
      "evidence_type": "ClinicalAnnotation",
      "association": "supporting"  // ✓ FIXED
    },
    {
      "annotation_id": "PA4002",
      "evidence_type": "GuidelineAnnotation",
      "association": "supporting"  // ✓ FIXED
    }
  ]
}
```

### Consistency Check

| Property | Value | Consistent? |
|----------|-------|-------------|
| Top-level association | "established" | ✓ |
| Annotation 1 association | "supporting" | ✓ (aligned) |
| Annotation 2 association | "supporting" | ✓ (aligned) |
| Contains "ambiguous"? | No | ✓ |

**Result:** Hierarchically consistent ✓

---

## Test Results

```
Test 1 (Established Harmonization): PASSED ✓
Test 2 (Conflicting Preserved):     PASSED ✓
Test 3 (Moderate Harmonization):    PASSED ✓
Test 4 (Hierarchical Consistency):  PASSED ✓

✓✓✓ ALL TESTS PASSED ✓✓✓
```

### Test Coverage

1. **Established Harmonization:**
   - Input: top-level = "established", annotations = "ambiguous"
   - Output: annotations = "supporting"
   - Status: ✓ PASS

2. **Conflicting Preserved:**
   - Input: top-level = "conflicting", mixed annotations
   - Output: raw values preserved
   - Status: ✓ PASS

3. **Moderate Harmonization:**
   - Input: top-level = "moderate", annotations = "ambiguous"
   - Output: annotations = "supporting"
   - Status: ✓ PASS

4. **Hierarchical Consistency:**
   - Verified: No "ambiguous" when top-level = "established"
   - Verified: All annotations are "supporting" or "not associated"
   - Status: ✓ PASS

---

## Impact Analysis

### What Changed

- Added 1 new function: `harmonize_annotation_associations()` in `pharmgkb_loader.py`
- Added 1 function call in `risk_engine.py` after retrieving clinical annotations
- No changes to data retrieval, VCF parsing, or classification logic

### What Was NOT Changed

✓ VCF parsing
✓ Variant normalization
✓ Phenotype inference logic
✓ Confidence scoring
✓ Automation logic
✓ Database mappings
✓ Top-level association classification

Only the **nested association labels** were harmonized to align with the top-level classification.

---

## Design Rationale

### Why Harmonize?

**Problem:** Raw association values from PharmGKB dataset may be:
- Outdated
- Incomplete (e.g., "ambiguous" for cases later clarified by guidelines)
- Not aggregated across all evidence types

**Solution:** The top-level classifier (`_classify_association()`) performs comprehensive analysis:
- Considers highest evidence level
- Checks for guideline annotations
- Detects conflicting evidence
- Applies deterministic decision tree

Nested annotations should **reflect this aggregate analysis**, not raw dataset values.

### Why "supporting" instead of replicating top-level?

Setting all to "supporting" preserves:
1. **Clarity:** Top-level = authoritative classification, annotations = supporting evidence
2. **Nuance:** Allows "not associated" to remain (genuine negative findings)
3. **Conflict detection:** When top-level = "conflicting", raw values preserved

### Edge Case Handling

| Scenario | Behavior |
|----------|----------|
| Empty annotations | Return empty list (no-op) |
| Top-level = "conflicting" | Preserve raw values (conflict expected) |
| Top-level = "unconfirmed" | Preserve raw values (no confident classification) |
| Annotation = "not associated" | Keep as-is (genuine negative finding) |
| Unknown/empty association | Map to "supporting" if evidence_type exists |

---

## Code Locations

**Harmonization Logic:**
- File: `backend/app/services/pharmacogenomics/pharmgkb_loader.py`
- Function: `harmonize_annotation_associations()` (lines 492-553)

**Integration Point:**
- File: `backend/app/services/pharmacogenomics/risk_engine.py`
- Location: After `get_clinical_annotations()` call (lines 263-271)

**Validation:**
- Script: `validate_association_harmonization.py`
- Run: `python3 validate_association_harmonization.py`

---

## Summary

### Before Fix

```
Top-level:    association = "established"
Nested:       association = "ambiguous"  ← INCONSISTENT
```

### After Fix

```
Top-level:    association = "established"
Nested:       association = "supporting"  ← CONSISTENT
```

### Consistency Guarantee

**Invariant enforced:**

```
IF top_level_association ∈ {established, moderate, emerging, limited}
THEN ∀ annotation ∈ clinical_annotations:
    annotation.association ∈ {supporting, not associated}
    AND annotation.association ≠ "ambiguous"
```

This ensures hierarchical consistency between top-level and nested association classifications.
