# Pharmacogenomics Risk Assessment Fix Summary

## Overview

Fixed two critical logical inconsistencies in the risk scoring and association classification logic:

1. **Confidence Score Inconsistency**: Overall confidence = 1.0 while phenotype_confidence = 0 and automation blocked
2. **Association Label Inconsistency**: Association = "ambiguous" with evidence_level = 1A, confirmed = true, and guideline evidence

---

## Fix 1: Confidence Scoring Model

### File Modified
`backend/app/services/pharmacogenomics/confidence.py`

### Location
Lines 137-165 (modified `final` property)

### Problem
The `classification_confidence` formula produced confidence = 1.0 when:
- phenotype_confidence = 0 (unresolved phenotype)
- automation_status.allowed = false (automation blocked)

This violated mathematical consistency: overall confidence cannot be maximal when phenotype determination failed and automation gates are blocked.

### Solution
Added deterministic caps to the `final` property:

```python
@property
def final(self) -> float:
    """
    Final confidence score with deterministic caps.

    Caps applied:
    - phenotype_confidence = 0 → cap at 0.50
    - automation_status.allowed = false → cap at 0.70
    """
    base_confidence = self.classification_confidence

    # Apply phenotype cap
    phenotype_cap = 0.50 if self.phenotype_confidence == 0 else 1.0

    # Apply automation cap
    automation_status = self.get_automation_status()
    automation_cap = 0.70 if not automation_status.get("allowed", True) else 1.0

    # Final confidence = min(base, phenotype_cap, automation_cap)
    final_confidence = min(base_confidence, phenotype_cap, automation_cap)

    return max(0.0, min(1.0, round(final_confidence, 4)))
```

### Validation
**Test Case:**
- phenotype_confidence = 0
- automation_status.allowed = false
- knowledge_confidence = 1.0

**Before Fix:**
- confidence_score = 1.0 ❌

**After Fix:**
- confidence_score = 0.50 ✓ (capped by phenotype_cap)

---

## Fix 2: Association Classification Model

### File Modified
`backend/app/services/pharmacogenomics/pharmgkb_loader.py`

### Locations
1. Lines 430-482: New `_classify_association()` helper function
2. Lines 185-253: Modified `confirm_gene_drug_pair()` method

### Problem
The association classification used simple keyword matching on raw association strings:
```python
# OLD (INCORRECT) LOGIC
if "associated" in associations:
    best_assoc = "associated"
elif "ambiguous" in associations:
    best_assoc = "ambiguous"  # ← Could be chosen even with 1A evidence!
else:
    best_assoc = "not associated"
```

This ignored:
- Evidence level quality (1A, 1B, etc.)
- Confirmation status
- Evidence types (guideline annotations)

### Solution
Implemented deterministic decision tree in `_classify_association()`:

```python
def _classify_association(
    confirmed: bool,
    evidence_level: str,
    evidence_types: List[str],
    raw_associations: set,
) -> str:
    """
    Deterministic association classification.

    Decision tree:
    1. If not confirmed → "unconfirmed"
    2. If conflicting evidence → "conflicting"
    3. If level 1A/1B + has guideline → "established"
    4. If level 2A/2B → "moderate"
    5. If level 3 and multiple evidence types → "emerging"
    6. Otherwise → "limited"
    """
    # Priority 1: Unconfirmed
    if not confirmed:
        return "unconfirmed"

    # Priority 2: Conflicting evidence
    has_associated = "associated" in raw_associations
    has_not_associated = "not associated" in raw_associations
    if has_associated and has_not_associated:
        return "conflicting"

    # Priority 3: High-certainty established association
    has_guideline = any("Guideline" in et for et in evidence_types)
    if evidence_level in ("1A", "1B") and has_guideline:
        return "established"  # ← FIX: Never "ambiguous" with 1A + guideline

    # Priority 4: Moderate evidence
    if evidence_level in ("2A", "2B"):
        return "moderate"

    # Priority 5: Emerging evidence
    if evidence_level == "3" and len(evidence_types) >= 3:
        return "emerging"

    # Default: Limited evidence
    return "limited"
```

Updated `confirm_gene_drug_pair()` to use this logic:
```python
# Get evidence level for this gene-drug pair
evidence_info = self.get_evidence_level(gene_key, drug_key)
evidence_level = evidence_info.get("level", "none")

# Classify association using deterministic rule set
association = _classify_association(
    confirmed=True,
    evidence_level=evidence_level,
    evidence_types=list(evidence_types),
    raw_associations=associations,
)
```

### Validation
**Test Case:**
- confirmed = true
- evidence_level = "1A"
- evidence_types = ["ClinicalAnnotation", "GuidelineAnnotation"]

**Before Fix:**
- association = "ambiguous" ❌

**After Fix:**
- association = "established" ✓

---

## Validation Results

All tests passed:

```
Test 1 (Confidence Scoring):      PASSED ✓
Test 2 (Association Classification): PASSED ✓
Test 3 (Edge Cases):                PASSED ✓

✓✓✓ ALL TESTS PASSED ✓✓✓
```

### Test Coverage
1. **Confidence Scoring**: Verified caps applied correctly
2. **Association Classification**: Verified 1A + guideline → "established"
3. **Edge Cases**: Tested conflicting, moderate, emerging, limited, unconfirmed

---

## Impact Summary

### What Changed
- Modified 1 property method (`final`) in `confidence.py`
- Added 1 new helper function (`_classify_association`) in `pharmgkb_loader.py`
- Modified 1 existing method (`confirm_gene_drug_pair`) in `pharmgkb_loader.py`

### What Was NOT Changed
- VCF parsing logic ✓
- Variant normalization ✓
- Quality metrics calculation ✓
- Phenotype inference logic ✓
- Database mappings ✓

Only the internal consistency of scoring and association classification was fixed.

---

## Mathematical Consistency Proof

### Before Fix
```
phenotype_confidence = 0
automation_status.allowed = false
→ confidence_score = 1.0  ← INCONSISTENT
```

### After Fix
```
phenotype_confidence = 0
automation_status.allowed = false
→ base_confidence = 1.0 (from classification_confidence)
→ phenotype_cap = 0.50
→ automation_cap = 0.70
→ confidence_score = min(1.0, 0.50, 0.70) = 0.50  ← CONSISTENT
```

### Semantic Consistency Proof

### Before Fix
```
evidence_level = "1A"
confirmed = true
has_guideline = true
→ association = "ambiguous"  ← SEMANTICALLY INCORRECT
```

### After Fix
```
evidence_level = "1A"
confirmed = true
has_guideline = true
→ Decision tree path: Priority 3
→ association = "established"  ← SEMANTICALLY CORRECT
```

---

## Files Modified

1. `backend/app/services/pharmacogenomics/confidence.py`
   - Modified `final` property to apply caps

2. `backend/app/services/pharmacogenomics/pharmgkb_loader.py`
   - Added `_classify_association()` helper function
   - Modified `confirm_gene_drug_pair()` method

3. `validate_fixes.py` (new)
   - Validation test suite

4. `FIX_SUMMARY.md` (this file)
   - Documentation

---

## References

### Code Locations

**Confidence Scoring:**
- File: `backend/app/services/pharmacogenomics/confidence.py`
- Method: `ConfidenceBreakdown.final` (property)
- Lines: 137-165

**Association Classification:**
- File: `backend/app/services/pharmacogenomics/pharmgkb_loader.py`
- Functions: `_classify_association()` (430-482), `confirm_gene_drug_pair()` (185-253)

### Validation
- Script: `validate_fixes.py`
- Run: `python3 validate_fixes.py`
