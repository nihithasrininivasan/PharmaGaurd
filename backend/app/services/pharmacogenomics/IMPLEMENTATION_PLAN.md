# CPIC-Aligned Pharmacogenomic Decision Engine Implementation Plan

## Goal Description
Implement a deterministic, unit-testable backend medical logic layer for Pharmacogenomic risk assessment based on CPIC guidelines. This system will ingest VCF data, resolve diplotypes/phenotypes, and generate risk assessments for specific drugs.

## User Review Required
> [!IMPORTANT]
> The Output Schema is STRICT and cannot be changed. All internal logic must map to this final JSON structure.

> [!WARNING]
> No LLMs are to be used for risk classification. All logic must be rule-based and deterministic.

## Proposed Architecture & Logic

### 1. Internal Input Data Model (Intermediate VCF Output)
This schema represents the localized, parsed data from the VCF, ready for the medical logic layer. It captures all necessary information for accurate calling.

```python
from pydantic import BaseModel
from typing import List, Optional, Dict

class VariantCall(BaseModel):
    chrom: str
    pos: int
    rsid: Optional[str]
    ref: str
    alt: str
    zygosity: str  # "HET", "HOM_REF", "HOM_ALT"
    quality: float
    filter: str
    # Allele depth or frequency if available
    ad: Optional[List[int]] = None 

class GenotypeData(BaseModel):
    sample_id: str
    gene_symbol: str
    variants: List[VariantCall]
    # Coverage metadata (e.g., mean depth across gene region)
    coverage_mean: Optional[float] 
    # List of specific positions covered (for confidence calculation)
    covered_positions: List[int]
```

### 2. CPIC Dataset Processing Strategy
**Source:** Local CPIC Excel files in `backend/data/cpic/{GENE}/*allele_definition_table.xlsx` and `*diplotype_phenotype_table.xlsx`.
**Transformation:**
- **ETL Script (`utils/cpic_etl.py`)**: Run periodically (or at build time) to convert raw Excel data into optimized JSON lookups.
   - **Input:**
     - `CYP2D6_allele_definition_table.xlsx`: Contains nucleotide definitions for each star allele.
     - `CYP2D6_diplotype_phenotype.xlsx`: (If present) Maps diplotypes to phenotypes.
     - `CYP2D6_CDS.xlsx` (Clinical Decision Support): Contains drug recommendations.
   - **Logic:**
     - Read `Alleles` sheet. Parse header rows to map column index to genomic position `g.XXXX`.
     - Read rows for each allele (e.g., `*1`, `*2`).
     - Handle `NaN`: Treat as reference (`*1`) match.
     - Handle Indels/CNVs if represented in specific notation.
   - **Output:** `backend/data/cpic_cache.json` (keyed by Gene).

- **Runtime Structures (`cpic_loader.py`)**:
  - Load `cpic_cache.json` at startup (singleton).
  - **Lookup 1: variant_to_allele**: `{ "GENE": { "POS:REF:ALT": ["*2", "*3"] } }` (Inverse map).
  - **Lookup 2: allele_definitions**: `{ "GENE": { "*2": ["POS:REF:ALT", ...] } }`.
  - **Lookup 3: phenotype_map**: `{ "GENE": { "*1/*1": "NM", "*1/*2": "IM", ... } }`.
  - **Lookup 4: drug_recommendations**: `{ "DRUG": { "PHENOTYPE": { "risk": "...", "severity": "...", "text": "..." } } }`.

### 3. Diplotype Resolution Plan
**Logic (`phenotype_mapper.py`):**
1.  **Identify Candidate Alleles:**
    -   Iterate through `GenotypeData.variants`.
    -   Find all star alleles that *could* be present based on identified variants.
    -   Maintain a score or count of matched variants for each star allele definition.
2.  **Select Best Match:**
    -   **Homozygous (`HOM_ALT`)** variants count towards both copies.
    -   **Heterozygous (`HET`)** variants count towards one.
    -   Algorithm: "Star Allele Calling" (similar to Stargazer/PyPGx logic but simplified):
        -   Assume `*1/*1` baseline.
        -   If variants match `*X` definition completely, propose `*X`.
        -   If variants match `*X` and `*Y`, check if they are mutually exclusive sites (same position) or distinct.
3.  **Phasing:**
    -   **Assumption:** If multiple HET variants are found that define distinct alleles (e.g., `varA` -> `*2`, `varB` -> `*3`), assume **trans** (compound heterozygote `*2/*3`) as it is the more conservative/risk-averse assumption for pharmacogenomics (often implies lower function), unless specific phasing info says otherwise or population frequency data suggests `*1/*4` (cis) is more likely. *Note: Strict trans assumption is standard for basic MVPs.*
4.  **Indeterminate Cases:**
    -   If variants are found but don't match any known star allele pattern: `Indeterminate`.
    -   If crucial variants for a candidate `*X` are missing (low coverage): Fallback to `*1` but flag confidence.

### 4. Diplotype → Phenotype Mapping Plan
**Logic (`phenotype_mapper.py`):**
-   **Input:** `Gene`, resolved `Diplotype` (e.g., `*1/*22`).
-   **Direct Lookup:** Query `phenotype_map` from CPIC data.
-   **Activity Score Fallback:**
    -   If direct diplotype lookup fails, use "Activity Value" method (if CPIC provides it):
        -   `Val(*1) = 1.0`, `Val(*22) = 0.5`.
        -   `Total Score = 1.5`.
        -   Map Score -> Phenotype (e.g., 1.0-2.0 -> `NM`).
-   **Output:** Standardized strings: `PM`, `IM`, `NM`, `RM`, `URM`.

### 5. Risk Engine Contract
**File:** `services/pharmacogenomics/risk_engine.py`

```python
from pydantic import BaseModel

class RiskAssessment(BaseModel):
    risk_label: str
    confidence_score: float
    severity: str # "none" | "low" | "moderate" | "high" | "critical"

class ClinicalRecommendation(BaseModel):
    text: str
    implication: str
    recommendation_url: Optional[str]

def evaluate_risk(
    drug: str,
    gene: str,
    phenotype: str,
    diplotype: str
) -> Tuple[RiskAssessment, ClinicalRecommendation]:
    """
    1. Validate drug is supported.
    2. Lookup (Drug, Phenotype) in DrugRecommendations.
    3. If no recommendation (e.g. NM for some drugs), return "Low Risk/Standard Dosing".
    4. Map CPIC "implication" string to severity level.
    """
```

**Severity Mapping:**
- "Actionable PGx" / "High Risk" -> `critical` or `high`.
- "Testing Recommended" -> `moderate`.
- "Informative" -> `low`.
- "No Recommendation" / "Normal" -> `none`.

### 6. Confidence Score Computation
**Formula:** `Confidence = Base * Coverage * Ambiguity`

-   **Base:** 1.0
-   **Coverage Penalty:**
    -   Identify "Key Variants" for the gene (variants defining common actionable alleles *2, *3, *4, etc.).
    -   Check if these positions are in `GenotypeData.covered_positions`.
    -   If missing: `Confidence *= 0.8` per missing key site.
-   **Ambiguity Penalty:**
    -   Unphased Heterozygotes defining two alleles: `Confidence *= 0.9` (Phase uncertainty).
    -   Indeterminate/Partial Match: `Confidence *= 0.5`.
    -   Rare/Unknown Alleles: `Confidence *= 0.7`.

### 7. Multi-Drug Processing Flow
1.  **Ingest VCF:** Parse once -> `Dict[Gene, GenotypeData]`.
2.  **Generate Profile:**
    -   For each `Gene` in CPIC database:
        -   Run `Diplotype Resolution` -> `Diplotype`.
        -   Run `Phenotype Mapping` -> `Phenotype`.
        -   Store in `PatientProfile`.
3.  **Assessment Loop:**
    -   For each `Drug` in input list:
        -   Lookup `Gene` associated with `Drug` (e.g., Warfarin -> CYP2C9, VKORC1).
        -   Retrieve `Phenotype` from profile.
        -   Run `evaluate_risk`.
        -   Construct final JSON object.

### 8. Edge Case Handling
-   **Missing Alleles:** If NO variants found, assume `*1/*1`. If metadata shows no coverage at key sites, return `Indeterminate`.
-   **Unsupported Drug:** Return specific error or "Unknown Drug" risk object (zero confidence).
-   **Partial Genotype Coverage:** If a defining variant for `*4` is missing/no-call, but others present, do NOT call `*4`. Call matching parent or `*1` (`Ambiguous`).
-   **Conflicting Variants:** Start with highest confidence/most specific allele. If conflict (e.g. `*4` and `*5` definitions clash on same position), report `Indeterminate`.

### 9. Unit Testing Scope
-   **`tests/test_phenotype_mapper.py`:**
    -   `test_resolve_homozygous_wildtype`: No variants -> `*1/*1`.
    -   `test_resolve_homozygous_variant`: `HOM_ALT` -> `*2/*2`.
    -   `test_resolve_compound_het`: `HET` var1, `HET` var2 -> `*2/*3`.
-   **`tests/test_risk_engine.py`:**
    -   `test_codeine_pm`: CYP2D6 `PM` -> "Avoid Codeine" (High Severity).
    -   `test_warfarin_sensitivity`: CYP2C9 `*2` -> "Lower Dose".
    -   `test_unsupported_drug`: Returns safe default/error.
-   **`tests/test_indeterminate`:**
    -   Test behavior when coverage is missing.
