# CPIC Implementation Improvements

## Relevant Skills
-   **Excel Processing**: `.agent/skills/xlsx/SKILL.md`
-   **Robust Testing**: `.agent/skills/test-driven-development/SKILL.md`

## Critical Priority (Must Fix for Accuracy)

### 1. Robust Variant Normalization & Parsing
**Current State:** Variants are matched by exact string: `POS:REF:ALT`.
**Problem:** VCFs use various representations for indels (left-alignment, padding bases) and structural variants (e.g., `<DEL>` for *5 gene deletions). Exact string matching will fail silently for valid variants, leading to incorrect *1/*1 calls.
**Improvement Plan:**
1.  **Library Selection**: Use `biocommons/hgvs` or `vcfpy` for rigorous VCF parsing. Do NOT write custom regex parsers for VCF lines.
2.  **Normalization Workflow**:
    -   Ingest VCF record.
    -   Convert to SPDI (Sequence Position Deletion Insertion) or similar canonical format for comparison.
    -   *Crucial*: Left-align all indels to match CPIC's likely alignment (usually based on limited reference context).
3.  **Structural Variant (SV) Handling**:
    -   Scan VCF `ALT` column for symbolic alleles: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`.
    -   Scan `INFO` field for `SVTYPE` tags.
    -   Map `<DEL>` covering the entire gene region (check `POS` and `END` tags) to the gene deletion allele (e.g., CYP2D6 `*5`).
    -   Map `<DUP>` to gene duplication (e.g., `xN` copy number).

**Relevant Skill:** `test-driven-development`
-   **Step 1 (Red)**: Create a test file `tests/test_variant_normalization.py`.
-   **Step 2 (Red)**: Add test cases for:
    -   Standard SNP: `chr22 42522611 . C G` (Should match)
    -   Padding Base Indel: `chr22 42522611 TC T` (Should match deletion of C at 42522612)
    -   Symbolic Deletion: `chr22 42522000 . N <DEL>` (Should match *5 if overlaps gene)
-   **Step 3 (Green)**: Implement the normalization logic until these match your internal CPIC representation.

### 2. Genomic Build Validation (GRCh37 vs GRCh38)
**Current State:** Assumes input VCF matches the build of the CPIC source files (usually GRCh38 nowadays, but legacy data exists).
**Problem:** Mixing builds leads to silent errors (variants at wrong positions).
**Improvement Plan:**
1.  **Header Check**: Parse VCF header for `##reference=file:///.../GRCh38.fa` or `##contig=<ID=chr1,...>` tags.
2.  **Validation**:
    -   If CPIC data is GRCh38 (check source metadata), enforce VCF is GRCh38.
    -   Raise explicit `ValueError: Incompatible genomic build` if mismatch.
3.  **Liftover (Optional but Recommended)**: Use `pyliftover` if you need to support GRCh37 inputs against GRCh38 data.

### 3. Data-Driven Phenotype & Drug Mapping
**Current State:** `cpic_etl.py` hardcodes phenotype mappings (e.g., `*1/*1` -> `NM`) and drug guidelines.
**Improvement Plan:**
1.  **Phenotypes**:
    -   Target File: `*diplotype_phenotype.xlsx` in CPIC data.
    -   Logic: Read `Diplotype` and `Phenotype` columns.
    -   *Caution*: Handle wildcards if CPIC uses them (e.g., `*1/*x` -> `NM`).
2.  **Drug Guidelines**:
    -   Target File: `*CDS.xlsx` (Clinical Decision Support).
    -   Logic:
        -   Map `Phenotype` column to `Recommendation` logic.
        -   Extract `Implication` text.
        -   Extract `Classification` for severity mapping (Strong, Moderate, Optional -> High, Moderate, Low).

**Relevant Skill:** `xlsx`
-   **Library**: Use `pandas`.
-   **Pattern**:
    ```python
    df = pd.read_excel(file_path, sheet_name='Diplotype_Phenotype')
    # Clean column names (strip whitespace, lower case)
    df.columns = [c.strip().lower() for c in df.columns]
    # Iterate and build efficient lookup dict
    phenotype_map = dict(zip(df['diplotype'], df['phenotype']))
    ```
-   **Validation**: Ensure no duplicate keys in the map. Raise error if multiple phenotypes map to same diplotype.

### 4. Activity Score Data
**Current State:** Activity scores are hardcoded in `cpic_loader.py` with simplified logic.
**Improvement Plan:**
1.  **Source**: `*allele_functionality_reference.xlsx` (if available) or `Allele Functionality` sheet in definition table.
2.  **Logic**:
    -   Map Allele -> Clinical Function (e.g., `*4` -> `No Function`).
    -   Map Function -> Score (e.g., `No Function` -> 0, `Decreased` -> 0.5, `Normal` -> 1.0, `Increased` -> >1.0).
    -   Implementation: `get_activity_score(allele)` should query this lookup.

**Relevant Skill:** `xlsx`
-   Use `pandas` to load this map during ETL and save to `cpic_cache.json`.

### 5. Quality Control & Validation
**Current State:** Basic Pydantic type checks only.
**Improvement Plan:**
1.  **Intermediate Model Update**: Add validators to `VariantCall`.
2.  **Logic**:
    -   `if quality < 20`: Filter out or flag as `LOW_QUALITY`.
    -   `if ad and sum(ad) < 10`: Warn "Low read depth".
    -   `if zygosity is None`: Raise error "Ambiguous genotype".

**Relevant Skill:** `test-driven-development`
-   **Test Case**: Create a "dirty" VCF snippet with `QUAL=10` and ensure your parser filters it or marks it.

### 6. Observability & Logging
**Current State:** Minimal logging.
**Improvement Plan:**
-   Use Python's `logging` module or `structlog`.
-   **Trace ID**: Generate a unique ID per `evaluate_risk` call.
-   **Log Events**:
    -   "Loaded CPIC cache v1.0"
    -   "Parsing VCF for sample X..."
    -   "Gene CYP2D6: Matched variants [rs123, rs456] -> Alleles [*1, *4]"
    -   "Phenotype Map: *1/*4 -> IM"
    -   "Risk Assessment: Codeine -> Moderate Risk (IM)"

## Medium Priority

### 7. Phasing Logic
**Current State:** Assumes *trans* phase for compound hets.
**Improvement Plan:**
1.  **Check PS Tag**: If identifiers match (e.g., `0|1:12345` and `1|0:12345`), variants are phased relative to each other.
2.  **Logic**:
    -   Same Phase Set + Same Haplotype (0|1 and 0|1) -> *Cis* (Same chromosome).
    -   Same Phase Set + Diff Haplotype (0|1 and 1|0) -> *Trans*.
3.  **Fallback**: If unphased (`0/1`), keep current *trans* assumption (conservative).

### 8. Granular Indeterminate States
**Current State:** Returns "Indeterminate".
**Improvement Plan:**
1.  **Enum Definition**: Define `ResolutionStatus` enum: `SUCCESS`, `NO_COVERAGE`, `AMBIGUOUS`, `NOVEL_VARIANTS`.
2.  **Logic**:
    -   If no variants found AND no coverage at key sites -> `NO_COVERAGE`.
    -   If variants match `*1/*4` AND `*2/*5` equally validly -> `AMBIGUOUS`.
    -   If variants found but don't match any star allele -> `NOVEL_VARIANTS` (default to Indeterminate Risk).

**Relevant Skill:** `test-driven-development`
-   Write tests expecting these specific Enum returns, not just a string "Indeterminate".

### 9. Configurable Penalties
**Current State:** Hardcoded scaling factors.
**Improvement Plan:**
-   Create `config.py` with `CONFIDENCE_PENALTY_NO_PHASING = 0.8`, `CONFIDENCE_PENALTY_PARTIAL_MATCH = 0.7`.
-   Import these constants in `risk_engine.py`.

## Low Priority

### 10. Auto-ETL on Startup
**Current State:** Manual run required.
**Improvement Plan:**
-   In `cpic_loader.py` `__init__`:
    -   `if not cache_file.exists():`
    -   `log.info("Cache missing, running ETL...")`
    -   `from app.utils.cpic_etl import main as run_etl`
    -   `run_etl()`