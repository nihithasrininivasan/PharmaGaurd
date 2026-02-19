# Pharmacogenomics Service

CPIC-aligned pharmacogenomic decision engine for deterministic drug risk assessment.

## Overview

This service implements a **deterministic, rule-based** backend medical logic layer for pharmacogenomic risk assessment based on CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines. It ingests VCF data, resolves diplotypes/phenotypes, and generates risk assessments for specific drugs.

### Key Features

- ✅ **Deterministic Logic**: No LLMs - all risk classification is rule-based
- ✅ **CPIC-Aligned**: Based on official CPIC allele definitions and guidelines
- ✅ **Unit Testable**: Comprehensive test coverage for medical logic
- ✅ **Confidence Scoring**: Transparent confidence metrics based on coverage and ambiguity
- ✅ **Multi-Gene Support**: Currently supports CYP2D6, CYP2C19, CYP2C9, TPMT, SLCO1B1

## Architecture

```
┌─────────────────┐
│   VCF Data      │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────┐
│  Phenotype Mapper           │
│  - Diplotype Resolution     │
│  - Star Allele Calling      │
│  - Phenotype Mapping        │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│  Risk Engine                │
│  - Drug-Gene Lookup         │
│  - CPIC Recommendations     │
│  - Confidence Scoring       │
└────────┬────────────────────┘
         │
         ▼
┌─────────────────────────────┐
│  DrugAssessment Output      │
│  - Risk Label               │
│  - Severity Level           │
│  - Clinical Recommendation  │
└─────────────────────────────┘
```

## Quick Start

### 1. Run ETL to Process CPIC Data

```bash
cd backend
python3 app/utils/cpic_etl.py
```

This generates `data/cpic_cache.json` from the Excel files in `data/cpic/`.

### 2. Use the Service

```python
from app.services.pharmacogenomics import (
    VariantCall,
    GenotypeData,
    PhenotypeMapper,
    RiskEngine,
    PatientProfile
)

# Initialize services
mapper = PhenotypeMapper()
engine = RiskEngine()

# Create genotype data from VCF
genotype = GenotypeData(
    sample_id="PATIENT001",
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
        )
    ],
    coverage_mean=55.0,
    covered_positions=[42126611]
)

# Resolve diplotype and phenotype
diplotype_result = mapper.process_genotype(genotype)
print(f"Diplotype: {diplotype_result.diplotype}")
print(f"Phenotype: {diplotype_result.phenotype}")
print(f"Confidence: {diplotype_result.confidence}")

# Evaluate drug risk
risk, recommendation = engine.evaluate_risk(
    drug="codeine",
    gene="CYP2D6",
    phenotype=diplotype_result.phenotype,
    diplotype=diplotype_result.diplotype,
    diplotype_confidence=diplotype_result.confidence
)

print(f"Risk: {risk.risk_label}")
print(f"Severity: {risk.severity}")
print(f"Recommendation: {recommendation.text}")
```

## Data Models

### Input Models

**VariantCall**: Single variant from VCF
```python
VariantCall(
    chrom="chr22",
    pos=42126611,
    rsid="rs1135840",
    ref="C",
    alt="G",
    zygosity="HET",  # HET, HOM_REF, HOM_ALT
    quality=99.0,
    filter="PASS"
)
```

**GenotypeData**: Gene-specific genotype information
```python
GenotypeData(
    sample_id="PATIENT001",
    gene_symbol="CYP2D6",
    variants=[...],
    coverage_mean=50.0,
    covered_positions=[...]
)
```

### Output Models

**DiplotypeResult**: Resolved diplotype and phenotype
```python
DiplotypeResult(
    gene="CYP2D6",
    diplotype="*1/*4",
    phenotype="IM",
    confidence=0.85,
    is_indeterminate=False,
    notes="Heterozygous with wildtype"
)
```

**DrugAssessment**: Complete risk assessment
```python
DrugAssessment(
    drug="codeine",
    gene="CYP2D6",
    diplotype="*4/*4",
    phenotype="PM",
    risk=RiskAssessment(...),
    recommendation=ClinicalRecommendation(...)
)
```

## Diplotype Resolution Logic

### Algorithm

1. **Identify Candidate Alleles**: Match observed variants to star allele definitions
2. **Score Alleles**: Calculate match scores based on zygosity (HET=1, HOM_ALT=2)
3. **Select Best Diplotype**:
   - Homozygous variant → `*X/*X`
   - Single allele → `*1/*X` (heterozygous with wildtype)
   - Two alleles → `*X/*Y` (compound heterozygote, trans assumed)
4. **Map to Phenotype**: Use CPIC phenotype map or activity score method
5. **Calculate Confidence**: Adjust based on coverage and ambiguity

### Phenotype Codes

- **PM** - Poor Metabolizer (no/minimal enzyme activity)
- **IM** - Intermediate Metabolizer (reduced activity)
- **NM** - Normal Metabolizer (normal activity)
- **RM** - Rapid Metabolizer (increased activity)
- **UM** - Ultra-rapid Metabolizer (very high activity)

## Confidence Scoring

```
Confidence = Base * Coverage * Ambiguity
```

**Coverage Penalty**: Missing key variant positions → 0.8 per missing site

**Ambiguity Penalty**:
- Unphased heterozygotes → 0.9
- Indeterminate/partial match → 0.5
- Rare/unknown alleles → 0.7

## Severity Levels

| Severity | Meaning | Example |
|----------|---------|---------|
| **critical** | Contraindicated, life-threatening | Avoid codeine in PM |
| **high** | Major risk, therapeutic failure | Warfarin in PM |
| **moderate** | Dose adjustment needed | Clopidogrel in IM |
| **low** | Informative, minor impact | - |
| **none** | Standard dosing | Normal metabolizers |

## Supported Drugs

| Drug | Gene | Phenotypes | Guideline |
|------|------|-----------|-----------|
| Codeine | CYP2D6 | PM, IM, NM, UM | Avoid in PM/UM |
| Warfarin | CYP2C9 | PM, IM, NM | Reduce dose in PM/IM |
| Clopidogrel | CYP2C19 | PM, IM, NM, UM | Alternative in PM |

## Testing

### Run Tests

```bash
# Run all tests
python3 -m pytest tests/services/pharmacogenomics/ -v

# Run specific test file
python3 -m pytest tests/services/pharmacogenomics/test_phenotype_mapper.py -v

# Run with coverage
python3 -m pytest tests/services/pharmacogenomics/ --cov=app/services/pharmacogenomics
```

### Test Categories

1. **test_phenotype_mapper.py**: Diplotype resolution logic
2. **test_risk_engine.py**: Risk assessment and recommendations
3. **test_integration.py**: End-to-end workflows

## Edge Cases

### Missing Alleles
- No variants found → Assume `*1/*1` (wildtype)
- No coverage at key sites → Return `Indeterminate`

### Unsupported Drug
- Return error with `confidence_score = 0.0`

### Partial Coverage
- Defining variant missing → Don't call that allele
- Reduce confidence score based on missing positions

### Conflicting Variants
- Start with most specific allele
- Report `Indeterminate` if true conflict exists

## File Structure

```
app/services/pharmacogenomics/
├── __init__.py              # Package exports
├── models.py                # Pydantic data models
├── cpic_loader.py           # CPIC data singleton loader
├── phenotype_mapper.py      # Diplotype resolution
├── risk_engine.py           # Risk assessment engine
├── IMPLEMENTATION_PLAN.md   # Detailed implementation plan
└── README.md                # This file

app/utils/
└── cpic_etl.py              # ETL script for CPIC data

data/
├── cpic/                    # Raw CPIC Excel files
│   ├── CYP2D6/
│   ├── CYP2C19/
│   ├── CYP2C9/
│   ├── TPMT/
│   └── SLCO1B1/
└── cpic_cache.json          # Processed CPIC data

tests/services/pharmacogenomics/
├── test_phenotype_mapper.py
├── test_risk_engine.py
└── test_integration.py
```

## References

- [CPIC Guidelines](https://cpicpgx.org/guidelines/)
- [PharmVar](https://www.pharmvar.org/) - Star allele nomenclature
- [CPIC Allele Definition Tables](https://github.com/cpicpgx/cpic-data)

## Future Enhancements

- [ ] Add more genes (DPYD, UGT1A1, etc.)
- [ ] Support for CNVs and structural variants
- [ ] Activity score-based phenotype calculation for all genes
- [ ] VCF parser integration
- [ ] Phasing information support
- [ ] Population frequency priors for ambiguous calls

## License

This implementation follows CPIC guidelines which are freely available for clinical use.
