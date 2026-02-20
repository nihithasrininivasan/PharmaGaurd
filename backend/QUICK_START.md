# PharmaGaurd Risk Engine 2.0 - Quick Start Guide

## 🚀 5-Minute Quick Start

### Single Drug Risk Assessment

```python
from app.services.pharmacogenomics.risk_engine import RiskEngine

engine = RiskEngine()

risk, recommendation = engine.evaluate_risk(
    drug="clopidogrel",
    gene="CYP2C19",
    phenotype="PM",
    diplotype="*2/*2",
    diplotype_confidence=0.95
)

print(f"Risk Score: {risk.risk_score}/100")        # 82.3
print(f"Risk Level: {risk.risk_level}")            # "high"
print(f"Confidence: {risk.confidence_score}")      # 0.89
print(f"Action: {recommendation.text}")
```

---

### Multi-Drug Polypharmacy Analysis

```python
from app.services.pharmacogenomics.multi_drug_risk import analyze_multi_drug_risk

# Get individual drug assessments first
drug_assessments = engine.evaluate_multiple_drugs(
    drugs=["clopidogrel", "omeprazole", "warfarin"],
    patient_profile=patient_profile
)

# Analyze combined risk
multi_drug = analyze_multi_drug_risk(drug_assessments)

print(f"Combined Risk: {multi_drug.combined_risk_score}/100")
print(f"Interactions: {multi_drug.interaction_count}")
print(f"Recommendation: {multi_drug.polypharmacy_recommendation}")
```

---

### Submit Clinical Feedback

```bash
curl -X POST http://localhost:8000/api/feedback/ \
  -H "Content-Type: application/json" \
  -d '{
    "gene": "CYP2C19",
    "drug": "clopidogrel",
    "reported_diplotype": "*1/*2",
    "correct_diplotype": "*2/*2",
    "feedback_quality": 1.0,
    "comments": "Confirmed by Sanger sequencing"
  }'
```

**Response**:
```json
{
    "status": "success",
    "prior_update": {
        "previous_prior": 1.05,
        "new_prior": 1.095,
        "explanation": "Bayesian update applied"
    }
}
```

---

### Check Drug-Drug Interaction

```bash
curl -X POST http://localhost:8000/api/polypharmacy/check-pair \
  -H "Content-Type: application/json" \
  -d '{
    "drug_a": "warfarin",
    "drug_b": "fluconazole",
    "gene": "CYP2C9"
  }'
```

**Response**: Interaction severity + monitoring recommendations

---

## 📊 Key Outputs

### Risk Assessment Object

```python
{
    "risk_label": "Use Alternative",
    "risk_score": 82.3,              # 0-100 continuous
    "risk_level": "high",             # categorical
    "severity": "high",
    "confidence_score": 0.89,
    "confidence_breakdown": {...}
}
```

### Multi-Drug Assessment

```python
{
    "combined_risk_score": 87.3,
    "combined_risk_level": "major",
    "detected_interactions": [...],
    "critical_warnings": [...],
    "polypharmacy_recommendation": "...",
    "alternative_regimens": [...]
}
```

---

## 🎯 Common Use Cases

### 1. Pre-prescription Screening
```python
# Check if drug is safe for patient genotype
risk, rec = engine.evaluate_risk(...)
if risk.risk_score >= 90:
    print("CONTRAINDICATED - Choose alternative")
elif risk.risk_score >= 70:
    print("HIGH RISK - Review before prescribing")
```

### 2. Polypharmacy Review
```python
# Analyze patient's current medication list
multi_drug = analyze_multi_drug_risk(drug_assessments)
if multi_drug.critical_warnings:
    print("URGENT: Regimen modification required")
```

### 3. Continuous Learning
```python
# Submit correction after confirmation
POST /api/feedback/
# Automatically updates model via Bayesian learning
```

---

## 🔧 Configuration

### Enable/Disable Features

```python
engine = RiskEngine(
    enable_feedback_learning=True,   # Use learning_priors.json
    enable_calibration=True,          # Apply confidence calibration
)
```

### Custom Interaction Database

```python
from app.services.pharmacogenomics.multi_drug_risk import (
    MultiDrugRiskAnalyzer,
    DrugDrugInteraction
)

custom_interactions = [
    DrugDrugInteraction(
        drug_a="drug_x",
        drug_b="drug_y",
        gene="CYP3A5",
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=2.0,
        ...
    )
]

analyzer = MultiDrugRiskAnalyzer(interaction_db=custom_interactions)
```

---

## 📁 File Locations

```
data/
├── learning_priors.json        # Feedback learning state (auto-updated)
├── calibration_data.json       # Performance metrics
└── model_versions/             # Version snapshots
    ├── 2.0.0.json
    └── ...
```

---

## 🚨 Critical Warnings Detection

```python
# System automatically flags critical situations
if risk.risk_score >= 90:
    # CRITICAL - Contraindicated

if "CRITICAL INTERACTION" in multi_drug.critical_warnings:
    # CRITICAL - Drug-drug interaction

if multi_drug.monitoring_priority == "immediate_action_required":
    # URGENT - Immediate review needed
```

---

## 📖 API Endpoints

### Polypharmacy Analysis
- `POST /api/polypharmacy/analyze` — Full multi-drug analysis
- `POST /api/polypharmacy/check-pair` — Check specific drug pair
- `GET /api/polypharmacy/interactions` — Query database

### Feedback Learning
- `POST /api/feedback/` — Submit clinical correction

### Interaction Database
- `GET /api/polypharmacy/interaction-summary` — Statistics

---

## 🎓 Risk Score Interpretation

| Score | Level | Action |
|-------|-------|--------|
| 90-100 | Critical | **AVOID** — Contraindicated |
| 70-89 | High | Consider alternative first-line |
| 40-69 | Moderate | Adjust dosing, monitor closely |
| 20-39 | Low | Standard dosing with awareness |
| 0-19 | None | No pharmacogenomic concerns |

---

## ⚙️ Model Version Management

```python
from app.services.pharmacogenomics.model_versioning import ModelVersionManager

manager = ModelVersionManager()

# Create new version
new = manager.create_new_version(
    from_version="2.0.0",
    version_increment="minor",
    notes="Adjusted penalties based on calibration"
)

# Rollback if needed
manager.rollback(target_version="2.0.0")

# Compare versions
comparison = manager.compare_versions("2.0.0", "2.1.0")
```

---

## 🔍 Debugging & Monitoring

### Check Feedback Learning Status
```python
from app.services.pharmacogenomics.feedback_learning import load_learning_priors

priors = load_learning_priors()
print(f"Total feedback events: {priors.metadata['total_feedback_events']}")
print(f"Last updated: {priors.metadata['last_updated']}")
print(f"CYP2C19 *2/*2 prior: {priors.get_diplotype_prior('CYP2C19', '*2/*2')}")
```

### Check Model Calibration
```python
from app.services.pharmacogenomics.model_calibration import (
    ConfidenceCalibrator,
    DriftDetector
)

calibrator = ConfidenceCalibrator()
stats = calibrator.get_calibration_stats()
for s in stats:
    print(f"{s.bin_label}: {s.empirical_accuracy:.2%} accuracy")
```

---

## 🧪 Testing

```bash
# Run unit tests
pytest tests/services/pharmacogenomics/

# Test multi-drug analysis
pytest tests/services/pharmacogenomics/test_multi_drug_risk.py

# Test API endpoints
pytest tests/api/routes/test_polypharmacy.py
```

---

## 📚 Full Documentation

- **RISK_ENGINE_AUDIT** — Initial audit findings
- **MULTI_DRUG_RISK_INTEGRATION.md** — Complete multi-drug guide
- **IMPLEMENTATION_SUMMARY.md** — Full technical overview

---

## 🆘 Troubleshooting

### Feedback not applying?
```python
# Check if feedback learning is enabled
engine = RiskEngine(enable_feedback_learning=True)
```

### Interactions not detected?
```python
# Verify gene matches between drugs
# Check phenotype is in affected_phenotypes list
```

### Risk scores seem wrong?
```python
# Check confidence calibration is enabled
engine = RiskEngine(enable_calibration=True)
```

---

## 📞 Support

- Documentation: `/docs/`
- GitHub Issues: Report bugs and feature requests
- Email: pharmacogenomics@pharmaguard.io

---

**Version**: 2.0.0
**Last Updated**: 2026-02-19
