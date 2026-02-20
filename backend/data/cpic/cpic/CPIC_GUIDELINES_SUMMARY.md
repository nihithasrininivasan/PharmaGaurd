# CPIC Guidelines Summary for Pharmacogenomics
**Source**: Clinical Pharmacogenetics Implementation Consortium (CPIC)
**Date**: Feb 2026 (Live Retrieval)
**Purpose**: System context for LLM explanation generation.

---

## 1. Codeine & CYP2D6
**Gene**: `CYP2D6` (Drug Metabolizer)
**Drug**: Codeine (bioactivated to morphine)

| Phenotype | Activity Score | Implication | Recommendation | Strength |
|---|---|---|---|---|
| **Ultrarapid Metabolizer (UM)** | > 2.25 | Increased morphine formation. High risk of toxicity. | **Avoid**. Use non-tramadol alternative. | Strong |
| **Normal Metabolizer (NM)** | 1.25 – 2.25 | Normal bioactivation. | Use label recommended dosage. | Strong |
| **Intermediate Metabolizer (IM)** | 0.25 – 1.0 | Reduced bioactivation (less analgesia). | Use label dosage. Monitor for lack of efficacy. | Moderate |
| **Poor Metabolizer (PM)** | 0 | No bioactivation (no analgesia). | **Avoid**. Use non-tramadol alternative. | Strong |

---

## 2. Clopidogrel & CYP2C19
**Gene**: `CYP2C19` (Bioactivator)
**Drug**: Clopidogrel (Plavix) – Antiplatelet

| Phenotype | Implication | Recommendation | Strength |
|---|---|---|---|
| **Ultrarapid / Rapid (UM/RM)** | Increased bioactivation. | Use label recommended dosage. | Strong |
| **Normal Metabolizer (NM)** | Normal bioactivation. | Use label recommended dosage. | Strong |
| **Intermediate Metabolizer (IM)** | Reduced bioactivation. Low platelet inhibition. | **Avoid**. Use alternative (e.g., prasugrel, ticagrelor). | Strong |
| **Poor Metabolizer (PM)** | No bioactivation. High clotting risk. | **Avoid**. Use alternative (e.g., prasugrel, ticagrelor). | Strong |

---

## 3. Warfarin & CYP2C9
**Gene**: `CYP2C9` (Clearance) + *VKORC1* (Target)
**Drug**: Warfarin – Anticoagulant

*Note: Dosing is highly individual and often algorithm-based.*

| CYP2C9 Phenotype | Implication | Recommendation | Strength |
|---|---|---|---|
| **Normal (*1/*1)** | Normal clearance. | Use label dosage / standard algorithm. | Strong |
| **Intermediate (*1/*2, *1/*3)** | Reduced clearance (30-40% less). | **Reduce dose**. Calculate using genetics-based algorithm. | Strong |
| **Poor (*2/*2, *3/*3 etc)** | Severely reduced clearance (>80% less). High bleed risk. | **Greatly reduce dose** or consider alternative anticoagulant. | Strong |

---

## 4. Simvastatin & SLCO1B1
**Gene**: `SLCO1B1` (Transporter – OATP1B1)
**Drug**: Simvastatin – Statin (Cholesterol)

| Phenotype | Implication | Recommendation | Strength |
|---|---|---|---|
| **Normal Function** | Normal transport. | Desired starting dose. | Strong |
| **Decreased Function** | Reduced hepatic uptake. Increased systemic exposure (myopathy risk). | **Limit max dose to 20 mg/day**. Or use alternative statin (e.g., rosuvastatin). | Strong |
| **Poor Function** | High systemic exposure. High myopathy risk. | **Prescribe Alternative**. (e.g., rosuvastatin, pravastatin). | Strong |

---

## 5. Azathioprine & TPMT
**Gene**: `TPMT` (Clearance methyltransferase)
**Drug**: Azathioprine (Immunosuppressant)

| Phenotype | Implication | Recommendation | Strength |
|---|---|---|---|
| **Normal Metabolizer** | Normal clearance. | Start with normal dose (2–3 mg/kg/day). | Strong |
| **Intermediate Metabolizer** | Reduced clearance. High TGN metabolites. | **Reduce dose by 30–80%** (e.g. 0.6–2.4 mg/kg/day). | Strong |
| **Poor Metabolizer** | Zero clearance. Fatal myelosuppression risk. | **Avoid** or reduce dose by 10-fold and dose 3x/week. | Strong |

---

## 6. Fluorouracil (5-FU) & DPYD
**Gene**: `DPYD` (Clearance enzyme)
**Drug**: Fluorouracil / Capecitabine (Chemotherapy)

| Activity Score | Phenotype | Recommendation | Strength |
|---|---|---|---|
| **Score ≥ 2 (Normal)** | Normal clearance. | Use label dosage. | Strong |
| **Score 1 – 1.5 (Intermediate)** | Reduced clearance. 50% risk of severe toxicity. | **Reduce dose by 50%**. Titrate based on tolerance. | Strong |
| **Score 0 – 0.5 (Poor)** | Little/No clearance. | **Avoid**. Use non-fluoropyrimidine drug. | Strong |

---
**Disclaimer**: This document serves as a contextual guide for the AI Explanation Model. Clinical decisions should always be made by a qualified healthcare professional using the most up-to-date guidelines at [cpicpgx.org](https://cpicpgx.org).
