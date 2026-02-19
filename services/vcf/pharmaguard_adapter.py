from __future__ import annotations

import datetime as _dt
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .parser import VcfParseResult, parse_vcf
from .variant_extractor import ExtractedVariant, extract_pharmacogenes
from .cpic_tables import infer_phenotype_from_variants


# ---------------------------------------------------------------------------
# Drug → primary gene mapping
# ---------------------------------------------------------------------------
SUPPORTED_DRUGS_TO_GENE: Mapping[str, str] = {
    "CODEINE":     "CYP2D6",
    "WARFARIN":    "CYP2C9",
    "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1",
    "AZATHIOPRINE":"TPMT",
    "FLUOROURACIL":"DPYD",
}

# ---------------------------------------------------------------------------
# CPIC drug × phenotype recommendation table
# All phenotype strings must match those returned by cpic_tables.py
# ---------------------------------------------------------------------------
_URL = {
    "CYP2D6":  "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
    "CYP2C19": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
    "CYP2C9":  "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
    "SLCO1B1": "https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/",
    "TPMT":    "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/",
    "DPYD":    "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/",
}

# Structure: drug_upper → phenotype_string → {risk_label, confidence_score, severity, text, implication}
DRUG_PHENOTYPE_RECS: Dict[str, Dict[str, Dict]] = {

    # ── CODEINE / CYP2D6 ────────────────────────────────────────────────────
    "CODEINE": {
        "Ultrarapid Metabolizer": dict(
            risk_label="Avoid – Toxicity Risk",
            confidence_score=0.97,
            severity="critical",
            text="Avoid codeine. Ultrarapid CYP2D6 metabolizers convert codeine to morphine too rapidly, causing life-threatening opioid toxicity.",
            implication="Excessive morphine production; risk of respiratory depression and death.",
        ),
        "Normal Metabolizer": dict(
            risk_label="Use label-recommended dosage",
            confidence_score=0.95,
            severity="none",
            text="Initiate therapy with label-recommended codeine dosage.",
            implication="Normal codeine-to-morphine conversion expected.",
        ),
        "Intermediate Metabolizer": dict(
            risk_label="Use label dosage – monitor",
            confidence_score=0.80,
            severity="low",
            text="Initiate therapy with label-recommended dosage. Monitor for reduced efficacy.",
            implication="Reduced but likely adequate analgesic effect.",
        ),
        "Poor Metabolizer": dict(
            risk_label="Avoid – Ineffective",
            confidence_score=0.95,
            severity="moderate",
            text="Avoid codeine. Poor CYP2D6 metabolizers produce insufficient morphine; use an alternative analgesic.",
            implication="Inadequate pain relief; ineffective therapy.",
        ),
    },

    # ── WARFARIN / CYP2C9 ───────────────────────────────────────────────────
    "WARFARIN": {
        "Normal Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.93,
            severity="none",
            text="Initiate warfarin at standard doses and adjust based on INR monitoring.",
            implication="Normal warfarin clearance expected.",
        ),
        "Intermediate Metabolizer": dict(
            risk_label="Reduce starting dose",
            confidence_score=0.90,
            severity="moderate",
            text="Reduce starting warfarin dose by 25–50%. Monitor INR closely and adjust accordingly.",
            implication="Reduced CYP2C9 activity increases warfarin exposure and bleeding risk.",
        ),
        "Poor Metabolizer": dict(
            risk_label="Significantly reduce dose",
            confidence_score=0.95,
            severity="high",
            text="Reduce starting warfarin dose by ≥50%. Increase INR monitoring frequency. Consider alternative anticoagulant.",
            implication="Severely reduced CYP2C9 metabolism; major bleeding risk if standard doses used.",
        ),
    },

    # ── CLOPIDOGREL / CYP2C19 ───────────────────────────────────────────────
    "CLOPIDOGREL": {
        "Ultrarapid Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.88,
            severity="none",
            text="Initiate clopidogrel at standard dosage.",
            implication="Enhanced conversion to active metabolite; efficacy maintained or increased.",
        ),
        "Rapid Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.88,
            severity="none",
            text="Initiate clopidogrel at standard dosage.",
            implication="Normal or slightly increased clopidogrel activation expected.",
        ),
        "Normal Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.93,
            severity="none",
            text="Initiate clopidogrel at standard dosage.",
            implication="Normal clopidogrel activation expected.",
        ),
        "Intermediate Metabolizer": dict(
            risk_label="Consider alternative antiplatelet",
            confidence_score=0.85,
            severity="moderate",
            text="Consider an alternative antiplatelet agent (e.g., prasugrel, ticagrelor) where clinically appropriate. If clopidogrel used, monitor platelet function.",
            implication="Reduced CYP2C19 activation of clopidogrel; decreased antiplatelet effect.",
        ),
        "Poor Metabolizer": dict(
            risk_label="Avoid – use alternative antiplatelet",
            confidence_score=0.97,
            severity="high",
            text="Avoid clopidogrel. Use an alternative antiplatelet agent (e.g., prasugrel, ticagrelor) per clinical guidelines.",
            implication="Minimal clopidogrel activation; high risk of cardiovascular events.",
        ),
    },

    # ── SIMVASTATIN / SLCO1B1 ───────────────────────────────────────────────
    "SIMVASTATIN": {
        "Normal Function": dict(
            risk_label="Standard dosing",
            confidence_score=0.90,
            severity="none",
            text="Prescribe desired simvastatin dose per standard guidelines.",
            implication="Normal simvastatin hepatic uptake; standard myopathy risk.",
        ),
        "Decreased Function": dict(
            risk_label="Use lowest effective dose or switch statin",
            confidence_score=0.88,
            severity="moderate",
            text="Limit simvastatin dose to ≤20 mg/day or consider an alternative statin (e.g., rosuvastatin, pravastatin).",
            implication="Reduced SLCO1B1 transport increases simvastatin plasma levels; elevated myopathy risk.",
        ),
        "Poor Function": dict(
            risk_label="Avoid high-dose simvastatin",
            confidence_score=0.95,
            severity="high",
            text="Avoid doses >20 mg/day simvastatin. Consider an alternative statin with lower SLCO1B1 dependence.",
            implication="Severely reduced hepatic uptake; high risk of simvastatin-induced myopathy.",
        ),
        "Increased Function": dict(
            risk_label="Standard dosing",
            confidence_score=0.80,
            severity="none",
            text="Prescribe desired simvastatin dose per standard guidelines.",
            implication="Possibly enhanced hepatic uptake; no dose adjustment expected.",
        ),
    },

    # ── AZATHIOPRINE / TPMT ─────────────────────────────────────────────────
    "AZATHIOPRINE": {
        "Normal Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.93,
            severity="none",
            text="Initiate azathioprine at normal doses (1.5–3 mg/kg/day). Monitor CBC regularly.",
            implication="Normal thiopurine methylation; standard myelosuppression risk.",
        ),
        "Intermediate Metabolizer": dict(
            risk_label="Reduce dose by 30–70%",
            confidence_score=0.90,
            severity="moderate",
            text="Reduce azathioprine starting dose by 30–70% and titrate based on toxicity and efficacy. Monitor CBC weekly for first month.",
            implication="Reduced TPMT activity; increased thioguanine nucleotide accumulation and myelosuppression risk.",
        ),
        "Poor Metabolizer": dict(
            risk_label="Avoid or use drastically reduced dose",
            confidence_score=0.97,
            severity="critical",
            text="Avoid conventional azathioprine dosing. If use is necessary, consider extremely low doses (≤10% of normal) with close haematological monitoring. Consider alternative immunosuppressant.",
            implication="Absent TPMT activity; risk of life-threatening myelosuppression.",
        ),
    },

    # ── FLUOROURACIL / DPYD ─────────────────────────────────────────────────
    "FLUOROURACIL": {
        "Normal Metabolizer": dict(
            risk_label="Standard dosing",
            confidence_score=0.93,
            severity="none",
            text="Initiate fluorouracil at standard doses per protocol.",
            implication="Normal DPD-mediated fluoropyrimidine catabolism expected.",
        ),
        "Intermediate Metabolizer": dict(
            risk_label="Reduce starting dose by 25–50%",
            confidence_score=0.90,
            severity="moderate",
            text="Reduce starting fluorouracil dose by 25–50% followed by titration based on toxicity. Increase monitoring frequency.",
            implication="Reduced DPYD activity; increased fluoropyrimidine exposure and toxicity risk.",
        ),
        "Poor Metabolizer": dict(
            risk_label="Avoid fluorouracil",
            confidence_score=0.97,
            severity="critical",
            text="Avoid fluorouracil and other fluoropyrimidines. If therapy is essential, use ≤25% of normal dose under strict monitoring.",
            implication="Near-absent DPD activity; severe or fatal fluoropyrimidine toxicity risk.",
        ),
    },
}

# Default fallback for unrecognised phenotype
_FALLBACK_REC = dict(
    risk_label="No specific CPIC recommendation available",
    confidence_score=0.50,
    severity="low",
    text="No CPIC recommendation found for this drug–phenotype combination. Use clinical judgment.",
    implication="Phenotype–drug combination not covered by current CPIC guidelines.",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@dataclass
class RiskAssessment:
    risk_label: str
    confidence_score: float
    severity: str


def _normalise_drug_name(name: str) -> str:
    return name.strip().upper()


def _get_rec(drug_upper: str, phenotype: str) -> Dict:
    return (
        DRUG_PHENOTYPE_RECS.get(drug_upper, {}).get(phenotype)
        or _FALLBACK_REC
    )


def _build_detected_variants(variants: Sequence[ExtractedVariant]) -> List[Dict]:
    out: List[Dict] = []
    for v in variants:
        out.append(
            {
                "rsid":     v.rsid or None,
                "star":     v.star or None,
                "chrom":    v.chrom,
                "pos":      v.pos,
                "ref":      v.ref,
                "alt":      v.alt,
                "gt":       v.gt or None,
                "zygosity": v.zygosity,
                "gene":     v.gene,
                "info":     v.raw_info,
            }
        )
    return out


def _now_iso() -> str:
    return _dt.datetime.now(tz=_dt.timezone.utc).isoformat()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def analyze_vcf_for_drugs(
    file_bytes: bytes,
    drugs: Iterable[str],
) -> List[Dict]:
    """
    Parse VCF bytes and produce a CPIC-informed risk report per drug.

    Returns a list of dicts (one per drug) containing:
      patient_id, drug, timestamp, risk_assessment,
      pharmacogenomic_profile, clinical_recommendation,
      llm_generated_explanation, quality_metrics
    """
    parsed: VcfParseResult = parse_vcf(file_bytes)
    by_gene = extract_pharmacogenes(parsed.variants)

    results: List[Dict] = []
    timestamp = _now_iso()

    for raw_drug in drugs:
        drug_upper = _normalise_drug_name(raw_drug)
        primary_gene = SUPPORTED_DRUGS_TO_GENE.get(drug_upper, "")
        variants: Sequence[ExtractedVariant] = by_gene.get(primary_gene, [])

        # ── CPIC inference ──────────────────────────────────────────────────
        if primary_gene:
            diplotype, phenotype = infer_phenotype_from_variants(primary_gene, variants)
        else:
            diplotype, phenotype = "Unknown", "Unknown"

        rec = _get_rec(drug_upper, phenotype)
        gene_url = _URL.get(primary_gene)  # type: ignore[attr-defined]

        risk_assessment = {
            "risk_label":       rec["risk_label"],
            "confidence_score": rec["confidence_score"],
            "severity":         rec["severity"],
        }

        pharmacogenomic_profile = {
            "primary_gene":      primary_gene or "Unknown",
            "diplotype":         diplotype,
            "phenotype":         phenotype,
            "detected_variants": _build_detected_variants(variants),
        }

        clinical_recommendation = {
            "text":               rec["text"],
            "implication":        rec["implication"],
            "recommendation_url": gene_url,
        }

        llm_generated_explanation = {
            "summary": (
                f"Patient has {phenotype} status for {primary_gene}. "
                f"Diplotype inferred: {diplotype}. "
                f"{rec['implication']}"
            ),
            "supporting_facts": [
                f"{v.gene} variant at {v.chrom}:{v.pos} ({v.ref}→{v.alt}), "
                f"zygosity={v.zygosity}, star={v.star or 'unknown'}"
                for v in variants
            ],
        }

        quality_metrics = {
            **parsed.quality_metrics,
            "vcf_parsing_success": bool(parsed.quality_metrics.get("vcf_parsing_success", True)),
        }

        results.append(
            {
                "patient_id":              parsed.patient_id or "PATIENT_UNKNOWN",
                "drug":                    drug_upper,
                "timestamp":               timestamp,
                "risk_assessment":         risk_assessment,
                "pharmacogenomic_profile": pharmacogenomic_profile,
                "clinical_recommendation": clinical_recommendation,
                "llm_generated_explanation": llm_generated_explanation,
                "quality_metrics":         quality_metrics,
            }
        )

    return results
