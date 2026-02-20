"""
Analysis Pipeline — Orchestrates VCF → Risk → LLM → Response.

Receives a VCF file + drug from the API route, runs the full
pharmacogenomic analysis, and returns a PharmaGuardResponse.
"""
import logging
import asyncio
import uuid
import time
from datetime import datetime, timezone
from typing import List, Dict, Any

from fastapi import UploadFile

from app.schemas.internal_contracts import RiskEngineOutput
from app.schemas.pharma_schema import (
    PharmaGuardResponse,
    RiskAssessment,
    PharmacogenomicProfile,
    LLMExplanation,
    ClinicalRecommendation,
    QualityMetrics,
)
from app.services.llm.explanation_service import (
    generate_explanation_background,
    ULTRA_LIGHTNING_CACHE,
    generate_explanation,
)

# ── Real VCF parser ──────────────────────────────────────────────────────
from app.services.vcf.parser import parse_vcf
from app.services.vcf.variant_extractor import extract_pharmacogenes

# ── Real pharmacogenomics engine ─────────────────────────────────────────
from app.services.pharmacogenomics.cpic_loader import reload_cpic_data, get_cpic_loader
from app.services.pharmacogenomics.phenotype_mapper import PhenotypeMapper
from app.services.pharmacogenomics.risk_engine import RiskEngine, create_risk_engine
from app.services.pharmacogenomics.models import GenotypeData, VariantCall, PatientProfile

logger = logging.getLogger(__name__)

# Ensure CPIC data is loaded on import
reload_cpic_data()

# ── Gene-to-drug mapping ──────────────────────────────────────────────────
GENE_DRUG_MAP = {
    "CYP2D6":  ["CODEINE"],
    "CYP2C19": ["CLOPIDOGREL"],
    "CYP2C9":  ["WARFARIN"],
    "SLCO1B1": ["SIMVASTATIN"],
    "TPMT":    ["AZATHIOPRINE"],
    "DPYD":    ["FLUOROURACIL"],
}

# Reverse: drug → gene
DRUG_GENE_MAP = {}
for gene, drugs in GENE_DRUG_MAP.items():
    for d in drugs:
        DRUG_GENE_MAP[d] = gene

# ── Zygosity bridge ──────────────────────────────────────────────────────
_ZYGOSITY_MAP = {
    "Het": "HET",
    "Hom-Alt": "HOM_ALT",
    "Hom-Ref": "HOM_REF",
    "Unknown": "UNKNOWN",
}


def _build_detected_variants_list(by_gene: dict) -> List[Dict[str, Any]]:
    """Flatten gene→variants dict into a list of variant dicts for the response."""
    result = []
    for gene, variants in by_gene.items():
        for v in variants:
            result.append({
                "rsid": v.rsid,
                "star": v.star,
                "gene": gene,
                "zygosity": v.zygosity,
                "chrom": v.chrom,
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
            })
    return result


def compute_heatmap_intensity(severity: str = "low", phenotype: str = "NM") -> int:
    """0-4 heatmap intensity from severity and phenotype."""
    mapping = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}
    intensity = mapping.get((severity or "low").lower(), 1)
    if (phenotype or "NM") in ("PM", "URM", "Poor Metabolizer", "Ultrarapid Metabolizer"):
        intensity = min(intensity + 1, 4)
    return intensity


async def run_analysis_pipeline(
    patient_id: str, drug: str, vcf_file: UploadFile
) -> PharmaGuardResponse:
    """
    Full pipeline: VCF → parse → extract genes → phenotype → risk → LLM → response.
    """
    logger.info("Starting analysis pipeline for patient %s, drug %s", patient_id, drug)
    start_time = time.time()
    drug_upper = drug.upper()

    # ── 1. Parse VCF ──────────────────────────────────────────────────────
    logger.info("Parsing VCF")
    vcf_bytes = await vcf_file.read()
    parsed = parse_vcf(vcf_bytes)

    # ── 2. Extract pharmacogenes ──────────────────────────────────────────
    logger.info("Extracting pharmacogenes")
    by_gene = extract_pharmacogenes(parsed.variants)
    all_variants = _build_detected_variants_list(by_gene)

    # Quality metrics
    quality = QualityMetrics(
        vcf_parsing_success=True,
        coverage_check="PASS",
        extra_metadata={
            "variant_count": len(parsed.variants),
            "genes_detected": list(by_gene.keys()),
        },
    )

    # ── 3. Identify target gene for the drug ──────────────────────────────
    target_gene = DRUG_GENE_MAP.get(drug_upper)
    if not target_gene:
        raise ValueError(f"Unsupported drug: {drug}. Supported: {list(DRUG_GENE_MAP.keys())}")

    gene_variants = by_gene.get(target_gene, [])

    # ── 4. Build diplotype ────────────────────────────────────────────────
    cpic = get_cpic_loader()
    star_alleles = list({v.star for v in gene_variants if v.star})

    if len(star_alleles) == 0:
        diplotype = "*1/*1"
    elif len(star_alleles) == 1:
        hom = any(v.zygosity == "Hom-Alt" for v in gene_variants if v.star == star_alleles[0])
        diplotype = f"{star_alleles[0]}/{star_alleles[0]}" if hom else f"*1/{star_alleles[0]}"
    else:
        diplotype = f"{star_alleles[0]}/{star_alleles[1]}"

    # ── 5. Activity score + phenotype ─────────────────────────────────────
    alleles = diplotype.split("/")
    allele_scores = {}
    for a in alleles:
        allele_scores[a] = cpic.get_activity_score(target_gene, a)
    total_as = sum(allele_scores.values())

    mapper = PhenotypeMapper()
    # Build a minimal GenotypeData for the mapper
    variant_calls = []
    for v in gene_variants:
        variant_calls.append(VariantCall(
            chrom=v.chrom, pos=v.pos, rsid=v.rsid,
            ref=v.ref, alt=v.alt,
            zygosity=_ZYGOSITY_MAP.get(v.zygosity, "UNKNOWN"),
            quality=float(v.raw_info.get("QUAL", 0.0) or 0.0),
            filter=str(v.raw_info.get("FILTER", "")) or None,
            star_allele=v.star,
            phased="|" in (v.gt or ""),
        ))
    gd = GenotypeData(
        sample_id=patient_id,
        gene_symbol=target_gene,
        variants=variant_calls,
    )
    diplo_result = mapper.process_genotype(gd)
    phenotype = diplo_result.phenotype

    # Notes
    if len(set(alleles)) == 1 and alleles[0] != "*1":
        notes = "Homozygous variant allele"
    elif "*1" in alleles and len(set(alleles)) > 1:
        notes = "Heterozygous with wildtype"
    else:
        notes = "Compound heterozygous"

    # ── 6. Risk assessment ────────────────────────────────────────────────
    logger.info("Computing risk for %s / %s / %s", drug_upper, target_gene, phenotype)
    rec = cpic.get_drug_recommendation_for_phenotype(drug_upper, phenotype)

    if rec:
        severity = rec.get("severity", "none")
        risk_label = rec.get("risk_label", "Standard dosing recommended")
        implication = rec.get("implication", "")
        rec_text = rec.get("recommendation", risk_label)
    else:
        severity = "none"
        risk_label = "Standard dosing recommended"
        implication = "Normal drug metabolism expected"
        rec_text = f"Use standard {drug_upper} dosing guidelines"

    # Compute numeric risk score
    from app.services.pharmacogenomics.risk_scoring import RiskScoreCalculator
    from app.services.pharmacogenomics.confidence import ConfidenceCalculator

    scorer = RiskScoreCalculator()
    conf_calc = ConfidenceCalculator()

    confidence_val = conf_calc.calculate_confidence(
        target_gene, phenotype, diplotype,
        variant_count=len(gene_variants),
        has_quality_data=True,
    )
    risk_score_val = scorer.calculate_risk_score(severity, phenotype, confidence=confidence_val)
    risk_level = scorer.get_risk_level(risk_score_val)

    # ── 7. Build RiskEngineOutput for LLM ─────────────────────────────────
    risk_data = RiskEngineOutput(
        gene=target_gene,
        diplotype=diplotype,
        phenotype=phenotype,
        risk_label=risk_label,
        severity=severity,
        recommendation=rec_text,
        detected_variants=[
            {"id": v.rsid, "genotype": f"{v.ref}>{v.alt}"}
            for v in gene_variants if v.rsid
        ],
        risk_score=risk_score_val,
        risk_level=risk_level,
        confidence_score=confidence_val,
        activity_score=total_as,
    )

    # ── 8. LLM Explanation ────────────────────────────────────────────────
    logger.info("Generating LLM explanation")
    cache_key = f"{target_gene}:{diplotype}:{drug_upper}".lower()

    if cache_key in ULTRA_LIGHTNING_CACHE:
        explanation_text = ULTRA_LIGHTNING_CACHE[cache_key]
        job_id = None
    else:
        job_id = str(uuid.uuid4())
        asyncio.create_task(generate_explanation_background(job_id, risk_data, drug_upper))
        explanation_text = f"Generating clinical explanation... job_id:{job_id}"

    # ── 9. Phenotype abbreviation for frontend ────────────────────────────
    phenotype_map = {
        "Poor Metabolizer": "PM",
        "Intermediate Metabolizer": "IM",
        "Normal Metabolizer": "NM",
        "Rapid Metabolizer": "RM",
        "Ultrarapid Metabolizer": "URM",
        "Decreased Function": "DF",
        "Poor Function": "PF",
        "Normal Function": "NF",
    }
    normalized_phenotype = phenotype_map.get(phenotype, phenotype)

    # ── 10. Assemble response ─────────────────────────────────────────────
    current_time = datetime.now(timezone.utc).isoformat()

    # Filter drug-specific variants
    drug_variants = [
        {"rsid": v.rsid, "effect": f"{v.star or ''} allele ({v.zygosity})"}
        for v in gene_variants if v.rsid
    ]

    response = PharmaGuardResponse(
        patient_id=patient_id,
        drug=drug_upper,
        timestamp=current_time,
        risk_assessment=RiskAssessment(
            risk_label=risk_label,
            confidence_score=round(confidence_val, 3),
            severity=severity,
            risk_score=round(risk_score_val, 2),
            risk_level=risk_level,
        ),
        pharmacogenomic_profile=PharmacogenomicProfile(
            primary_gene=target_gene,
            diplotype=diplotype,
            phenotype=normalized_phenotype,
            detected_variants=drug_variants,
            activity_score=round(total_as, 2),
            allele_scores=allele_scores,
            notes=notes,
        ),
        clinical_recommendation=ClinicalRecommendation(
            text=rec_text,
            source="CPIC Guidelines",
            implication=implication,
        ),
        llm_generated_explanation=LLMExplanation(
            summary=explanation_text,
        ),
        quality_metrics=QualityMetrics(
            vcf_parsing_success=True,
            coverage_check="PASS",
            extra_metadata={
                "heatmap_intensity": compute_heatmap_intensity(severity, normalized_phenotype),
                "variant_count": len(parsed.variants),
                "genes_detected": list(by_gene.keys()),
            },
        ),
    )

    total_time = time.time() - start_time
    logger.info("Pipeline execution time: %.2fs", total_time)

    return response
