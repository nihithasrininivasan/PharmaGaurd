import logging
import asyncio
import uuid
import time
from datetime import datetime, timezone

from fastapi import UploadFile

from app.schemas.internal_contracts import RiskEngineOutput
from app.schemas.pharma_schema import (
    PharmaGuardResponse,
    RiskAssessment,
    PharmacogenomicProfile,
    LLMExplanation,
    ClinicalRecommendation,
    QualityMetrics
)
from app.services.llm.explanation_service import generate_explanation_background, ULTRA_LIGHTNING_CACHE, generate_explanation
from app.services.vcf.pharmaguard_adapter import analyze_vcf_for_drugs

logger = logging.getLogger(__name__)


def compute_heatmap_intensity(severity: str = "low", phenotype: str = "NM") -> int:
    """Computes a 0-4 heatmap intensity score from severity and phenotype."""
    mapping = {
        "none": 0,
        "low": 1,
        "moderate": 2,
        "high": 3,
        "critical": 4
    }

    intensity = mapping.get((severity or "low").lower(), 1)

    # Phenotype boost for extreme metabolizer statuses
    if (phenotype or "NM") in ["PM", "URM", "Poor Metabolizer", "Ultrarapid Metabolizer"]:
        intensity = min(intensity + 1, 4)

    return intensity


async def run_analysis_pipeline(patient_id: str, drug: str, vcf_file: UploadFile) -> PharmaGuardResponse:
    """
    Orchestrates the full pharmacogenomic analysis pipeline:
    1. Parse VCF via real adapter
    2. Extract risk from real risk engine
    3. Generate LLM Explanation
    4. Assemble Response
    """
    logger.info(f"Starting analysis pipeline for patient {patient_id}, drug {drug}")
    start_time = time.time()

    # 1. Parse VCF + compute risk via the real adapter
    logger.info("Parsing VCF and computing risk...")
    vcf_bytes = await vcf_file.read()

    # Validate drug is supported before processing
    from app.services.vcf.pharmaguard_adapter import SUPPORTED_DRUGS_TO_GENE
    drug_upper = drug.strip().upper()
    if drug_upper not in SUPPORTED_DRUGS_TO_GENE:
        supported = ", ".join(sorted(SUPPORTED_DRUGS_TO_GENE.keys()))
        raise ValueError(
            f"Drug '{drug}' is not supported for pharmacogenomic analysis. "
            f"Supported drugs: {supported}"
        )

    results = analyze_vcf_for_drugs(vcf_bytes, [drug_upper])

    if not results:
        target_gene = SUPPORTED_DRUGS_TO_GENE[drug_upper]
        raise ValueError(
            f"The uploaded VCF file does not contain any variants in the {target_gene} gene region, "
            f"which is required for {drug_upper} pharmacogenomic analysis. "
            f"Please upload a VCF file that includes {target_gene} variant data."
        )

    result = results[0]  # Single drug result

    # Extract fields from the real adapter output
    ra = result.get("risk_assessment", {})
    profile = result.get("pharmacogenomic_profile", {})
    rec = result.get("clinical_recommendation", {})

    gene = profile.get("primary_gene", "Unknown")
    diplotype = profile.get("diplotype", "Unknown")
    phenotype = profile.get("phenotype", "Unknown")
    risk_label = ra.get("risk_label", "Unknown")
    severity = ra.get("severity", "moderate")
    confidence_score = ra.get("confidence_score", 0.0)
    recommendation_text = rec.get("text", "Consult clinical pharmacist.")
    raw_variants = profile.get("detected_variants", [])

    # Filter out Hom-Ref (0/0) variants — patient doesn't carry these
    # Also transform to frontend-expected format: { rsid, effect }
    detected_variants = []
    for v in raw_variants:
        if v.get("zygosity") == "Hom-Ref" or v.get("gt") == "0/0":
            continue
        info = v.get("info", {})
        effect_parts = []
        if info.get("FUNC"):
            effect_parts.append(info["FUNC"].replace("_", " ").title())
        if info.get("CLNSIG"):
            effect_parts.append(info["CLNSIG"].replace("_", " "))
        if v.get("star"):
            effect_parts.append(f"Star allele *{v['star'].lstrip('*')}")
        effect = " — ".join(effect_parts) if effect_parts else f"{v.get('gene', '')} variant"
        detected_variants.append({
            "rsid": v.get("rsid", "unknown"),
            "effect": effect,
        })

    # Build RiskEngineOutput for LLM explanation
    risk_data = RiskEngineOutput(
        gene=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        risk_label=risk_label,
        severity=severity,
        recommendation=recommendation_text,
        detected_variants=[{"id": v["rsid"], "effect": v["effect"]} for v in detected_variants],
    )

    # 2. LLM Explanation (Non-blocking with job tracking)
    logger.info("LLM explanation generation started")

    cache_key = f"{gene}:{diplotype}:{drug}".lower()

    if cache_key in ULTRA_LIGHTNING_CACHE:
        explanation_text = ULTRA_LIGHTNING_CACHE[cache_key]
        job_id = None
    else:
        job_id = str(uuid.uuid4())
        asyncio.create_task(generate_explanation_background(job_id, risk_data, drug))
        explanation_text = f"Generating clinical explanation... job_id:{job_id}"

    # 3. Normalize phenotype for frontend display
    phenotype_map = {
        "Poor Metabolizer": "PM",
        "Intermediate Metabolizer": "IM",
        "Normal Metabolizer": "NM",
        "Rapid Metabolizer": "RM",
        "Ultra Rapid Metabolizer": "URM",
        "Ultrarapid Metabolizer": "URM",
        "Normal Function": "NF",
        "Increased Function": "IF",
        "Decreased Function": "DF",
    }
    normalized_phenotype = phenotype_map.get(phenotype, phenotype)

    # 4. Normalize risk label for frontend consistency
    label_map = {
        "Standard dosing recommended": "Safe",
        "Avoid": "Toxic",
        "Use Alternative": "Ineffective",
        "Adjust Dosage": "Adjust Dosage",
        "Toxic": "Toxic",
        "Ineffective": "Ineffective",
        "Safe": "Safe",
    }
    normalized_label = label_map.get(risk_label, "Adjust Dosage")

    current_time = datetime.now(timezone.utc).isoformat()

    response = PharmaGuardResponse(
        patient_id=patient_id,
        drug=drug.upper(),
        timestamp=current_time,
        risk_assessment=RiskAssessment(
            risk_label=normalized_label,
            confidence_score=confidence_score,
            severity=severity,
        ),
        pharmacogenomic_profile=PharmacogenomicProfile(
            primary_gene=gene,
            diplotype=diplotype,
            phenotype=normalized_phenotype,
            detected_variants=detected_variants,
        ),
        clinical_recommendation=ClinicalRecommendation(
            action=recommendation_text,
            source="CPIC Guidelines",
        ),
        llm_generated_explanation=LLMExplanation(
            summary=explanation_text,
        ),
        quality_metrics=QualityMetrics(
            vcf_parsing_success=True,
            coverage_check="PASS",
            extra_metadata={
                "heatmap_intensity": compute_heatmap_intensity(
                    severity, normalized_phenotype
                )
            },
        ),
    )

    total_time = time.time() - start_time
    logger.info(f"⚡ Pipeline execution time: {total_time:.2f} seconds")

    return response
