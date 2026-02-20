from __future__ import annotations

import datetime as _dt
import logging
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

from .parser import VcfParseResult, parse_vcf
from .variant_extractor import ExtractedVariant, extract_pharmacogenes

# Import core engine components
from app.services.pharmacogenomics.phenotype_mapper import PhenotypeMapper
from app.services.pharmacogenomics.risk_engine import RiskEngine
from app.services.pharmacogenomics.models import (
    VariantCall, GenotypeData, RiskAssessment, ClinicalRecommendation,
)
from app.services.pharmacogenomics.variant_normalizer import normalize_variants
from app.services.pharmacogenomics.confidence import ConfidenceCalculator

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Drug → primary gene mapping
# ---------------------------------------------------------------------------
SUPPORTED_DRUGS_TO_GENE: Mapping[str, str] = {
    "CODEINE":       "CYP2D6",
    "WARFARIN":      "CYP2C9",
    "CLOPIDOGREL":   "CYP2C19",
    "SIMVASTATIN":   "SLCO1B1",
    "AZATHIOPRINE":  "TPMT",
    "THIOGUANINE":   "TPMT",
    "FLUOROURACIL":  "DPYD",
    "5-FLUOROURACIL":"DPYD",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise_drug_name(name: str) -> str:
    return name.strip().upper()


def _map_zygosity(vcf_zygosity: str) -> Optional[str]:
    """
    Map VCF parser zygosity strings to Pharmacogenomics Engine model.
    Returns None for Unknown zygosity — caller must decide how to handle.
    """
    mapping = {
        "Hom-Ref": "HOM_REF",
        "Het": "HET",
        "Hom-Alt": "HOM_ALT",
    }
    return mapping.get(vcf_zygosity)  # Returns None for "Unknown" or unmapped


def _extract_allele_depth(sample_info: Dict) -> Optional[List[int]]:
    """Extract allele depth (AD) from VCF sample fields."""
    ad_str = sample_info.get("AD", "")
    if not ad_str or ad_str == ".":
        return None
    try:
        return [int(x) for x in ad_str.split(",")]
    except (ValueError, AttributeError):
        return None


def _detect_phasing(gt: Optional[str]) -> bool:
    """Detect phasing from GT field. '|' separator = phased."""
    if not gt:
        return False
    return "|" in gt


def _build_detected_variants(variants: Sequence[ExtractedVariant]) -> List[Dict]:
    """Preserve the original output structure for detected variants."""
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
    # 1. Parse VCF
    parsed: VcfParseResult = parse_vcf(file_bytes)
    by_gene = extract_pharmacogenes(parsed.variants)

    results: List[Dict] = []
    timestamp = _now_iso()
    
    # Initialize Engine Components
    mapper = PhenotypeMapper()
    engine = RiskEngine()
    confidence_calc = ConfidenceCalculator()

    for raw_drug in drugs:
        drug_upper = _normalise_drug_name(raw_drug)
        primary_gene = SUPPORTED_DRUGS_TO_GENE.get(drug_upper, "")
        
        # Get raw extracted variants for this gene
        extracted_variants: Sequence[ExtractedVariant] = by_gene.get(primary_gene, [])

        # 2. Map to Domain Models — propagate REAL quality data
        domain_variants: List[VariantCall] = []
        rejected_unknown_zyg: List[Dict] = []

        # Check Genome Build
        if parsed.genome_build == "GRCh37":
            # Critical warning: mismatch likely
            # We could reject entirely, but for now we flag it heavily.
            # In a real strict mode, we might return an empty result with error.
            logger.warning("Detected GRCh37/hg19 genome build. System expects GRCh38.")
            # For this implementation plan, we proceed but flag it. 
            # Or better: if we want to be SAFE, we should probably stop if we rely on positions.
            # But specific gene extractors might rely on gene names/RSIDs.
            pass

        for v in extracted_variants:
            # --- Hardening: Strict QC ---
            
            # 1. FILTER check
            # Accept PASS, ".", or None. Reject everything else.
            if v.filter and v.filter not in ("PASS", ".", ""):
                rejected_unknown_zyg.append({
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "reason": f"FILTER failed: {v.filter}",
                })
                continue

            # 2. QUAL check
            # Reject if raw QUAL < 20
            raw_qual = v.raw_info.get("_QUAL")
            quality = float(raw_qual) if raw_qual is not None else 0.0
            if raw_qual is not None and quality < 20.0:
                 rejected_unknown_zyg.append({
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "reason": f"Low QUAL: {quality}",
                })
                 continue

            # 3. Depth check
            # Reject if DP < 10 (calculated from AD sum)
            ad = _extract_allele_depth(v.raw_info)
            if ad and sum(ad) < 10:
                 rejected_unknown_zyg.append({
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "reason": f"Low Depth (DP < 10): {sum(ad)}",
                })
                 continue

            # Map zygosity — reject Unknown instead of guessing HET
            mapped_zyg = _map_zygosity(v.zygosity)
            if mapped_zyg is None:
                # Unknown zygosity → reject, do not guess
                rejected_unknown_zyg.append({
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "reason": f"Unknown zygosity: {v.zygosity}",
                })
                logger.warning(
                    f"Rejected variant at {v.chrom}:{v.pos} — unknown zygosity '{v.zygosity}'"
                )
                continue

            # Filter out HOM_REF (0/0) — patient does NOT carry the ALT allele.
            # These are evidence of wild-type at this position, not variant evidence.
            # Keeping them would prevent the resolver from entering the wildtype path.
            if mapped_zyg == "HOM_REF":
                continue

            raw_filter = v.raw_info.get("_FILTER")
            filt = str(raw_filter) if raw_filter else None

            # Detect phasing from GT
            phased = _detect_phasing(v.gt)

            domain_variants.append(
                VariantCall(
                    chrom=v.chrom,
                    pos=v.pos,
                    ref=v.ref,
                    alt=v.alt,
                    rsid=v.rsid or None,
                    star_allele=v.star or None,
                    zygosity=mapped_zyg,
                    quality=quality,
                    filter=filt,
                    ad=ad,
                    phased=phased,
                )
            )

        # 2b. Run normalisation pipeline on domain variants
        norm_result = normalize_variants(domain_variants)
        clean_variants = norm_result.clean_variants

        # 3. Process Phenotype (Diplotype Resolution)
        genotype_data = GenotypeData(
            sample_id=parsed.patient_id or "PATIENT_UNKNOWN",
            gene_symbol=primary_gene,
            variants=clean_variants,
            covered_positions=[],  # Not available from simple VCF parse
            coverage_mean=30.0,
            genome_build=parsed.genome_build or None,
        )
        
        if primary_gene:
            diplotype_result = mapper.process_genotype(genotype_data)
            diplotype = diplotype_result.diplotype
            phenotype = diplotype_result.phenotype
            confidence = diplotype_result.confidence
            confidence_breakdown = diplotype_result.confidence_breakdown
        else:
            diplotype = "Unknown"
            phenotype = "Unknown"
            confidence = 0.0
            confidence_breakdown = None

        # 4. Evaluate Risk (CPIC Guidelines)
        if primary_gene:
            risk, recommendation = engine.evaluate_risk(
                drug=drug_upper,
                gene=primary_gene,
                phenotype=phenotype,
                diplotype=diplotype,
                diplotype_confidence=confidence,
                diplotype_confidence_breakdown=confidence_breakdown,
            )
        else:
            risk = RiskAssessment(
                risk_label="Unknown Drug/Gene",
                confidence_score=0.0,
                severity="none",
                confidence_breakdown=None,
            )
            recommendation = ClinicalRecommendation(
                text="Drug not supported.",
                implication="None",
                recommendation_url=None,
            )

        # 5. Construct Output
        risk_assessment = {
            "risk_label":       risk.risk_label,
            "confidence_score": risk.confidence_score,
            "severity":         risk.severity,
            "confidence_breakdown": risk.confidence_breakdown,
        }

        # Attach PharmGKB data layers if present on the RiskAssessment
        if hasattr(risk, 'gene_drug_confirmation') and risk.gene_drug_confirmation:
            risk_assessment["gene_drug_confirmation"] = risk.gene_drug_confirmation
        if hasattr(risk, 'evidence_level') and risk.evidence_level:
            risk_assessment["evidence_level"] = risk.evidence_level
        if hasattr(risk, 'clinical_annotations') and risk.clinical_annotations:
            risk_assessment["clinical_annotations"] = risk.clinical_annotations

        # Automation status: single location on RiskAssessment model only
        if hasattr(risk, 'automation_status') and risk.automation_status:
            risk_assessment["automation_status"] = risk.automation_status

        # Variant-level annotations from PharmGKB
        variant_annotations = []
        if extracted_variants and primary_gene:
            try:
                from ..pharmacogenomics.pharmgkb_loader import get_pharmgkb_loader
                pgkb = get_pharmgkb_loader()
                # Collect rsIDs and star alleles for lookup
                var_ids = []
                for v in extracted_variants:
                    if v.rsid:
                        var_ids.append(v.rsid)
                    if v.star:
                        var_ids.append(v.star)
                if var_ids:
                    variant_annotations = pgkb.get_variant_annotations(primary_gene, var_ids)
            except Exception:
                pass  # Graceful degradation if loader unavailable

        pharmacogenomic_profile = {
            "primary_gene":      primary_gene or "Unknown",
            "diplotype":         diplotype,
            "phenotype":         phenotype,
            "detected_variants": _build_detected_variants(extracted_variants),
            "variant_annotations": variant_annotations if variant_annotations else (
                [] if extracted_variants else
                [{"note": "No actionable variants detected in required gene loci"}]
            ),
        }

        clinical_recommendation = {
            "text":               recommendation.text,
            "implication":        recommendation.implication,
            "recommendation_url": recommendation.recommendation_url,
        }

        llm_generated_explanation = {
            "summary": (
                f"Patient has {phenotype} status for {primary_gene}. "
                f"Diplotype inferred: {diplotype}. "
                f"{recommendation.implication}"
            ),
            "supporting_facts": [
                f"{v.gene} variant at {v.chrom}:{v.pos} ({v.ref}→{v.alt}), "
                f"zygosity={v.zygosity}, star={v.star or 'unknown'}"
                for v in extracted_variants
            ],
        }

        quality_metrics = {
            **parsed.quality_metrics,
            "vcf_parsing_success": bool(parsed.quality_metrics.get("vcf_parsing_success", True)),
            "variants_rejected_unknown_zygosity": len(rejected_unknown_zyg),
            "variants_rejected_quality": len(norm_result.rejected_variants),
            "variants_after_normalization": len(clean_variants),
            "duplicates_removed": norm_result.duplicates_removed,
            "chromosomes_normalized": norm_result.chromosome_normalized,
        }

        if norm_result.build_validation and not norm_result.build_validation.is_valid:
            quality_metrics["genome_build_warning"] = norm_result.build_validation.warnings

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
