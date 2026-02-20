from fastapi import APIRouter, HTTPException, Query
from datetime import datetime
from typing import List

from app.services.pharmacogenomics.models import (
    PharmacogenomicReport,
    PharmacogenomicProfile,
    GenotypeData,
    VariantCall
)
from app.services.pharmacogenomics.cpic_loader import get_cpic_loader
from app.services.pharmacogenomics.phenotype_mapper import PhenotypeMapper
from app.services.pharmacogenomics.risk_engine import RiskEngine

router = APIRouter()

@router.get("/report", response_model=PharmacogenomicReport)
async def get_pharmacogenomic_report(
    drug: str = Query(..., description="Drug name (e.g., clopidogrel)"),
    gene: str = Query(..., description="Gene symbol (e.g., CYP2C19)"),
    # In a real app, variants would come from a VCF file or database
    # For this implementation, we allow passing variant details as query params for testing
    # or assume we simulate a patient profile.
    # To match the manual_test.py flow without VCF upload, we'll accept a simplified variant list
    # or just use a dummy patient ID to trigger the logic.
    patient_id: str = Query("PATIENT_001", description="Patient Identifier")
):
    """
    Generate a Pharmacogenomic Risk Report for a specific drug/gene pair.
    
    Current implementation is a direct port of the logic in `manual_test.py`.
    It assumes a user has specific variants. Since GET requests make complex variant lists hard,
    we will simulate the input variants that triggered the specific cases in manual test
    OR we can switch to POST to accept a list of variants.
    
    Let's make it a POST endpoint to accept variants properly.
    """
    pass

@router.post("/report", response_model=PharmacogenomicReport)
async def create_pharmacogenomic_report(
    drug: str,
    gene: str,
    variants: List[VariantCall],
    patient_id: str = "PATIENT_001"
):
    """
    Generate a Pharmacogenomic Risk Report based on provided variants.
    """
    loader = get_cpic_loader()
    
    # 1. Process Genotype
    genotype_data = GenotypeData(
        sample_id=patient_id,
        gene_symbol=gene,
        variants=variants,
        covered_positions=[], 
        coverage_mean=30.0
    )
    
    mapper = PhenotypeMapper()
    diplotype_result = mapper.process_genotype(genotype_data)
    
    # 2. Evaluate Risk
    engine = RiskEngine()
    risk, recommendation = engine.evaluate_risk(
        drug=drug,
        gene=gene,
        phenotype=diplotype_result.phenotype,
        diplotype=diplotype_result.diplotype,
        diplotype_confidence=diplotype_result.confidence
    )
    
    # 3. Construct Report
    return PharmacogenomicReport(
        patient_id=patient_id,
        drug=drug,
        timestamp=datetime.utcnow().isoformat() + "Z",
        risk_assessment=risk,
        pharmacogenomic_profile=PharmacogenomicProfile(
            primary_gene=gene,
            diplotype=diplotype_result.diplotype,
            phenotype=diplotype_result.phenotype,
            detected_variants=variants
        ),
        clinical_recommendation=recommendation
    )
