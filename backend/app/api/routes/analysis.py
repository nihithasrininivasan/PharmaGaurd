from fastapi import APIRouter, UploadFile, File, Form, HTTPException, status
import logging

from app.schemas.pharma_schema import PharmaGuardResponse
from app.services.pipeline.analysis_pipeline import run_analysis_pipeline
from app.services.llm.explanation_service import EXPLANATION_STORE

router = APIRouter()
logger = logging.getLogger(__name__)

@router.post(
    "/analyze",
    response_model=PharmaGuardResponse,
    status_code=status.HTTP_200_OK,
    summary="Analyze Pharmacogenomic Risk",
    description="Upload a VCF file and specify a drug to receive a comprehensive pharmacogenomic risk assessment."
)
async def analyze_pharmacogenomics(
    drug: str = Form(..., description="The name of the drug to analyze (e.g., Clopidogrel)"),
    vcf: UploadFile = File(..., description="Patient's VCF file containing genetic variants"),
    patient_id: str = Form("anonymous", description="Optional patient identifier"),
) -> PharmaGuardResponse:
    """
    Endpoint to trigger the pharmacogenomic analysis pipeline.
    
    - **drug**: Target drug name
    - **vcf**: Genetic data file
    - **patient_id**: Optional identifier
    """
    if not vcf.filename.endswith(('.vcf', '.vcf.gz')):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid file format. Please upload a .vcf or .vcf.gz file."
        )

    try:
        response = await run_analysis_pipeline(patient_id, drug, vcf)
        return response
        
    except ValueError as ve:
        logger.error(f"Validation error in pipeline: {str(ve)}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(ve)
        )
    except Exception as e:
        logger.exception(f"Unexpected error in analysis pipeline: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="An error occurred during the analysis pipeline."
        )


@router.get("/explanation/{job_id}")
async def get_explanation(job_id: str):
    """Poll for async LLM explanation result by job_id."""
    if job_id in EXPLANATION_STORE:
        return {"summary": EXPLANATION_STORE[job_id]}
    return {"summary": None}
