from fastapi import APIRouter, UploadFile, File, Form, HTTPException, status
from pydantic import BaseModel
import logging

from app.schemas.pharma_schema import PharmaGuardResponse
from app.services.pipeline.analysis_pipeline import run_analysis_pipeline
from app.services.llm.explanation_service import EXPLANATION_STORE, is_supported_gene_drug
from app.services.llm.ollama_client import OllamaClient

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

    # STRICT VALIDATION: Check for valid VCF header
    # We peek at the first chunk without fully consuming if possible, or read/reset
    # But since run_analysis_pipeline reads it, we just need to ensure it's valid first.
    # For UploadFile, we can read, check, and seek(0).
    header_chunk = await vcf.read(1024)
    await vcf.seek(0)
    
    try:
        header_str = header_chunk.decode('utf-8', errors='ignore')
    except:
        raise HTTPException(status_code=400, detail="File is not a valid text-based VCF.")

    if not header_str.startswith("##fileformat=VCF"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid VCF: File must start with '##fileformat=VCF'. Please upload a verified ISGCR/VCF file."
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


class AskRequest(BaseModel):
    question: str
    gene: str = ""
    diplotype: str = ""
    phenotype: str = ""
    drug: str = ""


@router.post("/ask")
async def ask_pharmaguard(req: AskRequest):
    """Clinical AI chatbot — calls LLM directly (no cache) with clinical grounding."""

    # Gene-drug guardrail
    gene_note = ""
    if req.gene and req.drug and not is_supported_gene_drug(req.gene, req.drug):
        gene_note = "\nNOTE: Gene is not primary metabolism pathway for this drug."

    prompt = f"""SYSTEM:
You are a pharmacogenomics clinical assistant.
Follow CPIC guidance strictly.

Rules:
- Only use provided context.
- Do NOT invent biology.
- If gene is not primary for drug, say so clearly.
- Maximum 2 sentences.
- Maximum 45 words.
- Clinician tone only.{gene_note}

Gene={req.gene}
Diplotype={req.diplotype}
Phenotype={req.phenotype}
Drug={req.drug}

Question: {req.question}

Provide concise CPIC-grounded clinical interpretation.

ANSWER:"""

    try:
        client = OllamaClient()
        answer = await client.generate_chat_text(prompt)

        if not answer:
            answer = "I'm unable to generate a response right now. Please consult CPIC guidelines or your clinical pharmacist."

        return {"answer": answer.strip()}

    except Exception as e:
        logger.error(f"Ask PharmaGuard error: {str(e)}")
        return {"answer": "Clinical AI is temporarily unavailable. Please try again."}

