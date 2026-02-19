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

    # Build clinical context (source of truth)
    clinical_context = f"""PRIMARY_GENE: {req.gene}
DIPLOTYPE: {req.diplotype}
PHENOTYPE: {req.phenotype}
DRUG: {req.drug}"""

    # Gene-drug evidence check
    evidence_note = ""
    if req.gene and req.drug and not is_supported_gene_drug(req.gene, req.drug):
        evidence_note = "\nClinical note: The relationship between this gene and drug may be indirect or low-evidence. State this uncertainty explicitly."

    prompt = f"""SYSTEM:
You are a pharmacogenomics clinical reasoning assistant.

STRICT RULES:
- Only explain mechanisms directly supported by PRIMARY_GENE and DRUG below.
- If the gene-drug relationship is weak or uncertain, explicitly say:
  'Evidence linking this gene to this drug is limited.'
- NEVER invent alternative genes not listed in the context.
- NEVER introduce new variants not in the context.
- NEVER contradict phenotype meaning:
    RM/URM → faster metabolism
    PM/IM → slower metabolism
    NM → normal metabolism
- Prefer conservative CPIC-style language.
- Use phrases: 'may influence', 'is associated with', 'based on CPIC guidance'.
- Do NOT provide direct medical advice.
- Keep your answer under 3 sentences, maximum 80 words.{evidence_note}

CLINICAL CONTEXT (SOURCE OF TRUTH):
{clinical_context}

USER QUESTION:
{req.question}

ANSWER:"""

    try:
        client = OllamaClient()
        answer = await client.generate_text(prompt)

        if not answer:
            answer = "I'm unable to generate a response right now. Please consult CPIC guidelines or your clinical pharmacist."

        return {"answer": answer.strip()}

    except Exception as e:
        logger.error(f"Ask PharmaGuard error: {str(e)}")
        return {"answer": "Clinical AI is temporarily unavailable. Please try again."}

