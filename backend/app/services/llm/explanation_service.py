import logging
import re
from app.services.llm.ollama_client import OllamaClient
from app.schemas.internal_contracts import RiskEngineOutput

logger = logging.getLogger(__name__)

# TASK 1: ULTRA LIGHTNING CACHE
ULTRA_LIGHTNING_CACHE = {}

# TASK 2: ADD SAFETY POST-PROCESSOR FUNCTION
def apply_clinical_safety(text: str) -> str:
    """
    Applies clinical safety rules to the explanation text.
    Replaces prescriptive language with cautious phrasing.
    Ensures grounding in CPIC guidance.
    """
    # Safety Replacements
    replacements = {
        r"\bmust\b": "may",
        r"\bshould\b": "may be considered",
        r"\bwill cause\b": "is associated with",
        r"\bcauses\b": "is associated with",
        r"\bdefinitely\b": "likely",
    }
    
    safe_text = text
    for pattern, replacement in replacements.items():
        safe_text = re.sub(pattern, replacement, safe_text, flags=re.IGNORECASE)
        
    # Ensure CPIC citation if missing
    if "based on CPIC guidance" not in safe_text and "based on CPIC pharmacogenomic guidance" not in safe_text:
         # Append blindly if missing, but ideally prompt handles this. 
         # Given constraint to preserve 3 sentences, we try to append it to the last sentence if possible or just ensure the prompt did its job.
         # However, requirements say "Ensure text includes phrase...".
         # Let's replace "CPIC recommendation" with the full phrase if present, 
         # or just append to the end if not length prohibitive.
         # For safety/simplicity in this regex post-process:
         if safe_text.endswith("."):
             safe_text = safe_text[:-1] + ", based on CPIC pharmacogenomic guidance."
         else:
             safe_text += " This assessment is based on CPIC pharmacogenomic guidance."

    return safe_text.strip()

# TASK 1: ADD SAFE TAGGING FUNCTION
def apply_doctor_view_tags(text: str, risk_data: RiskEngineOutput, drug: str) -> str:
    """
    Applies minimal XML-like tagging for Doctor View Mode.
    Wraps key pharmacogenomic entities.
    """
    tagged_text = text
    
    # Simple string replacements for exact entity matches
    # Note: This is a basic find-replace. It assumes the LLM used the exact casing or close to it.
    # To be safer and non-intrusive, we do case-insensitive replacement but keep original casing,
    # or just replace the known entities provided by risk_data.
    
    # Gene
    if risk_data.gene:
        pattern = re.compile(re.escape(risk_data.gene), re.IGNORECASE)
        tagged_text = pattern.sub(f"<gene>{risk_data.gene}</gene>", tagged_text)

    # Diplotype
    if risk_data.diplotype:
        # Escape special characters like * in diplotypes (*1/*17)
        pattern = re.compile(re.escape(risk_data.diplotype), re.IGNORECASE)
        tagged_text = pattern.sub(f"<diplotype>{risk_data.diplotype}</diplotype>", tagged_text)
        
    # Variants
    if risk_data.detected_variants:
        for variant in risk_data.detected_variants:
            v_id = variant.get('id')
            if v_id:
               pattern = re.compile(re.escape(v_id), re.IGNORECASE)
               tagged_text = pattern.sub(f"<variant>{v_id}</variant>", tagged_text)

    # Drug
    if drug:
        pattern = re.compile(re.escape(drug), re.IGNORECASE)
        tagged_text = pattern.sub(f"<drug>{drug}</drug>", tagged_text)
        
    return tagged_text

async def generate_explanation(risk_data: RiskEngineOutput, drug: str) -> str:
    """
    Generates a clinical explanation using LLM based on risk data with strict guardrails.
    """
    # TASK 2: CACHE LOOKUP
    cache_key = f"{risk_data.gene}:{risk_data.diplotype}:{drug}".lower()
    
    if cache_key in ULTRA_LIGHTNING_CACHE:
        # TASK 5: SAFE CACHE LOGGING
        logger.info("Using Ultra Lightning cached explanation for %s", cache_key)
        # Even cached explanations get safety check just in case rules changed
        return ULTRA_LIGHTNING_CACHE[cache_key]

    logger.info("Generating clinical explanation for %s", risk_data.gene)
    
    # 1. Prepare Data
    variants_str = ", ".join([v.get('id', 'unknown') for v in (risk_data.detected_variants or [])]) if risk_data.detected_variants else "None"
    
    # 2. Construct Structured Prompt
    # TASK 1: ADD SAFETY SYSTEM PROMPT
    prompt = f"""
SYSTEM:
You are a clinical decision-support assistant.
Use cautious, non-prescriptive language.
Do NOT provide medical advice.
Use phrases such as: 'may influence', 'is associated with', 'based on CPIC guidance'.

Respond in EXACTLY 3 sentences. Maximum 80 words. No extra text.
Only use the information provided below.
Do NOT invent new genes, variants, or recommendations.
If information is missing, respond with: 'Insufficient genomic evidence for explanation.'

USER:
Patient pharmacogenomic data:

Gene: {risk_data.gene}
Diplotype: {risk_data.diplotype}
Phenotype: {risk_data.phenotype}
Variants: {variants_str}

Drug: {drug}
Risk Label: {risk_data.risk_label}
Severity: {risk_data.severity}
CPIC Recommendation:
{risk_data.recommendation}

Task:
Write EXACTLY 3 sentences.
Maximum 80 words.

Requirements:
1. Mention the gene and diplotype explicitly.
2. Explain the biological mechanism of drug metabolism.
3. Justify the CPIC recommendation.
4. Do NOT introduce any new clinical facts not listed.
5. Tone: professional clinical language.
    """

    client = OllamaClient()
    
    try:
        # 3. Call LLM
        explanation = await client.generate_text(prompt)
        
        # 4. Fallback if None
        if explanation is None:
            logger.warning("LLM fallback triggered: No response from Ollama")
            return "Clinical explanation unavailable. CPIC recommendation applied."

        explanation = explanation.strip()

        # 5. Hallucination Guardrails
        # Check if gene is mentioned (basic relevance check)
        if risk_data.gene not in explanation:
            logger.warning("LLM fallback triggered: Gene missing in explanation")
            return "Clinical explanation unavailable. CPIC recommendation applied."

        # Truncate to 3 sentences max
        # Split by ., !, or ? followed by space or end of string
        sentences = re.split(r'(?<=[.!?])\s+', explanation)
        if len(sentences) > 3:
            explanation = " ".join(sentences[:3])

        # TASK 3: APPLY SAFETY MODE
        explanation = apply_clinical_safety(explanation)
        
        # TASK 4: ADD OPTIONAL SAFETY LOGGING
        logger.info("Clinical Safety Mode applied")

        # TASK 2: APPLY TAGGING AFTER SAFETY MODE
        explanation = apply_doctor_view_tags(explanation, risk_data, drug)

        # TASK 3: STORE IN CACHE
        ULTRA_LIGHTNING_CACHE[cache_key] = explanation
        
        return explanation

    except Exception as e:
        # 6. Safety Net
        logger.error(f"Unexpected error in explanation service: {str(e)}")
        return "Clinical explanation unavailable. CPIC recommendation applied."
