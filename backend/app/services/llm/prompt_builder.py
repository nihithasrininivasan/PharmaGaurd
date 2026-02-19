def build_prompt(gene: str, diplotype: str, phenotype: str, drug: str, recommendation: str) -> str:
    """
    Constructs a prompt for the LLM to generate a clinical explanation.
    
    Args:
        gene: The gene symbol.
        diplotype: The detected diplotype.
        phenotype: The metabolizer status.
        drug: The drug name.
        recommendation: The clinical recommendation text.
        
    Returns:
        A formatted prompt string.
    """
    prompt = (
        f"Act as a clinical pharmacogenomics expert. "
        f"Patient data: Gene: {gene}, Diplotype: {diplotype}, Phenotype: {phenotype}, Drug: {drug}. "
        f"Clinical Recommendation: {recommendation}. "
        f"Write a concise 3 sentence clinical pharmacogenomic explanation mentioning enzyme metabolism and drug response."
    )
    return prompt
