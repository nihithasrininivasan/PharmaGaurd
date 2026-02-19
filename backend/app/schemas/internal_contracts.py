from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any

class RiskEngineOutput(BaseModel):
    """
    Internal contract between the CPIC Risk Engine (Person 2) and the Backend Analysis layer (Person 3).
    This model encapsulates the raw pharmacogenomic findings before they are enriched by the LLM
    or formatted for the final API response.
    """
    gene: str = Field(..., description="The gene symbol (e.g., CYP2C19)")
    diplotype: str = Field(..., description="The detected diplotype (e.g., *1/*2)")
    phenotype: str = Field(..., description="The metabolizer status (e.g., Intermediate Metabolizer)")
    risk_label: str = Field(..., description="High-level risk category (e.g., High Risk, Moderate Risk)")
    severity: str = Field(..., description="Clinical severity level (e.g., Warning, Contraindication)")
    recommendation: str = Field(..., description="Core clinical recommendation text from CPIC guidelines")
    detected_variants: Optional[List[Dict[str, Any]]] = Field(
        default=None, 
        description="List of specific variants detected (e.g., rs12345: G>A)"
    )
