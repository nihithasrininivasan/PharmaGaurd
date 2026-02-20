from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any

class RiskEngineOutput(BaseModel):
    """
    Internal contract between the CPIC Risk Engine and the Backend Analysis layer.
    Encapsulates raw pharmacogenomic findings before LLM enrichment.
    """
    gene: str = Field(..., description="Gene symbol (e.g., CYP2C19)")
    diplotype: str = Field(..., description="Detected diplotype (e.g., *1/*2)")
    phenotype: str = Field(..., description="Metabolizer status (e.g., Intermediate Metabolizer)")
    risk_label: str = Field(..., description="High-level risk category")
    severity: str = Field(..., description="Clinical severity level")
    recommendation: str = Field(..., description="Core clinical recommendation text")
    detected_variants: Optional[List[Dict[str, Any]]] = Field(
        default=None,
        description="List of detected variants (e.g., rs12345: G>A)"
    )
    risk_score: Optional[float] = Field(default=None, description="Numeric risk score 0-100")
    risk_level: Optional[str] = Field(default=None, description="Risk level classification")
    confidence_score: Optional[float] = Field(default=None, description="Confidence 0.0-1.0")
    activity_score: Optional[float] = Field(default=None, description="Total activity score")
