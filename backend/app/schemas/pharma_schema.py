from pydantic import BaseModel, Field, validator
from typing import Dict, List, Optional, Any
from datetime import datetime

class RiskAssessment(BaseModel):
    risk_label: str
    confidence_score: float
    severity: str

class PharmacogenomicProfile(BaseModel):
    primary_gene: str
    diplotype: str
    phenotype: str
    detected_variants: list = []

class LLMExplanation(BaseModel):
    summary: str
    

class ClinicalRecommendation(BaseModel):
    action: str
    source: str = "CPIC Guidelines"

class QualityMetrics(BaseModel):
    vcf_parsing_success: bool = True
    coverage_check: str = "PASS"
    extra_metadata: Optional[Dict[str, Any]] = None

class PharmaGuardResponse(BaseModel):
    patient_id: str
    drug: str
    timestamp: str 
    risk_assessment: RiskAssessment
    pharmacogenomic_profile: PharmacogenomicProfile
    clinical_recommendation: ClinicalRecommendation
    llm_generated_explanation: LLMExplanation
    quality_metrics: QualityMetrics

    @validator('timestamp')
    def validate_timestamp(cls, v):
        try:
            # simple check if it's ISO formatted, though string type is specified
            datetime.fromisoformat(v.replace('Z', '+00:00'))
            return v
        except ValueError:
            raise ValueError("Timestamp must be a valid ISO 8601 string")
