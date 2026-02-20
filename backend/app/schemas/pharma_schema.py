from pydantic import BaseModel, Field, validator
from typing import Dict, List, Optional, Any
from datetime import datetime

class RiskAssessment(BaseModel):
    risk_label: str
    confidence_score: float
    severity: str
    risk_score: Optional[float] = None
    risk_level: Optional[str] = None

class PharmacogenomicProfile(BaseModel):
    primary_gene: str
    diplotype: str
    phenotype: str
    detected_variants: list = []
    activity_score: Optional[float] = None
    allele_scores: Optional[Dict[str, float]] = None
    notes: Optional[str] = None

class LLMExplanation(BaseModel):
    summary: str
    supporting_facts: Optional[List[str]] = None
    mechanism: Optional[str] = None
    llm_available: Optional[bool] = None
    llm_model: Optional[str] = None

class ClinicalRecommendation(BaseModel):
    text: str
    source: str = "CPIC Guidelines"
    implication: Optional[str] = None
    recommendation_url: Optional[str] = None

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
            datetime.fromisoformat(v.replace('Z', '+00:00'))
            return v
        except ValueError:
            raise ValueError("Timestamp must be a valid ISO 8601 string")
