"""
Polypharmacy Analysis API - Multi-drug pharmacogenomic risk assessment endpoint.

Endpoints:
- POST /api/polypharmacy/analyze - Analyze multiple drugs for a patient
- GET /api/polypharmacy/interactions - Get interaction database
- POST /api/polypharmacy/check-pair - Check specific drug pair interaction
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import List, Dict, Optional
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent.parent))

from app.services.pharmacogenomics.risk_engine import RiskEngine
from app.services.pharmacogenomics.multi_drug_risk import (
    MultiDrugRiskAnalyzer,
    InteractionMatrix,
    get_interaction_database,
    DrugDrugInteraction,
    InteractionType,
    InteractionSeverity,
)
from app.services.pharmacogenomics.models import PatientProfile, DiplotypeResult


router = APIRouter()


# ============================================================================
# Request/Response Models
# ============================================================================

class PolypharmacyAnalysisRequest(BaseModel):
    """Request for polypharmacy analysis."""
    patient_id: str = Field(..., description="Patient identifier")
    drugs: List[str] = Field(..., min_items=2, description="List of drugs to analyze (minimum 2)")

    # Patient pharmacogenomic profile
    diplotypes: Dict[str, Dict] = Field(
        ...,
        description="Patient diplotypes by gene (gene -> {diplotype, phenotype, confidence})"
    )

    # Optional parameters
    include_alternatives: bool = Field(True, description="Include alternative regimen suggestions")
    interaction_severity_filter: Optional[str] = Field(
        None,
        description="Filter interactions by severity (critical/major/moderate/minor)"
    )


class DrugInteractionInfo(BaseModel):
    """Drug-drug interaction information."""
    drug_a: str
    drug_b: str
    gene: str
    interaction_type: str
    severity: str
    risk_multiplier: float
    mechanism: str
    clinical_implication: str
    monitoring_recommendation: str
    affected_phenotypes: List[str]


class DrugRiskSummary(BaseModel):
    """Summary of individual drug risk."""
    drug: str
    gene: str
    diplotype: str
    phenotype: str
    risk_score: float
    risk_level: str
    confidence_score: float
    recommendation: str


class PolypharmacyAnalysisResponse(BaseModel):
    """Response for polypharmacy analysis."""
    patient_id: str
    drugs_analyzed: List[str]
    analysis_timestamp: str

    # Combined risk
    combined_risk_score: float = Field(..., ge=0.0, le=100.0)
    combined_risk_level: str
    combined_confidence: float = Field(..., ge=0.0, le=1.0)

    # Individual drug risks
    individual_drug_risks: List[DrugRiskSummary]

    # Interactions
    detected_interactions: List[DrugInteractionInfo]
    interaction_count: int
    highest_interaction_severity: str

    # Priority and warnings
    highest_priority_drug: str
    critical_warnings: List[str]

    # Recommendations
    polypharmacy_recommendation: str
    monitoring_priority: str
    alternative_regimens: Optional[List[Dict]] = None

    # Risk contributions
    risk_contributions: List[Dict]


class DrugPairCheckRequest(BaseModel):
    """Request to check specific drug pair interaction."""
    drug_a: str
    drug_b: str
    gene: Optional[str] = None
    phenotype: Optional[str] = None


class DrugPairCheckResponse(BaseModel):
    """Response for drug pair interaction check."""
    drug_a: str
    drug_b: str
    has_interaction: bool
    interactions: List[DrugInteractionInfo]
    recommendation: str


# ============================================================================
# Endpoints
# ============================================================================

@router.post("/analyze", response_model=PolypharmacyAnalysisResponse)
async def analyze_polypharmacy(request: PolypharmacyAnalysisRequest):
    """
    Analyze combined pharmacogenomic risk from multiple drugs.

    Detects drug-drug-gene interactions, calculates combined risk scores,
    and provides polypharmacy-specific recommendations.
    """
    from datetime import datetime

    try:
        # Initialize risk engine
        risk_engine = RiskEngine()

        # Build patient profile from diplotypes
        patient_profile = PatientProfile(
            sample_id=request.patient_id,
            diplotypes={}
        )

        for gene, data in request.diplotypes.items():
            patient_profile.diplotypes[gene] = DiplotypeResult(
                gene=gene,
                diplotype=data.get("diplotype", "Unknown"),
                phenotype=data.get("phenotype", "Unknown"),
                confidence=data.get("confidence", 0.0),
                is_indeterminate=data.get("is_indeterminate", False),
            )

        # Evaluate each drug individually
        drug_assessments = risk_engine.evaluate_multiple_drugs(
            drugs=request.drugs,
            patient_profile=patient_profile
        )

        if not drug_assessments:
            raise HTTPException(
                status_code=400,
                detail="No valid drug assessments could be generated. Check drug names and patient profile."
            )

        # Analyze multi-drug risk
        analyzer = MultiDrugRiskAnalyzer()
        multi_drug_assessment = analyzer.analyze_multi_drug_risk(drug_assessments)

        # Build individual drug risk summaries
        individual_risks = []
        for assess in drug_assessments:
            individual_risks.append(DrugRiskSummary(
                drug=assess.drug,
                gene=assess.gene,
                diplotype=assess.diplotype,
                phenotype=assess.phenotype,
                risk_score=assess.risk.risk_score or 0.0,
                risk_level=assess.risk.risk_level or assess.risk.severity,
                confidence_score=assess.risk.confidence_score,
                recommendation=assess.recommendation.text,
            ))

        # Build interaction info
        interactions_info = []
        for interaction in multi_drug_assessment.detected_interactions:
            # Filter by severity if requested
            if request.interaction_severity_filter:
                if interaction.severity.value != request.interaction_severity_filter:
                    continue

            interactions_info.append(DrugInteractionInfo(
                drug_a=interaction.drug_a,
                drug_b=interaction.drug_b,
                gene=interaction.gene,
                interaction_type=interaction.interaction_type.value,
                severity=interaction.severity.value,
                risk_multiplier=interaction.risk_multiplier,
                mechanism=interaction.mechanism,
                clinical_implication=interaction.clinical_implication,
                monitoring_recommendation=interaction.monitoring_recommendation,
                affected_phenotypes=interaction.affected_phenotypes,
            ))

        # Build risk contributions
        risk_contributions = []
        for contrib in multi_drug_assessment.individual_risks:
            # Calculate contribution percentage
            total_score = sum(
                a.risk.risk_score or 0.0
                for a in drug_assessments
            )
            contrib_pct = (
                (contrib.risk.risk_score or 0.0) / total_score * 100
                if total_score > 0 else 0.0
            )

            risk_contributions.append({
                "drug": contrib.drug,
                "individual_risk_score": contrib.risk.risk_score or 0.0,
                "contribution_percentage": round(contrib_pct, 1),
                "risk_level": contrib.risk.risk_level or contrib.risk.severity,
            })

        # Sort by contribution
        risk_contributions.sort(
            key=lambda x: x["contribution_percentage"],
            reverse=True
        )

        return PolypharmacyAnalysisResponse(
            patient_id=request.patient_id,
            drugs_analyzed=multi_drug_assessment.drugs,
            analysis_timestamp=datetime.now().isoformat(),
            combined_risk_score=multi_drug_assessment.combined_risk_score,
            combined_risk_level=multi_drug_assessment.combined_risk_level,
            combined_confidence=multi_drug_assessment.combined_confidence,
            individual_drug_risks=individual_risks,
            detected_interactions=interactions_info,
            interaction_count=len(interactions_info),
            highest_interaction_severity=multi_drug_assessment.highest_interaction_severity,
            highest_priority_drug=multi_drug_assessment.highest_priority_drug,
            critical_warnings=multi_drug_assessment.critical_warnings,
            polypharmacy_recommendation=multi_drug_assessment.polypharmacy_recommendation,
            monitoring_priority=multi_drug_assessment.monitoring_priority,
            alternative_regimens=(
                multi_drug_assessment.alternative_regimens
                if request.include_alternatives else None
            ),
            risk_contributions=risk_contributions,
        )

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error analyzing polypharmacy: {str(e)}"
        )


@router.get("/interactions", response_model=List[DrugInteractionInfo])
async def get_interactions(
    gene: Optional[str] = None,
    severity: Optional[str] = None,
    drug: Optional[str] = None,
):
    """
    Get drug-drug-gene interaction database.

    Optional filters:
    - gene: Filter by specific gene (e.g., CYP2D6)
    - severity: Filter by severity (critical/major/moderate/minor)
    - drug: Filter by specific drug
    """
    interaction_db = get_interaction_database()

    # Apply filters
    filtered = interaction_db

    if gene:
        filtered = [i for i in filtered if i.gene == gene]

    if severity:
        try:
            sev = InteractionSeverity(severity.lower())
            filtered = [i for i in filtered if i.severity == sev]
        except ValueError:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid severity: {severity}. Use critical/major/moderate/minor"
            )

    if drug:
        drug_lower = drug.lower()
        filtered = [
            i for i in filtered
            if drug_lower in (i.drug_a.lower(), i.drug_b.lower())
        ]

    # Convert to response model
    interactions_info = []
    for interaction in filtered:
        interactions_info.append(DrugInteractionInfo(
            drug_a=interaction.drug_a,
            drug_b=interaction.drug_b,
            gene=interaction.gene,
            interaction_type=interaction.interaction_type.value,
            severity=interaction.severity.value,
            risk_multiplier=interaction.risk_multiplier,
            mechanism=interaction.mechanism,
            clinical_implication=interaction.clinical_implication,
            monitoring_recommendation=interaction.monitoring_recommendation,
            affected_phenotypes=interaction.affected_phenotypes,
        ))

    return interactions_info


@router.post("/check-pair", response_model=DrugPairCheckResponse)
async def check_drug_pair(request: DrugPairCheckRequest):
    """
    Check for interaction between a specific drug pair.

    Optionally filter by gene and phenotype.
    """
    matrix = InteractionMatrix(get_interaction_database())

    # Get interactions for this pair
    interactions = matrix.get_interactions_for_pair(
        request.drug_a,
        request.drug_b
    )

    # Filter by gene if specified
    if request.gene:
        interactions = [i for i in interactions if i.gene == request.gene]

    # Filter by phenotype if specified
    if request.phenotype:
        interactions = [
            i for i in interactions
            if (not i.affected_phenotypes or
                request.phenotype in i.affected_phenotypes)
        ]

    has_interaction = len(interactions) > 0

    # Build response
    interactions_info = []
    for interaction in interactions:
        interactions_info.append(DrugInteractionInfo(
            drug_a=interaction.drug_a,
            drug_b=interaction.drug_b,
            gene=interaction.gene,
            interaction_type=interaction.interaction_type.value,
            severity=interaction.severity.value,
            risk_multiplier=interaction.risk_multiplier,
            mechanism=interaction.mechanism,
            clinical_implication=interaction.clinical_implication,
            monitoring_recommendation=interaction.monitoring_recommendation,
            affected_phenotypes=interaction.affected_phenotypes,
        ))

    # Generate recommendation
    if has_interaction:
        highest_severity = max(
            interactions,
            key=lambda i: {
                InteractionSeverity.CRITICAL: 4,
                InteractionSeverity.MAJOR: 3,
                InteractionSeverity.MODERATE: 2,
                InteractionSeverity.MINOR: 1,
            }.get(i.severity, 0)
        )

        if highest_severity.severity == InteractionSeverity.CRITICAL:
            recommendation = (
                f"AVOID combination: {highest_severity.mechanism}. "
                f"{highest_severity.monitoring_recommendation}"
            )
        elif highest_severity.severity == InteractionSeverity.MAJOR:
            recommendation = (
                f"Use with caution: {highest_severity.mechanism}. "
                f"{highest_severity.monitoring_recommendation}"
            )
        else:
            recommendation = (
                f"Monitor: {highest_severity.mechanism}. "
                f"{highest_severity.monitoring_recommendation}"
            )
    else:
        recommendation = "No known pharmacogenomic interaction detected for this drug pair."

    return DrugPairCheckResponse(
        drug_a=request.drug_a,
        drug_b=request.drug_b,
        has_interaction=has_interaction,
        interactions=interactions_info,
        recommendation=recommendation,
    )


@router.get("/interaction-summary")
async def get_interaction_summary():
    """Get summary statistics of interaction database."""
    interaction_db = get_interaction_database()

    # Count by severity
    severity_counts = {
        "critical": 0,
        "major": 0,
        "moderate": 0,
        "minor": 0,
    }

    for interaction in interaction_db:
        severity_counts[interaction.severity.value] += 1

    # Count by gene
    gene_counts = {}
    for interaction in interaction_db:
        gene_counts[interaction.gene] = gene_counts.get(interaction.gene, 0) + 1

    # Count by interaction type
    type_counts = {}
    for interaction in interaction_db:
        itype = interaction.interaction_type.value
        type_counts[itype] = type_counts.get(itype, 0) + 1

    # Get unique drugs
    drugs = set()
    for interaction in interaction_db:
        drugs.add(interaction.drug_a)
        drugs.add(interaction.drug_b)

    return {
        "total_interactions": len(interaction_db),
        "unique_drugs": len(drugs),
        "genes_covered": list(gene_counts.keys()),
        "severity_distribution": severity_counts,
        "gene_distribution": gene_counts,
        "interaction_type_distribution": type_counts,
    }
