"""
Recommendation Engine - Generate structured, context-aware clinical recommendations.

Features:
- Multi-tier recommendation logic (scaled by severity)
- Structured output with reasoning factors
- Drug-specific guidance and alternatives
- Monitoring priorities
- Deterministic recommendation mapping
"""

from typing import Dict, List, Optional
from pydantic import BaseModel, Field


# ============================================================================
# Data Models
# ============================================================================

class ClinicalGuidance(BaseModel):
    """Structured clinical guidance."""
    primary_action: str = Field(..., description="Primary recommended action")
    alternatives: List[str] = Field(default_factory=list, description="Alternative drugs")
    dosing_guidance: Optional[str] = Field(None, description="Specific dosing instructions")
    monitoring: Optional[str] = Field(None, description="Monitoring recommendations")
    urgency: str = Field(..., description="Urgency level: critical/high/moderate/low")


class StructuredRecommendation(BaseModel):
    """Complete structured recommendation with all components."""
    variant_id: str = Field(..., description="Variant identifier (gene_diplotype)")
    risk_score: float = Field(..., ge=0.0, le=100.0, description="Numeric risk score")
    risk_level: str = Field(..., description="Categorical risk level")
    confidence_score: float = Field(..., ge=0.0, le=1.0, description="Confidence score")
    confidence_interval: Optional[List[float]] = Field(None, description="[lower, upper] bounds")
    reasoning_factors: List[str] = Field(..., description="Explanation of risk factors")
    clinical_recommendation: ClinicalGuidance = Field(..., description="Clinical guidance")
    monitoring_priority: str = Field(..., description="Monitoring priority level")
    cpic_url: Optional[str] = Field(None, description="CPIC guideline URL")


# ============================================================================
# Drug-Specific Alternative Mappings
# ============================================================================

DRUG_ALTERNATIVES: Dict[str, List[str]] = {
    # Antiplatelets
    "clopidogrel": ["Prasugrel", "Ticagrelor"],

    # Opioids
    "codeine": ["Morphine", "Hydromorphone", "Oxycodone"],
    "tramadol": ["Morphine", "Hydromorphone"],

    # Antidepressants (CYP2D6)
    "amitriptyline": ["Citalopram", "Sertraline", "Venlafaxine"],
    "nortriptyline": ["Citalopram", "Sertraline"],
    "venlafaxine": ["Sertraline", "Escitalopram"],

    # Antidepressants (CYP2C19)
    "citalopram": ["Sertraline", "Venlafaxine"],
    "escitalopram": ["Sertraline", "Venlafaxine"],

    # Proton pump inhibitors
    "omeprazole": ["Pantoprazole", "Rabeprazole"],
    "lansoprazole": ["Pantoprazole", "Rabeprazole"],

    # Thiopurines
    "azathioprine": ["Mycophenolate mofetil", "Methotrexate"],
    "mercaptopurine": ["Methotrexate"],
    "thioguanine": ["Mercaptopurine (with dose adjustment)"],

    # Fluoropyrimidines
    "fluorouracil": ["Capecitabine (with dose adjustment)", "Raltitrexed"],
    "5-fluorouracil": ["Capecitabine (with dose adjustment)", "Raltitrexed"],
    "capecitabine": ["Raltitrexed", "Fluorouracil (with dose adjustment)"],

    # Anticoagulants
    "warfarin": ["Apixaban", "Rivaroxaban", "Dabigatran"],

    # Oncology
    "irinotecan": ["Oxaliplatin-based regimen"],
}


# ============================================================================
# Phenotype Explanations
# ============================================================================

def explain_phenotype_impact(phenotype: str, drug: str) -> str:
    """Generate explanation of how phenotype affects drug response."""

    explanations = {
        "PM": f"{drug} metabolism significantly reduced — risk of toxicity or lack of efficacy",
        "Poor Metabolizer": f"{drug} metabolism significantly reduced — risk of toxicity or lack of efficacy",

        "IM": f"{drug} metabolism moderately reduced — dose adjustment likely needed",
        "Intermediate Metabolizer": f"{drug} metabolism moderately reduced — dose adjustment likely needed",

        "RM": f"{drug} metabolism increased — risk of subtherapeutic levels or toxicity",
        "Rapid Metabolizer": f"{drug} metabolism increased — risk of subtherapeutic levels or toxicity",

        "UM": f"{drug} metabolism greatly increased — risk of treatment failure or toxicity",
        "Ultra-rapid Metabolizer": f"{drug} metabolism greatly increased — risk of treatment failure or toxicity",

        "NM": f"Normal {drug} metabolism expected",
        "Normal Metabolizer": f"Normal {drug} metabolism expected",
        "Extensive Metabolizer": f"Normal {drug} metabolism expected",

        "Low Activity": f"Reduced {drug} enzyme activity — high toxicity risk",
        "Intermediate Activity": f"Moderately reduced {drug} enzyme activity — dose adjustment needed",
        "Normal Activity": f"Normal {drug} enzyme activity",
        "High Activity": f"Increased {drug} enzyme activity — standard dosing appropriate",

        "Poor Function": f"Severely reduced {drug} metabolism — high toxicity risk",
        "Intermediate Function": f"Moderately reduced {drug} metabolism — dose adjustment needed",
        "Normal Function": f"Normal {drug} metabolism",

        "Indeterminate": f"{drug} metabolism phenotype could not be determined from available data",
        "Unknown": f"{drug} metabolism phenotype unknown",
    }

    return explanations.get(phenotype, f"{phenotype} phenotype")


# ============================================================================
# Recommendation Generator
# ============================================================================

class RecommendationEngine:
    """
    Generate structured, context-aware clinical recommendations.

    Maps risk scores to actionable guidance with:
    - Tiered recommendations (critical/high/moderate/low/none)
    - Drug-specific alternatives
    - Monitoring priorities
    - Detailed reasoning
    """

    def __init__(self, drug_alternatives: Optional[Dict[str, List[str]]] = None):
        self.drug_alternatives = drug_alternatives or DRUG_ALTERNATIVES

    def generate_recommendation(
        self,
        risk_score: float,
        risk_level: str,
        drug: str,
        gene: str,
        diplotype: str,
        phenotype: str,
        confidence_score: float,
        cpic_text: str = "",
        cpic_implication: str = "",
        cpic_url: Optional[str] = None,
    ) -> StructuredRecommendation:
        """
        Generate complete structured recommendation.

        Args:
            risk_score: Numeric risk score (0-100)
            risk_level: Categorical risk level
            drug: Drug name
            gene: Gene symbol
            diplotype: Diplotype (e.g., *1/*2)
            phenotype: Phenotype (e.g., IM)
            confidence_score: Confidence (0-1)
            cpic_text: CPIC recommendation text
            cpic_implication: CPIC clinical implication
            cpic_url: CPIC guideline URL

        Returns:
            StructuredRecommendation object
        """
        # Generate clinical guidance based on risk score tiers
        clinical_guidance = self._generate_clinical_guidance(
            risk_score=risk_score,
            risk_level=risk_level,
            drug=drug,
            phenotype=phenotype,
            cpic_text=cpic_text,
        )

        # Generate monitoring priority
        monitoring_priority = self._determine_monitoring_priority(risk_score)

        # Generate reasoning factors
        reasoning_factors = self._build_reasoning_factors(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            drug=drug,
            risk_score=risk_score,
            risk_level=risk_level,
            confidence_score=confidence_score,
            cpic_implication=cpic_implication,
        )

        # Estimate confidence interval (±0.05 for illustration)
        confidence_interval = [
            max(0.0, confidence_score - 0.05),
            min(1.0, confidence_score + 0.05),
        ]

        return StructuredRecommendation(
            variant_id=f"{gene}_{diplotype}",
            risk_score=round(risk_score, 2),
            risk_level=risk_level,
            confidence_score=round(confidence_score, 4),
            confidence_interval=[round(ci, 4) for ci in confidence_interval],
            reasoning_factors=reasoning_factors,
            clinical_recommendation=clinical_guidance,
            monitoring_priority=monitoring_priority,
            cpic_url=cpic_url,
        )

    def _generate_clinical_guidance(
        self,
        risk_score: float,
        risk_level: str,
        drug: str,
        phenotype: str,
        cpic_text: str,
    ) -> ClinicalGuidance:
        """Generate clinical guidance based on risk tier."""

        # Get alternatives for this drug
        alternatives = self.drug_alternatives.get(drug.lower(), [])

        # Tier-specific recommendations
        if risk_score >= 90:
            # CRITICAL (90-100)
            primary_action = f"AVOID {drug} — contraindicated or high risk of severe adverse event"
            urgency = "critical"
            monitoring = (
                f"If {drug} already prescribed, discontinue and switch to alternative immediately. "
                "Monitor for adverse events."
            )

        elif risk_score >= 70:
            # HIGH (70-89)
            primary_action = f"Consider alternative to {drug} as first-line therapy"
            urgency = "high"
            monitoring = (
                f"If {drug} is used, monitor closely for efficacy and adverse events. "
                "Consider therapeutic drug monitoring if available."
            )

        elif risk_score >= 40:
            # MODERATE (40-69)
            primary_action = f"Adjust {drug} dosing per pharmacogenomic guidelines"
            urgency = "moderate"
            monitoring = (
                f"Monitor {drug} response and adjust dose as needed. "
                "Standard monitoring protocols apply."
            )

        elif risk_score >= 20:
            # LOW (20-39)
            primary_action = f"Standard {drug} dosing acceptable with awareness of genotype"
            urgency = "low"
            monitoring = f"Standard clinical monitoring for {drug}"

        else:
            # NONE (0-19)
            primary_action = f"No pharmacogenomic concerns for {drug}"
            urgency = "low"
            monitoring = f"Standard clinical monitoring for {drug}"

        # Extract dosing guidance from CPIC text if available
        dosing_guidance = self._extract_dosing_guidance(cpic_text)

        return ClinicalGuidance(
            primary_action=primary_action,
            alternatives=alternatives,
            dosing_guidance=dosing_guidance,
            monitoring=monitoring,
            urgency=urgency,
        )

    def _determine_monitoring_priority(self, risk_score: float) -> str:
        """Map risk score to monitoring priority."""
        if risk_score >= 90:
            return "immediate_action_required"
        elif risk_score >= 70:
            return "review_before_prescribing"
        elif risk_score >= 40:
            return "standard_monitoring"
        else:
            return "routine_follow_up"

    def _build_reasoning_factors(
        self,
        gene: str,
        diplotype: str,
        phenotype: str,
        drug: str,
        risk_score: float,
        risk_level: str,
        confidence_score: float,
        cpic_implication: str,
    ) -> List[str]:
        """Build list of reasoning factors explaining the risk assessment."""

        factors = []

        # 1. Genotype/phenotype
        factors.append(f"{gene} genotype: {diplotype} ({phenotype})")

        # 2. Phenotype explanation
        explanation = explain_phenotype_impact(phenotype, drug)
        factors.append(explanation)

        # 3. Risk score
        factors.append(f"Risk score: {risk_score:.1f}/100 ({risk_level} severity)")

        # 4. Confidence
        confidence_pct = int(confidence_score * 100)
        factors.append(f"Confidence: {confidence_pct}% based on variant quality and data coverage")

        # 5. CPIC implication (if available)
        if cpic_implication and cpic_implication != "Standard considerations apply":
            factors.append(f"Clinical implication: {cpic_implication}")

        return factors

    def _extract_dosing_guidance(self, cpic_text: str) -> Optional[str]:
        """Extract dosing-specific guidance from CPIC text."""
        if not cpic_text:
            return None

        # Look for dosing keywords
        dosing_keywords = [
            "dose", "dosing", "mg", "%", "reduction", "increase",
            "initiate", "starting dose", "maintenance",
        ]

        text_lower = cpic_text.lower()
        if any(kw in text_lower for kw in dosing_keywords):
            # Return CPIC text as dosing guidance
            return cpic_text

        return None


# ============================================================================
# Convenience Functions
# ============================================================================

def generate_recommendation(
    risk_score: float,
    risk_level: str,
    drug: str,
    gene: str,
    diplotype: str,
    phenotype: str,
    confidence_score: float,
    cpic_text: str = "",
    cpic_implication: str = "",
    cpic_url: Optional[str] = None,
) -> StructuredRecommendation:
    """
    Convenience function to generate recommendation.

    See RecommendationEngine.generate_recommendation for details.
    """
    engine = RecommendationEngine()
    return engine.generate_recommendation(
        risk_score=risk_score,
        risk_level=risk_level,
        drug=drug,
        gene=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        confidence_score=confidence_score,
        cpic_text=cpic_text,
        cpic_implication=cpic_implication,
        cpic_url=cpic_url,
    )
