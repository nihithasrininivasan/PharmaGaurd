"""
Multi-Drug Risk Integration - Polypharmacy risk assessment and drug-drug interactions.

Features:
- Combined risk aggregation from multiple drugs
- Drug-drug interaction detection (pharmacogenomic)
- Synergistic/compounding risk identification
- Polypharmacy-specific recommendations
- Interaction severity scoring
- Unified risk prioritization
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, field
from enum import Enum
import math

from .models import DrugAssessment, RiskAssessment, ClinicalRecommendation


# ============================================================================
# Interaction Types
# ============================================================================

class InteractionType(str, Enum):
    """Types of drug-drug-gene interactions."""
    SYNERGISTIC_TOXICITY = "synergistic_toxicity"          # Combined toxicity > sum
    SYNERGISTIC_INEFFICACY = "synergistic_inefficacy"      # Combined inefficacy > sum
    COMPETITIVE_METABOLISM = "competitive_metabolism"       # Compete for same enzyme
    ENZYME_INDUCTION = "enzyme_induction"                  # One induces enzyme for other
    ENZYME_INHIBITION = "enzyme_inhibition"                # One inhibits enzyme for other
    ADDITIVE_RISK = "additive_risk"                        # Risks simply add
    PATHWAY_OVERLAP = "pathway_overlap"                    # Share metabolic pathway
    COUNTERACTING = "counteracting"                        # One reduces risk of other


class InteractionSeverity(str, Enum):
    """Severity of drug-drug interaction."""
    CRITICAL = "critical"      # Life-threatening, contraindicated
    MAJOR = "major"            # Significant risk, avoid if possible
    MODERATE = "moderate"      # Monitor closely, may need adjustment
    MINOR = "minor"            # Informational, standard monitoring


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class DrugDrugInteraction:
    """Represents a pharmacogenomic drug-drug interaction."""
    drug_a: str
    drug_b: str
    gene: str
    interaction_type: InteractionType
    severity: InteractionSeverity

    # Risk modification
    risk_multiplier: float = 1.0        # Multiplicative factor for combined risk
    confidence_penalty: float = 0.0     # Reduce confidence due to interaction

    # Clinical details
    mechanism: str = ""
    clinical_implication: str = ""
    monitoring_recommendation: str = ""

    # Phenotype-specific
    affected_phenotypes: List[str] = field(default_factory=list)


@dataclass
class MultiDrugRiskAssessment:
    """Combined risk assessment for multiple drugs."""
    drugs: List[str]
    individual_risks: List[DrugAssessment]

    # Aggregated scores
    combined_risk_score: float              # 0-100 continuous score
    combined_risk_level: str                # critical/major/moderate/minor/none
    combined_confidence: float              # 0-1

    # Interactions
    detected_interactions: List[DrugDrugInteraction]
    interaction_count: int
    highest_interaction_severity: str

    # Prioritization
    highest_priority_drug: str
    critical_warnings: List[str]

    # Recommendations
    polypharmacy_recommendation: str
    monitoring_priority: str
    alternative_regimens: List[Dict]


@dataclass
class RiskContribution:
    """Individual drug's contribution to combined risk."""
    drug: str
    individual_risk_score: float
    contribution_to_combined: float     # Percentage contribution
    interaction_multiplier: float       # Multiplier from interactions
    priority_rank: int                  # 1 = highest priority


# ============================================================================
# Drug-Drug-Gene Interaction Database
# ============================================================================

# Known pharmacogenomic drug-drug interactions
KNOWN_INTERACTIONS: List[DrugDrugInteraction] = [

    # CYP2D6 interactions
    DrugDrugInteraction(
        drug_a="codeine",
        drug_b="fluoxetine",  # CYP2D6 inhibitor
        gene="CYP2D6",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=1.5,
        confidence_penalty=0.1,
        mechanism="Fluoxetine inhibits CYP2D6, reducing codeine conversion to morphine",
        clinical_implication="Reduced analgesic efficacy in CYP2D6 normal/extensive metabolizers",
        monitoring_recommendation="Monitor pain control; may need alternative analgesic",
        affected_phenotypes=["NM", "EM", "UM"],
    ),

    DrugDrugInteraction(
        drug_a="tramadol",
        drug_b="paroxetine",  # CYP2D6 inhibitor
        gene="CYP2D6",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=1.4,
        confidence_penalty=0.1,
        mechanism="Paroxetine inhibits CYP2D6, reducing tramadol activation",
        clinical_implication="Reduced analgesic efficacy",
        monitoring_recommendation="Consider alternative pain management",
        affected_phenotypes=["NM", "EM"],
    ),

    DrugDrugInteraction(
        drug_a="amitriptyline",
        drug_b="bupropion",  # CYP2D6 inhibitor
        gene="CYP2D6",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.MODERATE,
        risk_multiplier=1.3,
        confidence_penalty=0.05,
        mechanism="Bupropion inhibits CYP2D6, increasing amitriptyline levels",
        clinical_implication="Increased risk of tricyclic toxicity",
        monitoring_recommendation="Monitor for anticholinergic effects, consider dose reduction",
        affected_phenotypes=["NM", "EM"],
    ),

    # CYP2C19 interactions
    DrugDrugInteraction(
        drug_a="clopidogrel",
        drug_b="omeprazole",  # CYP2C19 inhibitor
        gene="CYP2C19",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=1.8,
        confidence_penalty=0.15,
        mechanism="Omeprazole inhibits CYP2C19, reducing clopidogrel activation",
        clinical_implication="Reduced antiplatelet efficacy, increased cardiovascular risk",
        monitoring_recommendation="Use alternative PPI (pantoprazole) or H2 blocker",
        affected_phenotypes=["NM", "RM", "UM"],
    ),

    DrugDrugInteraction(
        drug_a="clopidogrel",
        drug_b="esomeprazole",  # CYP2C19 inhibitor
        gene="CYP2C19",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=1.7,
        confidence_penalty=0.15,
        mechanism="Esomeprazole inhibits CYP2C19, reducing clopidogrel activation",
        clinical_implication="Reduced antiplatelet efficacy",
        monitoring_recommendation="Avoid concurrent use; use pantoprazole if PPI needed",
        affected_phenotypes=["NM", "RM", "UM"],
    ),

    DrugDrugInteraction(
        drug_a="citalopram",
        drug_b="omeprazole",
        gene="CYP2C19",
        interaction_type=InteractionType.COMPETITIVE_METABOLISM,
        severity=InteractionSeverity.MODERATE,
        risk_multiplier=1.3,
        confidence_penalty=0.05,
        mechanism="Both metabolized by CYP2C19; competitive inhibition possible",
        clinical_implication="May increase citalopram levels in poor metabolizers",
        monitoring_recommendation="Monitor for QT prolongation and serotonergic effects",
        affected_phenotypes=["PM", "IM"],
    ),

    # CYP2C9 interactions
    DrugDrugInteraction(
        drug_a="warfarin",
        drug_b="fluconazole",  # CYP2C9 inhibitor
        gene="CYP2C9",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.CRITICAL,
        risk_multiplier=2.5,
        confidence_penalty=0.2,
        mechanism="Fluconazole strongly inhibits CYP2C9, increasing warfarin levels",
        clinical_implication="Significantly increased bleeding risk",
        monitoring_recommendation="Reduce warfarin dose 25-50%, monitor INR closely (every 2-3 days)",
        affected_phenotypes=["NM", "IM", "PM"],
    ),

    DrugDrugInteraction(
        drug_a="warfarin",
        drug_b="amiodarone",  # CYP2C9 inhibitor
        gene="CYP2C9",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.CRITICAL,
        risk_multiplier=2.2,
        confidence_penalty=0.2,
        mechanism="Amiodarone inhibits CYP2C9, increasing warfarin exposure",
        clinical_implication="Severe bleeding risk",
        monitoring_recommendation="Reduce warfarin dose 30-50%, increase INR monitoring",
        affected_phenotypes=["NM", "IM", "PM"],
    ),

    # TPMT interactions (synergistic toxicity)
    DrugDrugInteraction(
        drug_a="azathioprine",
        drug_b="allopurinol",  # Xanthine oxidase inhibitor
        gene="TPMT",
        interaction_type=InteractionType.SYNERGISTIC_TOXICITY,
        severity=InteractionSeverity.CRITICAL,
        risk_multiplier=3.0,
        confidence_penalty=0.25,
        mechanism="Allopurinol blocks alternative pathway, shunting to TPMT pathway",
        clinical_implication="Severe myelosuppression risk, especially in TPMT poor metabolizers",
        monitoring_recommendation="Reduce azathioprine dose by 75%, monitor CBC weekly initially",
        affected_phenotypes=["PM", "IM"],
    ),

    DrugDrugInteraction(
        drug_a="mercaptopurine",
        drug_b="allopurinol",
        gene="TPMT",
        interaction_type=InteractionType.SYNERGISTIC_TOXICITY,
        severity=InteractionSeverity.CRITICAL,
        risk_multiplier=3.0,
        confidence_penalty=0.25,
        mechanism="Allopurinol increases 6-MP levels 3-4 fold",
        clinical_implication="Life-threatening myelosuppression",
        monitoring_recommendation="Reduce mercaptopurine dose by 65-75%, monitor CBC closely",
        affected_phenotypes=["PM", "IM"],
    ),

    # DPYD interactions
    DrugDrugInteraction(
        drug_a="fluorouracil",
        drug_b="sorivudine",  # DPD inhibitor
        gene="DPYD",
        interaction_type=InteractionType.ENZYME_INHIBITION,
        severity=InteractionSeverity.CRITICAL,
        risk_multiplier=4.0,
        confidence_penalty=0.3,
        mechanism="Sorivudine irreversibly inhibits DPD, causing 5-FU accumulation",
        clinical_implication="Fatal toxicity reported; absolutely contraindicated",
        monitoring_recommendation="AVOID combination; wait 4 weeks after sorivudine before 5-FU",
        affected_phenotypes=["NM", "IM", "PM"],
    ),

    # Multi-gene interactions (warfarin + CYP2C9 + VKORC1)
    DrugDrugInteraction(
        drug_a="warfarin",
        drug_b="aspirin",
        gene="CYP2C9",  # Primary gene
        interaction_type=InteractionType.ADDITIVE_RISK,
        severity=InteractionSeverity.MAJOR,
        risk_multiplier=1.6,
        confidence_penalty=0.1,
        mechanism="Additive anticoagulant effects, especially in CYP2C9 poor metabolizers",
        clinical_implication="Increased bleeding risk",
        monitoring_recommendation="Monitor INR more frequently, consider lower warfarin dose",
        affected_phenotypes=["PM", "IM"],
    ),
]


def get_interaction_database() -> List[DrugDrugInteraction]:
    """Get the interaction database."""
    return KNOWN_INTERACTIONS


# ============================================================================
# Multi-Drug Risk Analyzer
# ============================================================================

class MultiDrugRiskAnalyzer:
    """
    Analyze combined pharmacogenomic risk from multiple drugs.

    Features:
    - Detect drug-drug-gene interactions
    - Calculate combined risk scores
    - Identify synergistic/compounding effects
    - Generate polypharmacy recommendations
    - Prioritize risks for clinical action
    """

    def __init__(self, interaction_db: Optional[List[DrugDrugInteraction]] = None):
        self.interaction_db = interaction_db or KNOWN_INTERACTIONS

    def analyze_multi_drug_risk(
        self,
        drug_assessments: List[DrugAssessment]
    ) -> MultiDrugRiskAssessment:
        """
        Analyze combined risk from multiple drugs.

        Args:
            drug_assessments: List of individual drug assessments

        Returns:
            MultiDrugRiskAssessment with combined analysis
        """
        if not drug_assessments:
            raise ValueError("At least one drug assessment required")

        # Extract drug names
        drugs = [a.drug for a in drug_assessments]

        # Detect interactions
        interactions = self._detect_interactions(drug_assessments)

        # Calculate combined risk score
        combined_score, risk_contributions = self._calculate_combined_risk(
            drug_assessments,
            interactions
        )

        # Determine combined risk level
        combined_level = self._get_risk_level(combined_score)

        # Calculate combined confidence
        combined_confidence = self._calculate_combined_confidence(
            drug_assessments,
            interactions
        )

        # Identify highest priority drug
        highest_priority = self._identify_highest_priority(risk_contributions)

        # Generate critical warnings
        critical_warnings = self._generate_critical_warnings(
            drug_assessments,
            interactions
        )

        # Determine interaction severity
        interaction_severity = self._get_highest_interaction_severity(interactions)

        # Generate polypharmacy recommendation
        poly_recommendation = self._generate_polypharmacy_recommendation(
            drug_assessments,
            interactions,
            combined_score
        )

        # Determine monitoring priority
        monitoring_priority = self._determine_monitoring_priority(
            combined_score,
            interactions
        )

        # Generate alternative regimens
        alternatives = self._generate_alternative_regimens(
            drug_assessments,
            interactions
        )

        return MultiDrugRiskAssessment(
            drugs=drugs,
            individual_risks=drug_assessments,
            combined_risk_score=combined_score,
            combined_risk_level=combined_level,
            combined_confidence=combined_confidence,
            detected_interactions=interactions,
            interaction_count=len(interactions),
            highest_interaction_severity=interaction_severity,
            highest_priority_drug=highest_priority,
            critical_warnings=critical_warnings,
            polypharmacy_recommendation=poly_recommendation,
            monitoring_priority=monitoring_priority,
            alternative_regimens=alternatives,
        )

    def _detect_interactions(
        self,
        drug_assessments: List[DrugAssessment]
    ) -> List[DrugDrugInteraction]:
        """Detect drug-drug-gene interactions."""
        interactions = []

        # Check all drug pairs
        for i, assess_a in enumerate(drug_assessments):
            for assess_b in drug_assessments[i+1:]:
                # Check both orderings in interaction database
                interaction = self._find_interaction(
                    assess_a.drug,
                    assess_b.drug,
                    assess_a.gene
                )

                if interaction:
                    # Check if phenotype is affected
                    phenotype_a = assess_a.phenotype
                    if (not interaction.affected_phenotypes or
                        phenotype_a in interaction.affected_phenotypes):
                        interactions.append(interaction)

        return interactions

    def _find_interaction(
        self,
        drug_a: str,
        drug_b: str,
        gene: str
    ) -> Optional[DrugDrugInteraction]:
        """Find interaction in database."""
        drug_a_lower = drug_a.lower()
        drug_b_lower = drug_b.lower()

        for interaction in self.interaction_db:
            # Check both orderings
            if (
                (interaction.drug_a.lower() == drug_a_lower and
                 interaction.drug_b.lower() == drug_b_lower and
                 interaction.gene == gene) or
                (interaction.drug_a.lower() == drug_b_lower and
                 interaction.drug_b.lower() == drug_a_lower and
                 interaction.gene == gene)
            ):
                return interaction

        return None

    def _calculate_combined_risk(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction]
    ) -> Tuple[float, List[RiskContribution]]:
        """
        Calculate combined risk score with interaction effects.

        Uses weighted geometric mean with interaction multipliers.
        """
        if not drug_assessments:
            return 0.0, []

        # Get individual risk scores
        individual_scores = []
        for assess in drug_assessments:
            score = assess.risk.risk_score if assess.risk.risk_score else 0.0
            individual_scores.append(score)

        # Apply interaction multipliers
        interaction_multipliers = self._calculate_interaction_multipliers(
            drug_assessments,
            interactions
        )

        # Calculate combined score using weighted approach
        # For synergistic interactions: use max with multiplier
        # For additive: use weighted sum
        # For competitive: use modified average

        if interactions:
            # Has interactions - use interaction-aware scoring
            combined_score = self._calculate_interaction_aware_score(
                individual_scores,
                interaction_multipliers,
                interactions
            )
        else:
            # No interactions - use simple weighted maximum
            # (highest risk dominates, but other risks contribute)
            max_score = max(individual_scores)
            avg_score = sum(individual_scores) / len(individual_scores)
            combined_score = 0.7 * max_score + 0.3 * avg_score

        # Build risk contributions
        contributions = []
        total_contribution = sum(individual_scores)

        for i, assess in enumerate(drug_assessments):
            contribution_pct = (
                (individual_scores[i] / total_contribution * 100)
                if total_contribution > 0 else 0.0
            )

            contributions.append(RiskContribution(
                drug=assess.drug,
                individual_risk_score=individual_scores[i],
                contribution_to_combined=contribution_pct,
                interaction_multiplier=interaction_multipliers.get(assess.drug, 1.0),
                priority_rank=0,  # Will be set later
            ))

        # Rank by contribution
        contributions.sort(key=lambda c: c.contribution_to_combined, reverse=True)
        for i, contrib in enumerate(contributions):
            contrib.priority_rank = i + 1

        return combined_score, contributions

    def _calculate_interaction_multipliers(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction]
    ) -> Dict[str, float]:
        """Calculate risk multipliers for each drug due to interactions."""
        multipliers = {assess.drug: 1.0 for assess in drug_assessments}

        for interaction in interactions:
            # Apply multiplier to both drugs
            drug_a = interaction.drug_a
            drug_b = interaction.drug_b

            if drug_a in multipliers:
                multipliers[drug_a] = max(
                    multipliers[drug_a],
                    interaction.risk_multiplier
                )

            if drug_b in multipliers:
                multipliers[drug_b] = max(
                    multipliers[drug_b],
                    interaction.risk_multiplier
                )

        return multipliers

    def _calculate_interaction_aware_score(
        self,
        individual_scores: List[float],
        interaction_multipliers: Dict[str, float],
        interactions: List[DrugDrugInteraction]
    ) -> float:
        """Calculate combined score considering interaction types."""

        # Categorize interactions
        synergistic = [i for i in interactions if i.interaction_type in (
            InteractionType.SYNERGISTIC_TOXICITY,
            InteractionType.SYNERGISTIC_INEFFICACY
        )]

        inhibitory = [i for i in interactions if i.interaction_type in (
            InteractionType.ENZYME_INHIBITION,
            InteractionType.COMPETITIVE_METABOLISM
        )]

        # Base score: weighted maximum
        max_score = max(individual_scores)
        avg_score = sum(individual_scores) / len(individual_scores)
        base_score = 0.7 * max_score + 0.3 * avg_score

        # Apply synergistic multiplier if present
        if synergistic:
            max_synergy_multiplier = max(i.risk_multiplier for i in synergistic)
            base_score = base_score * max_synergy_multiplier

        # Apply inhibitory multiplier if present
        elif inhibitory:
            max_inhibit_multiplier = max(i.risk_multiplier for i in inhibitory)
            base_score = base_score * max_inhibit_multiplier

        # Clamp to [0, 100]
        return max(0.0, min(100.0, base_score))

    def _calculate_combined_confidence(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction]
    ) -> float:
        """Calculate combined confidence with interaction penalties."""

        # Get individual confidences
        confidences = [assess.risk.confidence_score for assess in drug_assessments]

        # Base combined confidence: minimum (conservative)
        base_confidence = min(confidences)

        # Apply interaction penalties
        total_penalty = sum(i.confidence_penalty for i in interactions)

        combined = base_confidence - total_penalty

        return max(0.0, min(1.0, combined))

    def _identify_highest_priority(
        self,
        contributions: List[RiskContribution]
    ) -> str:
        """Identify highest priority drug for intervention."""
        if not contributions:
            return ""

        # Highest contribution = highest priority
        return contributions[0].drug

    def _generate_critical_warnings(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction]
    ) -> List[str]:
        """Generate critical warnings for high-risk situations."""
        warnings = []

        # Critical individual risks
        for assess in drug_assessments:
            if assess.risk.risk_score and assess.risk.risk_score >= 90:
                warnings.append(
                    f"CRITICAL: {assess.drug} contraindicated for {assess.gene} "
                    f"{assess.phenotype} (risk score: {assess.risk.risk_score:.0f}/100)"
                )

        # Critical interactions
        for interaction in interactions:
            if interaction.severity == InteractionSeverity.CRITICAL:
                warnings.append(
                    f"CRITICAL INTERACTION: {interaction.drug_a} + {interaction.drug_b} "
                    f"in {interaction.gene} — {interaction.mechanism}"
                )

        return warnings

    def _get_highest_interaction_severity(
        self,
        interactions: List[DrugDrugInteraction]
    ) -> str:
        """Get highest severity among all interactions."""
        if not interactions:
            return "none"

        severity_order = {
            InteractionSeverity.CRITICAL: 4,
            InteractionSeverity.MAJOR: 3,
            InteractionSeverity.MODERATE: 2,
            InteractionSeverity.MINOR: 1,
        }

        max_severity = max(
            interactions,
            key=lambda i: severity_order.get(i.severity, 0)
        )

        return max_severity.severity.value

    def _generate_polypharmacy_recommendation(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction],
        combined_score: float
    ) -> str:
        """Generate comprehensive polypharmacy recommendation."""

        if combined_score >= 90:
            return (
                "CRITICAL POLYPHARMACY RISK: Current drug combination is contraindicated. "
                "Immediate regimen modification required. Consult pharmacist/clinical pharmacologist."
            )

        elif combined_score >= 70:
            if interactions:
                return (
                    f"HIGH POLYPHARMACY RISK: {len(interactions)} significant drug-drug-gene "
                    f"interaction(s) detected. Consider alternative regimen or intensive monitoring. "
                    f"Pharmacogenomic consultation recommended."
                )
            else:
                return (
                    "HIGH POLYPHARMACY RISK: Multiple high-risk drugs prescribed. "
                    "Review regimen and consider alternatives where possible."
                )

        elif combined_score >= 40:
            return (
                f"MODERATE POLYPHARMACY RISK: Monitor patient closely for efficacy and adverse effects. "
                f"{len(interactions)} interaction(s) detected. Dose adjustments may be needed."
            )

        else:
            if interactions:
                return (
                    f"LOW POLYPHARMACY RISK: {len(interactions)} minor interaction(s) detected. "
                    "Standard monitoring appropriate."
                )
            else:
                return "No significant polypharmacy risks detected. Standard monitoring recommended."

    def _determine_monitoring_priority(
        self,
        combined_score: float,
        interactions: List[DrugDrugInteraction]
    ) -> str:
        """Determine monitoring priority level."""

        # Check for critical interactions
        critical_interactions = [
            i for i in interactions
            if i.severity == InteractionSeverity.CRITICAL
        ]

        if critical_interactions or combined_score >= 90:
            return "immediate_action_required"
        elif combined_score >= 70:
            return "review_before_prescribing"
        elif combined_score >= 40:
            return "enhanced_monitoring_required"
        else:
            return "standard_monitoring"

    def _generate_alternative_regimens(
        self,
        drug_assessments: List[DrugAssessment],
        interactions: List[DrugDrugInteraction]
    ) -> List[Dict]:
        """Suggest alternative drug regimens to reduce risk."""
        alternatives = []

        # For each high-risk drug or interaction, suggest alternative
        for interaction in interactions:
            if interaction.severity in (InteractionSeverity.CRITICAL, InteractionSeverity.MAJOR):
                # Suggest removing one drug or replacing
                alternatives.append({
                    "type": "replace_drug",
                    "problematic_drugs": [interaction.drug_a, interaction.drug_b],
                    "reason": f"{interaction.interaction_type.value}: {interaction.mechanism}",
                    "suggestion": f"Consider alternative to {interaction.drug_b} (or {interaction.drug_a})",
                    "benefit": f"Eliminates {interaction.severity.value} interaction",
                })

        # For high individual risk drugs
        for assess in drug_assessments:
            if assess.risk.risk_score and assess.risk.risk_score >= 80:
                alternatives.append({
                    "type": "replace_drug",
                    "problematic_drugs": [assess.drug],
                    "reason": f"High individual risk (score: {assess.risk.risk_score:.0f}/100)",
                    "suggestion": assess.recommendation.text,
                    "benefit": "Reduces overall polypharmacy risk",
                })

        return alternatives

    def _get_risk_level(self, risk_score: float) -> str:
        """Map risk score to categorical level."""
        if risk_score >= 90:
            return "critical"
        elif risk_score >= 70:
            return "major"
        elif risk_score >= 40:
            return "moderate"
        elif risk_score >= 20:
            return "minor"
        else:
            return "none"


# ============================================================================
# Interaction Matrix Builder
# ============================================================================

class InteractionMatrix:
    """Build and query drug-drug-gene interaction matrix."""

    def __init__(self, interactions: List[DrugDrugInteraction]):
        self.interactions = interactions
        self._build_matrix()

    def _build_matrix(self):
        """Build indexed matrix for fast lookup."""
        self.by_drug_pair: Dict[Tuple[str, str], List[DrugDrugInteraction]] = {}
        self.by_gene: Dict[str, List[DrugDrugInteraction]] = {}
        self.by_severity: Dict[str, List[DrugDrugInteraction]] = {}

        for interaction in self.interactions:
            # Index by drug pair (both orderings)
            pair1 = (interaction.drug_a.lower(), interaction.drug_b.lower())
            pair2 = (interaction.drug_b.lower(), interaction.drug_a.lower())

            if pair1 not in self.by_drug_pair:
                self.by_drug_pair[pair1] = []
            self.by_drug_pair[pair1].append(interaction)

            if pair2 not in self.by_drug_pair:
                self.by_drug_pair[pair2] = []
            self.by_drug_pair[pair2].append(interaction)

            # Index by gene
            if interaction.gene not in self.by_gene:
                self.by_gene[interaction.gene] = []
            self.by_gene[interaction.gene].append(interaction)

            # Index by severity
            severity = interaction.severity.value
            if severity not in self.by_severity:
                self.by_severity[severity] = []
            self.by_severity[severity].append(interaction)

    def get_interactions_for_pair(
        self,
        drug_a: str,
        drug_b: str
    ) -> List[DrugDrugInteraction]:
        """Get all interactions for a drug pair."""
        pair = (drug_a.lower(), drug_b.lower())
        return self.by_drug_pair.get(pair, [])

    def get_interactions_for_gene(self, gene: str) -> List[DrugDrugInteraction]:
        """Get all interactions involving a gene."""
        return self.by_gene.get(gene, [])

    def get_interactions_by_severity(
        self,
        severity: InteractionSeverity
    ) -> List[DrugDrugInteraction]:
        """Get all interactions of a given severity."""
        return self.by_severity.get(severity.value, [])


# ============================================================================
# Convenience Functions
# ============================================================================

def analyze_multi_drug_risk(
    drug_assessments: List[DrugAssessment]
) -> MultiDrugRiskAssessment:
    """
    Convenience function for multi-drug risk analysis.

    Args:
        drug_assessments: List of individual drug assessments

    Returns:
        MultiDrugRiskAssessment
    """
    analyzer = MultiDrugRiskAnalyzer()
    return analyzer.analyze_multi_drug_risk(drug_assessments)


def get_interaction_matrix() -> InteractionMatrix:
    """Get interaction matrix for querying."""
    return InteractionMatrix(KNOWN_INTERACTIONS)
