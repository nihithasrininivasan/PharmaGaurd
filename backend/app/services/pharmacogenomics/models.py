"""
Internal data models for pharmacogenomics service.
These models represent the intermediate data structures used for processing
VCF data and resolving diplotypes/phenotypes.
"""

from pydantic import BaseModel, Field
from typing import List, Optional, Dict
from enum import Enum
from datetime import datetime

class IndeterminateReason(str, Enum):
    """Specific reasons for indeterminate diplotype calls."""
    NONE = "none"  # Not indeterminate
    NO_COVERAGE = "no_coverage"  # Key positions missing coverage
    AMBIGUOUS = "ambiguous"  # Multiple valid configurations
    NOVEL_VARIANTS = "novel_variants"  # Variants don't match known alleles
    PARTIAL_MATCH = "partial_match"  # Incomplete match to allele definition
    LOW_QUALITY = "low_quality"  # Low quality variant calls
    UNSUPPORTED_GENE = "unsupported_gene"  # Gene not in database


class VariantCall(BaseModel):
    """Represents a single variant call from VCF data."""
    chrom: str = Field(..., description="Chromosome identifier")
    pos: int = Field(..., description="Position on chromosome")
    rsid: Optional[str] = Field(None, description="dbSNP reference ID")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    zygosity: str = Field(..., description="Zygosity: HET, HOM_REF, or HOM_ALT")
    quality: float = Field(default=0.0, description="Variant quality score")
    filter: Optional[str] = Field(None, description="Filter status (PASS, etc.)")
    ad: Optional[List[int]] = Field(None, description="Allele depth (reference, alternate)")
    star_allele: Optional[str] = Field(None, description="Star allele annotation from VCF INFO")

    # Phasing information
    phased: bool = Field(default=False, description="Whether this variant is phased")
    phase_set: Optional[str] = Field(None, description="Phase set ID (PS tag from VCF)")
    haplotype: Optional[int] = Field(None, description="Haplotype assignment (0 or 1)")

    def variant_key(self) -> str:
        """Generate a unique key for this variant: POS:REF:ALT"""
        return f"{self.pos}:{self.ref}:{self.alt}"


class GenotypeData(BaseModel):
    """Genotype data for a specific gene from a single sample."""
    sample_id: str = Field(..., description="Sample identifier")
    gene_symbol: str = Field(..., description="Gene symbol (e.g., CYP2D6)")
    variants: List[VariantCall] = Field(default_factory=list, description="List of variant calls")
    coverage_mean: Optional[float] = Field(None, description="Mean coverage depth across gene region")
    covered_positions: List[int] = Field(default_factory=list, description="Specific positions with adequate coverage")
    genome_build: Optional[str] = Field(None, description="Genome build (GRCh37, GRCh38, or None/Unknown)")

    def get_variant_keys(self) -> List[str]:
        """Get all variant keys for this genotype."""
        return [v.variant_key() for v in self.variants]


class VariantAnnotation(BaseModel):
    """Functional annotation for a specific variant from PharmGKB clinicalVariants.tsv."""
    variant_id: str = Field(..., description="rsID or star allele identifier")
    gene: str = Field(..., description="Gene symbol")
    annotation_type: str = Field(..., description="Type: Metabolism/PK, Toxicity, Efficacy, Dosage")
    evidence_level: str = Field(..., description="Evidence level: 1A, 1B, 2A, 2B, 3, 4")
    associated_chemicals: List[str] = Field(default_factory=list, description="Associated drug names")
    associated_phenotypes: List[str] = Field(default_factory=list, description="Associated phenotypes")


class GeneDrugConfirmation(BaseModel):
    """Deterministic gene-drug pair confirmation from PharmGKB relationships.tsv."""
    gene: str = Field(..., description="Gene symbol")
    drug: str = Field(..., description="Drug name")
    confirmed: bool = Field(..., description="Whether gene-drug pair exists in dataset")
    evidence_types: List[str] = Field(default_factory=list, description="Evidence types found")
    association: str = Field("not found", description="Strongest association: associated/ambiguous/not associated")
    source_pmids: List[str] = Field(default_factory=list, description="Supporting PubMed IDs")


class ClinicalAnnotationLink(BaseModel):
    """Clinical annotation link from PharmGKB relationships.tsv."""
    annotation_id: str = Field(..., description="PharmGKB annotation ID")
    gene: str = Field(..., description="Gene symbol")
    drug: str = Field(..., description="Drug name")
    evidence_type: str = Field(..., description="Evidence type")
    association: str = Field(..., description="Association direction")
    pmids: List[str] = Field(default_factory=list, description="Supporting PubMed IDs")


class EvidenceLevel(BaseModel):
    """Evidence level scoring for a gene-drug pair."""
    level: str = Field(..., description="Evidence level: 1A, 1B, 2A, 2B, 3, 4, none")
    confidence_weight: float = Field(..., ge=0.0, le=1.0, description="Confidence weight multiplier")
    allows_automated_recommendation: bool = Field(..., description="Whether automated CDS is permitted")


class RiskAssessment(BaseModel):
    """Risk assessment result for a drug-gene interaction."""
    risk_label: str = Field(..., description="Risk classification label")
    confidence_score: float = Field(..., ge=0.0, le=1.0, description="Confidence score (0-1)")
    severity: str = Field(..., description="Severity: none, low, moderate, high, critical, undetermined")
    confidence_breakdown: Optional[Dict] = Field(None, description="5-layer confidence component scores")
    risk_score: Optional[float] = Field(None, ge=0.0, le=100.0, description="Numeric risk score (0-100)")
    risk_level: Optional[str] = Field(None, description="Risk level: none/low/moderate/high/critical")
    gene_drug_confirmation: Optional[Dict] = Field(None, description="PharmGKB gene-drug pair confirmation")
    evidence_level: Optional[Dict] = Field(None, description="Evidence level scoring")
    clinical_annotations: Optional[List[Dict]] = Field(None, description="Clinical annotation links (deduplicated)")
    automation_status: Optional[Dict] = Field(None, description="Automation gate status: {allowed, blocked_reasons}")


class ClinicalRecommendation(BaseModel):
    """Clinical recommendation based on pharmacogenomic data."""
    text: str = Field(..., description="Recommendation text")
    implication: str = Field(..., description="Clinical implication")
    recommendation_url: Optional[str] = Field(None, description="URL to CPIC guideline")


class DiplotypeResult(BaseModel):
    """Result of diplotype resolution for a gene."""
    gene: str = Field(..., description="Gene symbol")
    diplotype: str = Field(..., description="Diplotype (e.g., *1/*2)")
    phenotype: str = Field(..., description="Phenotype (e.g., IM, NM, PM)")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Confidence score")
    is_indeterminate: bool = Field(False, description="Whether diplotype is indeterminate")
    indeterminate_reason: IndeterminateReason = Field(
        default=IndeterminateReason.NONE,
        description="Specific reason for indeterminate call"
    )
    notes: Optional[str] = Field(None, description="Additional notes about resolution")
    phased: bool = Field(default=False, description="Whether diplotype is based on phased data")
    confidence_breakdown: Optional[Dict] = Field(None, description="Itemised confidence component scores")


class PatientProfile(BaseModel):
    """Patient pharmacogenomic profile containing diplotypes for multiple genes."""
    sample_id: str = Field(..., description="Patient/sample identifier")
    diplotypes: Dict[str, DiplotypeResult] = Field(default_factory=dict, description="Diplotype results keyed by gene")


class DrugAssessment(BaseModel):
    """Complete assessment for a specific drug."""
    drug: str = Field(..., description="Drug name")
    gene: str = Field(..., description="Associated gene")
    diplotype: str = Field(..., description="Patient's diplotype")
    phenotype: str = Field(..., description="Patient's phenotype")
    risk: RiskAssessment = Field(..., description="Risk assessment")
    recommendation: ClinicalRecommendation = Field(..., description="Clinical recommendation")


class PharmacogenomicProfile(BaseModel):
    """Profile details for the report."""
    primary_gene: str = Field(..., description="Gene symbol")
    diplotype: str = Field(..., description="Diplotype result")
    phenotype: str = Field(..., description="Phenotype description")
    detected_variants: List[VariantCall] = Field(..., description="List of variants found")
    variant_annotations: List[VariantAnnotation] = Field(
        default_factory=list, description="PharmGKB variant annotations"
    )


class PharmacogenomicReport(BaseModel):
    """Final report structure matching user output requirements."""
    patient_id: str = Field(..., description="Patient Identifier")
    drug: str = Field(..., description="Drug Name")
    timestamp: str = Field(..., description="ISO8601 Timestamp")
    risk_assessment: RiskAssessment = Field(..., description="Risk details")
    pharmacogenomic_profile: PharmacogenomicProfile = Field(..., description="PGx Profile")
    clinical_recommendation: ClinicalRecommendation = Field(..., description="Clinical Recommendation")

