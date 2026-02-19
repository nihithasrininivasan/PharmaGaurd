"""
Configuration for pharmacogenomics service.
Centralizes tunable parameters for diplotype resolution and confidence scoring.
"""

from typing import Dict
from pydantic import BaseModel, Field


class ConfidencePenalties(BaseModel):
    """Confidence penalty factors for various scenarios."""

    # Coverage penalties
    missing_key_position: float = Field(
        default=0.8,
        ge=0.0,
        le=1.0,
        description="Penalty multiplier per missing key position (0.8 = 20% reduction)"
    )

    # Ambiguity penalties
    unphased_heterozygote: float = Field(
        default=0.9,
        ge=0.0,
        le=1.0,
        description="Penalty for unphased heterozygous variants (phase uncertainty)"
    )

    partial_allele_match: float = Field(
        default=0.7,
        ge=0.0,
        le=1.0,
        description="Penalty for partial matches to allele definitions"
    )

    indeterminate_call: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Penalty for indeterminate diplotype calls"
    )

    rare_unknown_allele: float = Field(
        default=0.7,
        ge=0.0,
        le=1.0,
        description="Penalty for rare or unknown alleles"
    )

    no_coverage_data: float = Field(
        default=0.9,
        ge=0.0,
        le=1.0,
        description="Small penalty when no coverage data is provided"
    )


class DiplotypeResolutionConfig(BaseModel):
    """Configuration for diplotype resolution algorithm."""

    # Scoring thresholds
    homozygous_score_threshold: float = Field(
        default=2.0,
        ge=0.0,
        description="Minimum score to call homozygous variant (HOM requires 2.0)"
    )

    heterozygous_score_threshold: float = Field(
        default=1.0,
        ge=0.0,
        description="Minimum score to call heterozygous variant"
    )

    compound_het_min_score: float = Field(
        default=1.0,
        ge=0.0,
        description="Minimum score for each allele in compound heterozygote"
    )

    # Allele completeness requirement
    require_complete_match: bool = Field(
        default=False,
        description="If True, require all defining variants to be present for allele call"
    )

    completeness_threshold: float = Field(
        default=0.8,
        ge=0.0,
        le=1.0,
        description="Minimum fraction of defining variants required (if not requiring complete match)"
    )

    # Phase assumption
    default_phase_assumption: str = Field(
        default="trans",
        description="Default phase assumption for unphased variants: 'trans' or 'cis'"
    )


class ActivityScoreConfig(BaseModel):
    """Configuration for activity score-based phenotype mapping."""

    # Activity score to phenotype thresholds
    # Based on CPIC standard ranges

    poor_metabolizer_max: float = Field(
        default=0.5,
        description="Maximum activity score for Poor Metabolizer (PM)"
    )

    intermediate_metabolizer_max: float = Field(
        default=1.5,
        description="Maximum activity score for Intermediate Metabolizer (IM)"
    )

    normal_metabolizer_max: float = Field(
        default=2.5,
        description="Maximum activity score for Normal Metabolizer (NM)"
    )

    # Above normal_metabolizer_max is UM (Ultra-rapid Metabolizer)

    # Gene-specific activity scores (can be extended)
    gene_specific_scores: Dict[str, Dict[str, float]] = Field(
        default_factory=lambda: {
            "CYP2D6": {
                "*1": 1.0,
                "*2": 1.0,
                "*3": 0.0,
                "*4": 0.0,
                "*5": 0.0,  # Gene deletion
                "*6": 0.0,
                "*7": 0.0,
                "*8": 0.0,
                "*9": 0.5,
                "*10": 0.5,
                "*17": 0.5,
                "*29": 0.5,
                "*41": 0.5,
                # Duplications
                "*1x2": 2.0,
                "*2x2": 2.0,
            },
            "CYP2C19": {
                "*1": 1.0,
                "*2": 0.0,
                "*3": 0.0,
                "*17": 1.5,  # Increased function
            },
            "CYP2C9": {
                "*1": 1.0,
                "*2": 0.5,
                "*3": 0.0,
            },
            "TPMT": {
                "*1": 1.0,
                "*2": 0.0,
                "*3A": 0.0,
                "*3C": 0.0,
            }
        },
        description="Gene-specific allele activity scores"
    )


class PharmacogenomicsConfig(BaseModel):
    """Main configuration for pharmacogenomics service."""

    confidence_penalties: ConfidencePenalties = Field(
        default_factory=ConfidencePenalties,
        description="Confidence penalty configuration"
    )

    diplotype_resolution: DiplotypeResolutionConfig = Field(
        default_factory=DiplotypeResolutionConfig,
        description="Diplotype resolution configuration"
    )

    activity_scores: ActivityScoreConfig = Field(
        default_factory=ActivityScoreConfig,
        description="Activity score configuration"
    )

    # Data paths
    cpic_cache_path: str = Field(
        default="data/cpic_cache.json",
        description="Path to CPIC cache file (relative to backend root)"
    )

    cpic_data_dir: str = Field(
        default="data/cpic",
        description="Path to CPIC raw data directory"
    )

    # Auto-ETL settings
    auto_run_etl: bool = Field(
        default=True,
        description="Automatically run ETL if cache is missing but data files exist"
    )

    # Logging
    verbose_logging: bool = Field(
        default=True,
        description="Enable verbose logging for debugging"
    )


# Global configuration instance
_config: PharmacogenomicsConfig = PharmacogenomicsConfig()


def get_config() -> PharmacogenomicsConfig:
    """Get the global configuration instance."""
    return _config


def update_config(**kwargs):
    """Update configuration parameters."""
    global _config
    current_dict = _config.model_dump()

    # Update nested parameters
    for key, value in kwargs.items():
        if '.' in key:
            # Handle nested keys like 'confidence_penalties.missing_key_position'
            parts = key.split('.')
            current = current_dict
            for part in parts[:-1]:
                current = current[part]
            current[parts[-1]] = value
        else:
            current_dict[key] = value

    _config = PharmacogenomicsConfig(**current_dict)
    return _config


def load_config_from_file(filepath: str):
    """Load configuration from a JSON file."""
    import json
    global _config

    with open(filepath, 'r') as f:
        config_dict = json.load(f)

    _config = PharmacogenomicsConfig(**config_dict)
    return _config


def save_config_to_file(filepath: str):
    """Save current configuration to a JSON file."""
    import json

    with open(filepath, 'w') as f:
        json.dump(_config.model_dump(), f, indent=2)


# Convenience accessors
def get_confidence_penalties() -> ConfidencePenalties:
    """Get confidence penalties configuration."""
    return _config.confidence_penalties


def get_diplotype_config() -> DiplotypeResolutionConfig:
    """Get diplotype resolution configuration."""
    return _config.diplotype_resolution


def get_activity_scores() -> ActivityScoreConfig:
    """Get activity score configuration."""
    return _config.activity_scores
