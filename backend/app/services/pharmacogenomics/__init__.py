"""
Pharmacogenomics Service

CPIC-aligned pharmacogenomic decision engine for drug risk assessment.
Provides deterministic, rule-based diplotype resolution and clinical recommendations.
"""

from .models import (
    VariantCall,
    GenotypeData,
    RiskAssessment,
    ClinicalRecommendation,
    DiplotypeResult,
    PatientProfile,
    DrugAssessment,
    IndeterminateReason
)
from .cpic_loader import get_cpic_loader, reload_cpic_data
from .phenotype_mapper import DiplotypeResolver, PhenotypeMapper
from .risk_engine import RiskEngine, create_risk_engine
from .config import (
    get_config,
    update_config,
    load_config_from_file,
    save_config_to_file,
    get_confidence_penalties,
    get_diplotype_config,
    get_activity_scores
)
from .population_data import (
    get_population_frequencies,
    Population,
    PopulationFrequencyData
)

__all__ = [
    # Models
    'VariantCall',
    'GenotypeData',
    'RiskAssessment',
    'ClinicalRecommendation',
    'DiplotypeResult',
    'PatientProfile',
    'DrugAssessment',

    # Loader
    'get_cpic_loader',
    'reload_cpic_data',

    # Phenotype Mapping
    'DiplotypeResolver',
    'PhenotypeMapper',

    # Risk Engine
    'RiskEngine',
    'create_risk_engine',
]
