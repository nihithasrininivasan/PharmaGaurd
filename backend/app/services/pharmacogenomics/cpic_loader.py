"""
CPIC Data Loader - Runtime singleton for loading and accessing CPIC data.
Provides optimized lookup structures for diplotype resolution and risk assessment.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Set
from functools import lru_cache


# ============================================================================
# Hardcoded CPIC Fallback Activity Scores
# Used when cpic_cache.json has empty activity_scores for a gene.
# Sources: CPIC guidelines for each gene.
# ============================================================================
_FALLBACK_ACTIVITY = {
    "CYP2D6": {
        "*1": 1.0, "*1x2": 2.0, "*1xN": 2.0, "*2": 1.0, "*2x2": 2.0,
        "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0, "*7": 0.0, "*8": 0.0,
        "*9": 0.5, "*10": 0.5, "*14": 0.0, "*17": 0.5, "*29": 0.5,
        "*35": 1.0, "*36": 0.0, "*41": 0.5, "*2A": 1.0,
    },
    "CYP2C19": {
        "*1": 1.0, "*2": 0.0, "*3": 0.0, "*4": 0.0, "*5": 0.0,
        "*6": 0.0, "*7": 0.0, "*8": 0.0, "*9": 0.5, "*10": 0.5,
        "*17": 1.5, "*27": 0.5, "*35": 0.5,
    },
    "CYP2C9": {
        "*1": 1.0, "*2": 0.5, "*3": 0.0, "*4": 0.0, "*5": 0.0,
        "*6": 0.0, "*8": 0.5, "*11": 0.5, "*12": 0.5, "*13": 0.0,
        "*14": 0.5,
    },
    "SLCO1B1": {
        "*1a": 1.0, "*1b": 1.0, "*1B": 1.0, "*1": 1.0,
        "*2": 0.5, "*3": 0.5, "*4": 0.5, "*5": 0.0, "*6": 0.0,
        "*9": 0.5, "*10": 0.5, "*14": 0.5, "*15": 0.0, "*16": 0.0,
        "*17": 0.5, "*19": 0.0, "*20": 0.5, "*22": 0.5, "*24": 0.0,
        "*25": 0.5, "*28": 0.5, "*31": 0.5, "*37": 0.0, "*38": 0.5,
        "*45": 0.5,
    },
    "TPMT": {
        "*1": 1.0, "*2": 0.0, "*3A": 0.0, "*3B": 0.0, "*3C": 0.0,
        "*3D": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0, "*7": 0.0,
        "*8": 0.5, "*9": 0.5, "*10": 0.5, "*11": 0.5, "*12": 0.5,
        "*18": 0.0, "*19": 0.5, "*20": 0.0, "*21": 0.5, "*22": 0.5,
        "*23": 0.5, "*24": 0.5, "*25": 0.5, "*26": 0.0, "*28": 0.0,
        "*29": 0.0, "*38": 0.5,
    },
    "DPYD": {
        "*1": 1.0, "*2A": 0.0, "*3": 0.5, "*4": 0.5, "*5": 0.5,
        "*6": 0.5, "*7": 0.0, "*8": 0.5, "*9A": 0.5, "*9B": 0.5,
        "*10": 0.5, "*11": 0.0, "*12": 0.0, "*13": 0.0,
    },
}


class CPICDataLoader:
    """
    Singleton loader for CPIC pharmacogenomic data.
    Loads data from cpic_cache.json at startup and provides efficient lookup methods.
    """

    _instance: Optional['CPICDataLoader'] = None
    _initialized: bool = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not CPICDataLoader._initialized:
            self._data: Dict = {}
            self._genes: Dict = {}
            self._drugs: Dict = {}
            self._load_data()
            CPICDataLoader._initialized = True

    def _load_data(self):
        """Load CPIC data from cache file."""
        from .config import get_config

        config = get_config()

        # Determine path to cpic_cache.json
        backend_dir = Path(__file__).parent.parent.parent.parent
        cache_file = backend_dir / config.cpic_cache_path

        # Check for cache freshness
        cpic_data_dir = backend_dir / config.cpic_data_dir
        if cache_file.exists() and cpic_data_dir.exists():
             cache_mtime = cache_file.stat().st_mtime
             
             # Check if any source file is newer than cache
             needs_refresh = False
             for source_file in cpic_data_dir.rglob("*.xlsx"):
                 if source_file.stat().st_mtime > cache_mtime:
                     print(f"Detected update in {source_file.name}. Refreshing cache...")
                     needs_refresh = True
                     break
             
             if needs_refresh:
                 self._run_etl(cpic_data_dir, cache_file)

        if not cache_file.exists():
            # Check if auto-ETL is enabled
            if config.auto_run_etl:
                if cpic_data_dir.exists():
                    print(f"CPIC cache not found. Running ETL from {cpic_data_dir}...")
                    self._run_etl(cpic_data_dir, cache_file)
                else:
                    raise FileNotFoundError(
                        f"CPIC cache file not found at {cache_file} and "
                        f"CPIC data directory not found at {cpic_data_dir}. "
                        "Please provide CPIC data files."
                    )
            else:
                raise FileNotFoundError(
                    f"CPIC cache file not found at {cache_file}. "
                    "Please run cpic_etl.py to generate the cache, or enable auto_run_etl in config."
                )

        with open(cache_file, 'r') as f:
            self._data = json.load(f)

        self._genes = self._data.get('genes', {})
        self._drugs = self._data.get('drugs', {})

        # Manual Patch for Warfarin (CYP2C9) if missing from cache
        if 'warfarin' not in self._drugs:
             self._drugs['warfarin'] = {
                 'gene': 'CYP2C9',
                 'id': 'warfarin',
                 'name': 'Warfarin',
                 'recommendations': {
                     'Normal Metabolizer': {
                         'risk': 'Standard dosing recommended',
                         'implication': 'Normal CYP2C9 metabolism. Dose based on clinical guidelines.',
                         'severity': 'none',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/'
                     },
                     'Intermediate Metabolizer': {
                         'risk': 'Consider lower dose',
                         'implication': 'Reduced CYP2C9 metabolism. Increased risk of bleeding.',
                         'severity': 'moderate',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/'
                     },
                     'Poor Metabolizer': {
                         'risk': 'Consider lower dose',
                         'implication': 'Poor CYP2C9 metabolism. Increased risk of bleeding.',
                         'severity': 'high',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/'
                     }
                 }
             }

        # Manual Patch for Azathioprine (TPMT) if missing from cache
        if 'azathioprine' not in self._drugs:
             self._drugs['azathioprine'] = {
                 'gene': 'TPMT',
                 'id': 'azathioprine',
                 'name': 'Azathioprine',
                 'recommendations': {
                     'Normal Metabolizer': {
                         'risk': 'Standard dosing recommended',
                         'implication': 'Normal TPMT activity. Standard dose expected to be effective and safe.',
                         'severity': 'none',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/'
                     },
                     'Intermediate Metabolizer': {
                         'risk': 'Reduce dose by 30-70%',
                         'implication': 'Reduced TPMT activity. Higher risk of myelosuppression at standard doses.',
                         'severity': 'moderate',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/'
                     },
                     'Poor Metabolizer': {
                         'risk': 'Avoid or drastically reduce dose',
                         'implication': 'Absent TPMT activity. Very high risk of life-threatening myelosuppression.',
                         'severity': 'critical',
                         'url': 'https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/'
                     }
                 }
             }

        if config.verbose_logging:
            print(f"CPIC Data Loader initialized: {len(self._genes)} genes, {len(self._drugs)} drugs")

    def _run_etl(self, cpic_data_dir: Path, output_file: Path):
        """Run ETL process to generate cache from CPIC data files."""
        try:
            # Import and run ETL
            import sys
            sys.path.insert(0, str(Path(__file__).parent.parent.parent))

            from utils.cpic_etl import CPICDataProcessor, add_manual_phenotype_mappings, add_drug_guidelines

            processor = CPICDataProcessor(str(cpic_data_dir), str(output_file))
            processor.process_all_genes()
            add_manual_phenotype_mappings(processor)
            add_drug_guidelines(processor)
            processor.save()

            print(f"✓ ETL completed successfully")
        except Exception as e:
            raise RuntimeError(f"Failed to run ETL: {e}")

    # ===== Gene Data Access =====

    def get_gene_data(self, gene: str) -> Optional[Dict]:
        """Get all data for a specific gene."""
        return self._genes.get(gene)

    def get_allele_definitions(self, gene: str) -> Dict[str, List[str]]:
        """
        Get allele definitions for a gene.
        Returns: {allele_name: [variant_keys], ...}
        Example: {"*2": ["42126611:C:G"], "*4": [...], ...}
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('allele_definitions', {})

    def get_variant_to_allele_map(self, gene: str) -> Dict[str, List[str]]:
        """
        Get variant-to-allele mapping for a gene.
        Returns: {variant_key: [allele_names], ...}
        Example: {"42126611:C:G": ["*2", "*8"], ...}
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('variant_to_allele', {})

    def get_phenotype_map(self, gene: str) -> Dict[str, str]:
        """
        Get diplotype-to-phenotype mapping for a gene.
        Returns: {diplotype: phenotype, ...}
        Example: {"*1/*1": "NM", "*1/*2": "IM", ...}
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('phenotype_map', {})

    def get_key_positions(self, gene: str) -> List[int]:
        """
        Get key genomic positions for a gene (positions that define common alleles).
        These are used for coverage confidence scoring.
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return []

        positions_data = gene_data.get('positions', {})
        return [pos_info['pos'] for pos_info in positions_data.values() if 'pos' in pos_info]

    def get_supported_genes(self) -> List[str]:
        """Get list of all supported gene symbols."""
        return list(self._genes.keys())

    # ===== Drug Data Access =====

    def get_drug_data(self, drug: str) -> Optional[Dict]:
        """Get all data for a specific drug."""
        return self._drugs.get(drug.lower())

    def get_drug_gene(self, drug: str) -> Optional[str]:
        """Get the primary gene associated with a drug."""
        drug_data = self.get_drug_data(drug)
        if not drug_data:
            return None
        return drug_data.get('gene')

    def get_drug_recommendations(self, drug: str) -> Dict[str, Dict]:
        """
        Get clinical recommendations for a drug by phenotype.
        Returns: {phenotype: {risk, severity, implication, url}, ...}
        """
        drug_data = self.get_drug_data(drug)
        if not drug_data:
            return {}
        return drug_data.get('recommendations', {})


    # Phenotype key aliases so both "NM" and "Normal Metabolizer" etc. match cache keys
    _PHENOTYPE_ALIASES: Dict[str, List[str]] = {
        "NM": ["NM", "Normal Metabolizer", "Extensive Metabolizer", "Normal Function", "Normal Activity"],
        "Normal Metabolizer": ["Normal Metabolizer", "NM", "Extensive Metabolizer"],
        "IM": ["IM", "Intermediate Metabolizer", "Intermediate Activity", "Intermediate Function"],
        "Intermediate Metabolizer": ["Intermediate Metabolizer", "IM"],
        "PM": ["PM", "Poor Metabolizer", "Poor Function", "Low Activity"],
        "Poor Metabolizer": ["Poor Metabolizer", "PM"],
        "UM": ["UM", "Ultra-rapid Metabolizer"],
        "Ultra-rapid Metabolizer": ["Ultra-rapid Metabolizer", "UM"],
        "RM": ["RM", "Rapid Metabolizer"],
        "Rapid Metabolizer": ["Rapid Metabolizer", "RM"],
    }

    # =====================================================================
    # Comprehensive CPIC-based manual recommendation patches
    # Keyed by (DRUG_UPPER, phenotype_string)
    # =====================================================================
    _MANUAL_RECOMMENDATIONS = {
        # ── CODEINE (CYP2D6) ──────────────────
        ("CODEINE", "Poor Metabolizer"): {
            "risk_label": "Avoid – Ineffective",
            "confidence_score": 0.95,
            "severity": "critical",
            "risk": "Avoid codeine (Ineffective)",
            "text": "Avoid codeine. Poor CYP2D6 metabolizers produce insufficient morphine; use an alternative analgesic.",
            "implication": "Inadequate pain relief; ineffective therapy.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        ("CODEINE", "Intermediate Metabolizer"): {
            "risk_label": "Use label dosage – monitor",
            "confidence_score": 0.85,
            "severity": "low",
            "risk": "Monitor for reduced efficacy",
            "text": "Initiate therapy with label-recommended dosage. Monitor for reduced efficacy as CYP2D6 activity is reduced.",
            "implication": "Reduced morphine production; potential for decreased analgesic effect.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        ("CODEINE", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use codeine per standard dosing guidelines. Normal CYP2D6 metabolism expected.",
            "implication": "Expected normal morphine formation and analgesic effect.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        ("CODEINE", "Ultrarapid Metabolizer"): {
            "risk_label": "Avoid – Toxicity Risk",
            "confidence_score": 0.95,
            "severity": "critical",
            "risk": "Avoid codeine (Toxicity risk)",
            "text": "Avoid codeine. Ultra-rapid CYP2D6 metabolizers produce excessive morphine; risk of respiratory depression.",
            "implication": "Excessive morphine formation; life-threatening respiratory depression possible.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },

        # ── CLOPIDOGREL (CYP2C19) ─────────────
        ("CLOPIDOGREL", "Poor Metabolizer"): {
            "risk_label": "Avoid – use alternative antiplatelet",
            "confidence_score": 0.97,
            "severity": "high",
            "risk": "Avoid clopidogrel",
            "text": "Avoid clopidogrel. Use an alternative antiplatelet agent (e.g., prasugrel, ticagrelor).",
            "implication": "Minimal clopidogrel activation; high risk of cardiovascular events.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        ("CLOPIDOGREL", "Intermediate Metabolizer"): {
            "risk_label": "Consider alternative antiplatelet",
            "confidence_score": 0.90,
            "severity": "moderate",
            "risk": "Consider alternative antiplatelet",
            "text": "Consider an alternative antiplatelet agent (e.g., prasugrel, ticagrelor). If clopidogrel is used, be aware of reduced platelet inhibition.",
            "implication": "Reduced clopidogrel activation; diminished antiplatelet response.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        ("CLOPIDOGREL", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use clopidogrel per standard dosing guidelines. Normal CYP2C19 metabolism expected.",
            "implication": "Normal clopidogrel activation; expected therapeutic antiplatelet response.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        ("CLOPIDOGREL", "Rapid Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.85,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use clopidogrel per standard dosing guidelines. Rapid metabolism may increase activation.",
            "implication": "Increased clopidogrel activation; standard therapy expected to be effective.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        ("CLOPIDOGREL", "Ultrarapid Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.85,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use clopidogrel per standard dosing guidelines. Ultra-rapid metabolism not associated with adverse effects.",
            "implication": "Increased clopidogrel activation; standard dosing appropriate.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },

        # ── WARFARIN (CYP2C9) ─────────────────
        ("WARFARIN", "Poor Metabolizer"): {
            "risk_label": "Significantly reduce dose",
            "confidence_score": 0.92,
            "severity": "high",
            "risk": "Significantly reduce warfarin dose",
            "text": "Reduce warfarin dose significantly. Poor CYP2C9 metabolizers have greatly increased bleeding risk.",
            "implication": "Significantly reduced warfarin clearance; high risk of over-anticoagulation and bleeding.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/"
        },
        ("WARFARIN", "Intermediate Metabolizer"): {
            "risk_label": "Consider lower dose",
            "confidence_score": 0.88,
            "severity": "moderate",
            "risk": "Consider lower warfarin dose",
            "text": "Consider a lower initial warfarin dose. Intermediate CYP2C9 metabolism increases bleeding risk.",
            "implication": "Reduced warfarin clearance; increased sensitivity and bleeding risk.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/"
        },
        ("WARFARIN", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use standard warfarin dosing algorithm. Normal CYP2C9 metabolism expected.",
            "implication": "Normal warfarin metabolism; standard dose adjustment per INR monitoring.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9/"
        },

        # ── SIMVASTATIN (SLCO1B1) ─────────────
        ("SIMVASTATIN", "Poor Function"): {
            "risk_label": "Avoid or use low dose",
            "confidence_score": 0.92,
            "severity": "high",
            "risk": "Avoid simvastatin or use ≤20mg",
            "text": "Avoid simvastatin or use a low dose (≤20 mg). Consider an alternative statin (e.g., pravastatin, rosuvastatin).",
            "implication": "Significantly increased simvastatin exposure; high risk of myopathy/rhabdomyolysis.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/"
        },
        ("SIMVASTATIN", "Decreased Function"): {
            "risk_label": "Prescribe ≤20mg or alternative",
            "confidence_score": 0.88,
            "severity": "moderate",
            "risk": "Use ≤20mg simvastatin or alternative statin",
            "text": "Prescribe simvastatin ≤20 mg or consider an alternative statin. Increased myopathy risk.",
            "implication": "Increased simvastatin plasma levels; elevated risk of myopathy.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/"
        },
        ("SIMVASTATIN", "Normal Function"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use simvastatin per standard dosing guidelines. Normal SLCO1B1 transporter function expected.",
            "implication": "Normal simvastatin hepatic uptake; standard myopathy risk.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/"
        },

        # ── FLUOROURACIL (DPYD) ───────────────
        ("FLUOROURACIL", "Poor Metabolizer"): {
            "risk_label": "Avoid – life-threatening toxicity",
            "confidence_score": 0.97,
            "severity": "critical",
            "risk": "Avoid fluorouracil",
            "text": "Avoid fluoropyrimidines. Complete DPD deficiency causes life-threatening toxicity.",
            "implication": "No DPD activity; fatal toxicity risk with standard doses.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/"
        },
        ("FLUOROURACIL", "Intermediate Metabolizer"): {
            "risk_label": "Reduce dose by 50%",
            "confidence_score": 0.92,
            "severity": "high",
            "risk": "Reduce fluorouracil dose by 50%",
            "text": "Reduce fluoropyrimidine dose by 50%. Partial DPD deficiency increases severe toxicity risk.",
            "implication": "Reduced DPD activity; increased risk of severe/fatal fluoropyrimidine toxicity.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/"
        },
        ("FLUOROURACIL", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use fluoropyrimidines per standard dosing. Normal DPD activity expected.",
            "implication": "Normal DPD metabolism; standard toxicity risk.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/"
        },

        # ── AZATHIOPRINE (TPMT) ───────────────
        ("AZATHIOPRINE", "Poor Metabolizer"): {
            "risk_label": "Avoid or drastically reduce dose",
            "confidence_score": 0.95,
            "severity": "critical",
            "risk": "Avoid azathioprine or use 10% of dose",
            "text": "Consider alternative agent or drastically reduce azathioprine dose (to 10% of standard). Monitor for myelosuppression.",
            "implication": "Absent TPMT activity; very high risk of life-threatening myelosuppression.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },
        ("AZATHIOPRINE", "Intermediate Metabolizer"): {
            "risk_label": "Reduce dose by 30-70%",
            "confidence_score": 0.90,
            "severity": "moderate",
            "risk": "Reduce azathioprine dose by 30-70%",
            "text": "Reduce azathioprine dose by 30-70%. Monitor closely for myelosuppression.",
            "implication": "Reduced TPMT activity; higher risk of myelosuppression at standard doses.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },
        ("AZATHIOPRINE", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use azathioprine per standard dosing guidelines. Normal TPMT activity expected.",
            "implication": "Normal TPMT activity; standard myelosuppression risk.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },

        # ── THIOGUANINE (TPMT) — alias of azathioprine pathway ──
        ("THIOGUANINE", "Poor Metabolizer"): {
            "risk_label": "Avoid or drastically reduce dose",
            "confidence_score": 0.95,
            "severity": "critical",
            "risk": "Avoid thioguanine or use 10% of dose",
            "text": "Consider alternative agent or drastically reduce thioguanine dose. Monitor for myelosuppression.",
            "implication": "Absent TPMT activity; very high risk of life-threatening myelosuppression.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },
        ("THIOGUANINE", "Intermediate Metabolizer"): {
            "risk_label": "Reduce dose by 30-70%",
            "confidence_score": 0.90,
            "severity": "moderate",
            "risk": "Reduce thioguanine dose by 30-70%",
            "text": "Reduce thioguanine dose by 30-70%. Monitor closely for myelosuppression.",
            "implication": "Reduced TPMT activity; higher risk of myelosuppression.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },
        ("THIOGUANINE", "Normal Metabolizer"): {
            "risk_label": "Standard dosing recommended",
            "confidence_score": 0.90,
            "severity": "none",
            "risk": "Standard dosing recommended",
            "text": "Use thioguanine per standard dosing guidelines. Normal TPMT activity expected.",
            "implication": "Normal TPMT activity; standard myelosuppression risk.",
            "recommendation_url": "https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/"
        },
    }

    def get_drug_recommendation_for_phenotype(self, drug: str, phenotype: str) -> Optional[Dict]:
        """Get specific recommendation for a drug-phenotype combination."""

        # 1. Check manual patches first (most reliable, CPIC-verified)
        manual = self._MANUAL_RECOMMENDATIONS.get((drug.upper(), phenotype))
        if manual:
            return dict(manual)  # Return a copy

        # 2. Try direct phenotype key lookup from CPIC cache recommendations
        recommendations = self.get_drug_recommendations(drug)
        rec = recommendations.get(phenotype)

        # 3. Try phenotype aliases
        if not rec and phenotype:
            for canonical, aliases in self._PHENOTYPE_ALIASES.items():
                if phenotype in aliases:
                    for key in aliases:
                        rec = recommendations.get(key)
                        if rec:
                            break
                    break

        # 4. Try activity-score-based key lookup
        #    Some CPIC drugs (codeine, fluorouracil) key recommendations by AS string
        if not rec and phenotype:
            gene = self.get_drug_gene(drug)
            if gene:
                total_as = self.calculate_total_activity_score(gene, self._last_diplotype if hasattr(self, '_last_diplotype') else '*1/*1')
                as_key = str(total_as) if total_as == int(total_as) else str(total_as)
                # Try exact score match
                rec = recommendations.get(as_key)
                # Try range keys (e.g., "≥5.0")
                if not rec:
                    for rkey in recommendations:
                        if rkey.startswith('≥'):
                            try:
                                threshold = float(rkey[1:])
                                if total_as >= threshold:
                                    rec = recommendations.get(rkey)
                                    break
                            except ValueError:
                                pass

        if not rec:
            return None
        return rec

    def get_supported_drugs(self) -> List[str]:
        """Get list of all supported drug names."""
        return list(self._drugs.keys())

    # ===== Utility Methods =====

    def is_gene_supported(self, gene: str) -> bool:
        """Check if a gene is supported."""
        return gene in self._genes

    def is_drug_supported(self, drug: str) -> bool:
        """Check if a drug is supported."""
        return drug.lower() in self._drugs

    def find_alleles_for_variant(self, gene: str, variant_key: str) -> List[str]:
        """
        Find which alleles contain a specific variant.
        variant_key format: "POS:REF:ALT"
        """
        variant_map = self.get_variant_to_allele_map(gene)
        return variant_map.get(variant_key, [])

    def get_allele_variants(self, gene: str, allele: str) -> List[str]:
        """Get all variants that define a specific allele."""
        allele_defs = self.get_allele_definitions(gene)
        return allele_defs.get(allele, [])

    @lru_cache(maxsize=128)
    def normalize_diplotype(self, diplotype: str) -> str:
        """
        Normalize diplotype notation (e.g., "*2/*1" -> "*1/*2").
        Ensures consistent ordering for lookup.
        """
        if '/' not in diplotype:
            return diplotype

        alleles = diplotype.split('/')
        # Sort alleles: *1 first, then numerically
        def sort_key(allele):
            if allele == '*1':
                return (0, 0)
            # Extract number from allele name
            import re
            match = re.search(r'\*(\d+)', allele)
            if match:
                return (1, int(match.group(1)))
            return (2, allele)

        sorted_alleles = sorted(alleles, key=sort_key)
        return '/'.join(sorted_alleles)

    def lookup_phenotype(self, gene: str, diplotype: str) -> Optional[str]:
        """
        Look up phenotype for a diplotype, trying both the original
        and normalized (reversed) forms.
        """
        phenotype_map = self.get_phenotype_map(gene)

        # Try direct lookup
        if diplotype in phenotype_map:
            return phenotype_map[diplotype]

        # Try normalized version
        normalized = self.normalize_diplotype(diplotype)
        if normalized in phenotype_map:
            return phenotype_map[normalized]

        return None

    # ===== NEW: Enriched Cache Accessors =====

    def get_rsid_map(self, gene: str) -> Dict[str, int]:
        """
        Get rsID → genomic position mapping for a gene.
        Returns: {"rs4244285": 94781859, ...}
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('rsid_map', {})

    def get_allele_function(self, gene: str, allele: str) -> Optional[str]:
        """
        Get the clinical functional status for an allele.
        Returns: e.g. "Normal function", "No function", "Decreased function"
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return None
        return gene_data.get('allele_functions', {}).get(allele)

    def get_all_allele_functions(self, gene: str) -> Dict[str, str]:
        """Get all allele → function status mappings for a gene."""
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('allele_functions', {})

    def get_frequencies(self, gene: str) -> Dict[str, Dict[str, float]]:
        """
        Get population allele frequencies for a gene.
        Returns: {"European": {"*1": 0.625, ...}, ...}
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return {}
        return gene_data.get('frequencies', {})

    # ===== Activity Scoring =====

    _FUNCTION_TO_SCORE = {
        "Normal function": 1.0,
        "Increased function": 1.5,
        "Decreased function": 0.5,
        "No function": 0.0,
        "Uncertain function": None,  # Cannot assign score
    }

    def get_activity_score(self, gene: str, allele: str) -> float:
        """
        Get activity score for an allele.
        Priority: direct cache score > fallback table > function-derived > default.
        """
        gene_data = self.get_gene_data(gene)

        # 1. Direct activity score from cache
        if gene_data:
            activity_scores = gene_data.get('activity_scores', {})
            if allele in activity_scores:
                return activity_scores[allele]

        # 2. Hardcoded fallback table (CPIC-sourced, covers missing cache data)
        fallback_table = _FALLBACK_ACTIVITY.get(gene, {})
        if allele in fallback_table:
            return fallback_table[allele]

        # 3. Derive from function status
        if gene_data:
            func_status = gene_data.get('allele_functions', {}).get(allele)
            if func_status:
                derived = self._FUNCTION_TO_SCORE.get(func_status)
                if derived is not None:
                    return derived

        # 4. Wildtype default
        if allele in ("*1", "Reference"):
            return 1.0

        return 1.0  # Conservative default

    def calculate_total_activity_score(self, gene: str, diplotype: str) -> float:
        """Calculate total activity score for a diplotype."""
        if '/' not in diplotype:
            return self.get_activity_score(gene, diplotype)

        alleles = diplotype.split('/')
        return sum(self.get_activity_score(gene, allele) for allele in alleles)


# Global singleton instance
_loader_instance: Optional[CPICDataLoader] = None


def get_cpic_loader() -> CPICDataLoader:
    """Get the global CPIC data loader instance."""
    global _loader_instance
    if _loader_instance is None:
        _loader_instance = CPICDataLoader()
    return _loader_instance


def reload_cpic_data():
    """Reload CPIC data (useful for testing or after ETL updates)."""
    global _loader_instance
    CPICDataLoader._initialized = False
    _loader_instance = CPICDataLoader()
    return _loader_instance
