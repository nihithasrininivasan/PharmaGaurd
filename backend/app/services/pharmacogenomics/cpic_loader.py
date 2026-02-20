"""
CPIC Data Loader - Runtime singleton for loading and accessing CPIC data.
Provides optimized lookup structures for diplotype resolution and risk assessment.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Set
from functools import lru_cache


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

    def get_drug_recommendation_for_phenotype(
        self, drug: str, phenotype: str, activity_score: Optional[float] = None
    ) -> Optional[Dict]:
        """
        Get specific recommendation for a drug-phenotype combination.

        Tries phenotype string match first (most drugs).
        Falls back to activity score lookup for drugs like codeine that use
        numeric score keys (e.g., '0.0', '≥3.0') instead of phenotype strings.
        """
        recommendations = self.get_drug_recommendations(drug)
        if not recommendations:
            return None

        # 1. Direct phenotype string lookup
        rec = recommendations.get(phenotype)
        if rec:
            return rec

        # 2. Activity score-based lookup (for codeine-style recommendations)
        if activity_score is not None:
            # Try exact score key first
            exact_key = str(float(activity_score))
            rec = recommendations.get(exact_key)
            if rec:
                return rec

            # Try threshold keys (≥N) — find the tightest applicable threshold
            threshold_keys = []
            for key in recommendations:
                if key.startswith('≥'):
                    try:
                        threshold = float(key[1:])
                        if activity_score >= threshold:
                            threshold_keys.append((threshold, key))
                    except ValueError:
                        pass

            if threshold_keys:
                # Use the highest applicable threshold (most specific match)
                threshold_keys.sort(key=lambda x: x[0], reverse=True)
                rec = recommendations.get(threshold_keys[0][1])
                if rec:
                    return rec

            # 3. CRITICAL SAFEGUARD: Do NOT guess closest numeric key.
            # If no exact match and no threshold match, return None.
            # This prevents 2.0 (Normal) falling back to 1.0 (Avoid) if 2.0 key is missing.
            pass

        return None

    def get_supported_drugs(self) -> List[str]:
        """Get list of all supported drug names."""
        return list(self._drugs.keys())

    # ===== Utility Methods =====

    def is_gene_supported(self, gene: str) -> bool:
        """Check if a gene is supported."""
        return gene in self._genes

    def is_drug_supported(self, drug: str) -> bool:
        """Check if a drug is supported.

        A drug is considered supported if:
        1. It exists in the local CPIC alert cache (cpic_cache.json), OR
        2. PharmGKB has Level 1A or 1B evidence for this drug across any gene.

        Case-insensitive. Aggregates across all genes.
        """
        drug_key = drug.strip().lower()
        # Check 1: CPIC alert files (local cache)
        if drug_key in self._drugs:
            return True
        # Check 2: Delegate to PharmGKB for Level 1A/1B evidence
        try:
            from .pharmgkb_loader import get_pharmgkb_loader
            pharmgkb = get_pharmgkb_loader()
            return pharmgkb.is_drug_supported(drug_key)
        except Exception:
            return False

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
        Priority: direct score > function-derived score > default.
        """
        gene_data = self.get_gene_data(gene)
        if not gene_data:
            return 1.0

        # 1. Direct activity score from cache
        activity_scores = gene_data.get('activity_scores', {})
        if allele in activity_scores:
            return activity_scores[allele]

        # 2. Derive from function status
        func_status = gene_data.get('allele_functions', {}).get(allele)
        if func_status:
            derived = self._FUNCTION_TO_SCORE.get(func_status)
            if derived is not None:
                return derived

        # 3. Wildtype default
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
