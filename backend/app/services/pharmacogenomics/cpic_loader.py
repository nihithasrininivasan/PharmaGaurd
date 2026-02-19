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

        if not cache_file.exists():
            # Check if auto-ETL is enabled
            if config.auto_run_etl:
                cpic_data_dir = backend_dir / config.cpic_data_dir
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

    def get_drug_recommendation_for_phenotype(self, drug: str, phenotype: str) -> Optional[Dict]:
        """Get specific recommendation for a drug-phenotype combination."""
        recommendations = self.get_drug_recommendations(drug)
        return recommendations.get(phenotype)

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

    def get_activity_score(self, gene: str, allele: str) -> float:
        """
        Get activity score for an allele (for activity score-based phenotype determination).
        Uses configuration-defined scores, with fallback to defaults.
        """
        from .config import get_activity_scores

        activity_config = get_activity_scores()

        # Check gene-specific scores first
        gene_scores = activity_config.gene_specific_scores.get(gene, {})
        if allele in gene_scores:
            return gene_scores[allele]

        # Default: wildtype has activity 1.0, unknown alleles default to 1.0
        if allele == "*1":
            return 1.0

        return 1.0  # Conservative default for unknown alleles

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
