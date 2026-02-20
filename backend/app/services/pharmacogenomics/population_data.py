"""
Population Frequency Data for Pharmacogenomic Alleles

Contains population-specific allele frequencies to inform phasing decisions
and provide priors for ambiguous diplotype calls.

Data sources:
- PharmGKB (https://www.pharmgkb.org)
- PharmVar (https://www.pharmvar.org)
- 1000 Genomes Project
- gnomAD

Note: These are approximate frequencies for common populations.
For production use, consider using live APIs or more detailed datasets.
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class AlleleFrequency:
    """Allele frequency data for a specific population."""
    allele: str
    frequency: float
    population: str


# Population codes
class Population:
    """Standard population codes."""
    GLOBAL = "global"  # Aggregated across populations
    EUR = "eur"  # European
    AFR = "afr"  # African
    EAS = "eas"  # East Asian
    SAS = "sas"  # South Asian
    AMR = "amr"  # Ad Mixed American


# ===== CYP2D6 Population Frequencies =====

CYP2D6_FREQUENCIES = {
    Population.EUR: {
        "*1": 0.35,   # Wild type
        "*2": 0.28,   # Normal function
        "*3": 0.01,   # Non-functional
        "*4": 0.20,   # Non-functional (most common variant allele in EUR)
        "*5": 0.03,   # Gene deletion
        "*6": 0.01,   # Non-functional
        "*9": 0.03,   # Decreased function
        "*10": 0.02,  # Decreased function
        "*17": 0.01,  # Decreased function
        "*41": 0.08,  # Decreased function
    },
    Population.EAS: {
        "*1": 0.35,
        "*2": 0.15,
        "*4": 0.01,   # Much less common in East Asian
        "*5": 0.06,
        "*10": 0.40,  # Very common in East Asian populations
        "*36": 0.02,
        "*41": 0.01,
    },
    Population.AFR: {
        "*1": 0.40,
        "*2": 0.20,
        "*4": 0.03,   # Less common in African
        "*5": 0.02,
        "*17": 0.20,  # Much more common in African populations
        "*29": 0.08,
        "*41": 0.02,
    }
}

# ===== CYP2C19 Population Frequencies =====

CYP2C19_FREQUENCIES = {
    Population.EUR: {
        "*1": 0.65,
        "*2": 0.15,   # Loss of function
        "*3": 0.001,  # Loss of function (rare)
        "*17": 0.20,  # Increased function
    },
    Population.EAS: {
        "*1": 0.35,
        "*2": 0.30,   # Much more common in East Asian
        "*3": 0.05,   # More common in East Asian
        "*17": 0.05,
    },
    Population.AFR: {
        "*1": 0.60,
        "*2": 0.18,
        "*3": 0.001,
        "*17": 0.20,
    }
}

# ===== CYP2C9 Population Frequencies =====

CYP2C9_FREQUENCIES = {
    Population.EUR: {
        "*1": 0.75,
        "*2": 0.13,   # Decreased function
        "*3": 0.08,   # Decreased function
    },
    Population.EAS: {
        "*1": 0.95,
        "*2": 0.001,  # Rare in East Asian
        "*3": 0.04,
    },
    Population.AFR: {
        "*1": 0.85,
        "*2": 0.01,   # Rare in African
        "*3": 0.01,
        "*5": 0.02,
        "*6": 0.02,
        "*8": 0.05,
        "*11": 0.03,
    }
}


class PopulationFrequencyData:
    """Access population frequency data for pharmacogenomic alleles."""

    def __init__(self):
        self.frequencies = {
            "CYP2D6": CYP2D6_FREQUENCIES,
            "CYP2C19": CYP2C19_FREQUENCIES,
            "CYP2C9": CYP2C9_FREQUENCIES,
        }

    def get_allele_frequency(
        self,
        gene: str,
        allele: str,
        population: str = Population.EUR
    ) -> float:
        """
        Get frequency of an allele in a specific population.

        Args:
            gene: Gene symbol (e.g., "CYP2D6")
            allele: Star allele name (e.g., "*4")
            population: Population code (default: EUR)

        Returns:
            Frequency (0.0 to 1.0), or 0.0 if unknown
        """
        gene_freqs = self.frequencies.get(gene, {})
        pop_freqs = gene_freqs.get(population, {})
        return pop_freqs.get(allele, 0.0)

    def get_diplotype_probability(
        self,
        gene: str,
        diplotype: str,
        population: str = Population.EUR,
        assume_hwe: bool = True
    ) -> float:
        """
        Calculate probability of a diplotype given population frequencies.

        Args:
            gene: Gene symbol
            diplotype: Diplotype (e.g., "*1/*4")
            population: Population code
            assume_hwe: Assume Hardy-Weinberg equilibrium

        Returns:
            Probability (0.0 to 1.0)
        """
        if '/' not in diplotype:
            return 0.0

        alleles = diplotype.split('/')
        freq1 = self.get_allele_frequency(gene, alleles[0], population)
        freq2 = self.get_allele_frequency(gene, alleles[1], population)

        if assume_hwe:
            if alleles[0] == alleles[1]:
                # Homozygous: p^2
                return freq1 ** 2
            else:
                # Heterozygous: 2pq
                return 2 * freq1 * freq2
        else:
            # Simple multiplication (less accurate)
            return freq1 * freq2

    def rank_diplotypes_by_probability(
        self,
        gene: str,
        diplotypes: List[str],
        population: str = Population.EUR
    ) -> List[tuple]:
        """
        Rank multiple possible diplotypes by their population probability.

        Args:
            gene: Gene symbol
            diplotypes: List of possible diplotypes
            population: Population code

        Returns:
            List of (diplotype, probability) tuples, sorted by probability (descending)
        """
        scored = []
        for diplotype in diplotypes:
            prob = self.get_diplotype_probability(gene, diplotype, population)
            scored.append((diplotype, prob))

        return sorted(scored, key=lambda x: x[1], reverse=True)

    def get_most_likely_phase(
        self,
        gene: str,
        allele1: str,
        allele2: str,
        population: str = Population.EUR
    ) -> tuple:
        """
        Determine most likely phasing for two alleles.

        Args:
            gene: Gene symbol
            allele1: First allele
            allele2: Second allele
            population: Population code

        Returns:
            (diplotype, probability, phase) where phase is "trans" or "cis"
        """
        # For compound heterozygotes, compare:
        # Trans: allele1/allele2
        # Cis: allele1/allele1 or allele2/allele2 might be more common

        trans_diplotype = f"{allele1}/{allele2}"
        trans_prob = self.get_diplotype_probability(gene, trans_diplotype, population)

        # Check if one homozygous form is more likely
        homo1 = f"{allele1}/{allele1}"
        homo2 = f"{allele2}/{allele2}"

        homo1_prob = self.get_diplotype_probability(gene, homo1, population)
        homo2_prob = self.get_diplotype_probability(gene, homo2, population)

        # Trans is typically assumed unless homozygous is much more likely
        if trans_prob >= max(homo1_prob, homo2_prob):
            return (trans_diplotype, trans_prob, "trans")
        elif homo1_prob > homo2_prob:
            return (homo1, homo1_prob, "cis")
        else:
            return (homo2, homo2_prob, "cis")

    def get_common_alleles(
        self,
        gene: str,
        population: str = Population.EUR,
        min_frequency: float = 0.01
    ) -> List[str]:
        """
        Get list of common alleles (frequency >= min_frequency) for a gene.

        Args:
            gene: Gene symbol
            population: Population code
            min_frequency: Minimum frequency threshold

        Returns:
            List of allele names
        """
        gene_freqs = self.frequencies.get(gene, {})
        pop_freqs = gene_freqs.get(population, {})

        common = [
            allele for allele, freq in pop_freqs.items()
            if freq >= min_frequency
        ]

        return sorted(common, key=lambda a: pop_freqs[a], reverse=True)


# Global instance
_pop_freq_instance: Optional[PopulationFrequencyData] = None


def get_population_frequencies() -> PopulationFrequencyData:
    """Get global population frequency data instance."""
    global _pop_freq_instance
    if _pop_freq_instance is None:
        _pop_freq_instance = PopulationFrequencyData()

# Alias for backward compatibility with RiskEngine
PopulationDataLoader = PopulationFrequencyData
