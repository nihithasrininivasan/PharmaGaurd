"""
Phenotype Mapper - Diplotype resolution and phenotype determination.
Implements star allele calling algorithm based on variant data.
"""

from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import re

from .models import GenotypeData, DiplotypeResult, VariantCall, IndeterminateReason
from .cpic_loader import get_cpic_loader
from .config import get_config, get_confidence_penalties, get_diplotype_config
from .population_data import get_population_frequencies, Population


class DiplotypeResolver:
    """Resolves diplotypes from variant calls using star allele calling logic."""

    def __init__(self, population: str = Population.EUR):
        self.loader = get_cpic_loader()
        self.config = get_config()
        self.penalties = get_confidence_penalties()
        self.diplotype_config = get_diplotype_config()
        self.pop_freq = get_population_frequencies()
        self.population = population  # Default to European population

    def resolve_diplotype(self, genotype_data: GenotypeData) -> DiplotypeResult:
        """
        Main entry point for diplotype resolution.
        Returns a DiplotypeResult with diplotype, phenotype, and confidence.
        """
        gene = genotype_data.gene_symbol

        # Check if gene is supported
        if not self.loader.is_gene_supported(gene):
            return DiplotypeResult(
                gene=gene,
                diplotype="Unknown",
                phenotype="Unknown",
                confidence=0.0,
                is_indeterminate=True,
                indeterminate_reason=IndeterminateReason.UNSUPPORTED_GENE,
                notes=f"Gene {gene} not supported"
            )

        # If no variants, assume wildtype
        if not genotype_data.variants:
            return self._resolve_wildtype(gene)

        # Get variant keys from genotype data
        variant_keys = genotype_data.get_variant_keys()

        # Identify candidate alleles
        candidate_alleles = self._identify_candidate_alleles(gene, genotype_data)

        # Check for phasing information
        has_phasing = any(v.phased for v in genotype_data.variants)

        # Resolve diplotype from candidates
        diplotype, confidence, notes, indet_reason = self._select_best_diplotype(
            gene, candidate_alleles, genotype_data, has_phasing
        )

        # Map diplotype to phenotype
        phenotype = self._map_phenotype(gene, diplotype)

        # Adjust confidence based on coverage
        final_confidence, coverage_issue = self._adjust_confidence_for_coverage(
            gene, confidence, genotype_data
        )

        # Determine final indeterminate status and reason
        is_indeterminate = diplotype == "Indeterminate" or final_confidence < 0.5

        # Use most specific indeterminate reason
        if coverage_issue and indet_reason == IndeterminateReason.NONE:
            indet_reason = IndeterminateReason.NO_COVERAGE
        elif is_indeterminate and indet_reason == IndeterminateReason.NONE:
            indet_reason = IndeterminateReason.LOW_QUALITY

        return DiplotypeResult(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            confidence=final_confidence,
            is_indeterminate=is_indeterminate,
            indeterminate_reason=indet_reason,
            notes=notes,
            phased=has_phasing
        )

    def _resolve_wildtype(self, gene: str) -> DiplotypeResult:
        """Handle case with no variants (wildtype)."""
        diplotype = "*1/*1"
        phenotype = self._map_phenotype(gene, diplotype)

        return DiplotypeResult(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            confidence=1.0,
            is_indeterminate=False,
            notes="No variants detected, assumed wildtype"
        )

    def _identify_candidate_alleles(
        self, gene: str, genotype_data: GenotypeData
    ) -> Dict[str, float]:
        """
        Identify candidate star alleles based on observed variants.
        Returns dict of {allele: match_score}
        """
        allele_definitions = self.loader.get_allele_definitions(gene)
        variant_to_allele = self.loader.get_variant_to_allele_map(gene)

        # Count matches for each allele
        allele_scores: Dict[str, float] = defaultdict(float)

        # Track which variants we've seen
        observed_variants = set(genotype_data.get_variant_keys())

        for variant_key in observed_variants:
            # Find alleles that contain this variant
            matching_alleles = variant_to_allele.get(variant_key, [])

            for allele in matching_alleles:
                # Get zygosity for this variant
                variant = next(
                    (v for v in genotype_data.variants if v.variant_key() == variant_key),
                    None
                )

                if variant:
                    if variant.zygosity == "HOM_ALT":
                        # Homozygous variant counts for both alleles
                        allele_scores[allele] += 2.0
                    elif variant.zygosity == "HET":
                        # Heterozygous variant counts for one allele
                        allele_scores[allele] += 1.0

        # Normalize scores by allele definition size (to favor complete matches)
        normalized_scores = {}
        for allele, score in allele_scores.items():
            allele_variants = set(self.loader.get_allele_variants(gene, allele))
            if allele_variants:
                # Check if all defining variants are present
                all_present = allele_variants.issubset(observed_variants)
                if all_present:
                    # Complete match - give full score
                    normalized_scores[allele] = score
                else:
                    # Partial match - penalize
                    completeness = len(allele_variants.intersection(observed_variants)) / len(allele_variants)
                    normalized_scores[allele] = score * completeness * 0.7

        return normalized_scores

    def _select_best_diplotype(
        self, gene: str, candidate_alleles: Dict[str, float], genotype_data: GenotypeData,
        has_phasing: bool = False
    ) -> Tuple[str, float, Optional[str], IndeterminateReason]:
        """
        Select the best diplotype from candidate alleles.
        Returns (diplotype, confidence, notes, indeterminate_reason)
        """
        if not candidate_alleles:
            # No matching alleles found
            if genotype_data.variants:
                return ("Indeterminate", 0.5, "Variants present but no matching alleles",
                        IndeterminateReason.NOVEL_VARIANTS)
            else:
                return ("*1/*1", 1.0, "No variants", IndeterminateReason.NONE)

        # Sort candidates by score
        sorted_candidates = sorted(
            candidate_alleles.items(), key=lambda x: x[1], reverse=True
        )

        # Analyze variant zygosity to determine diplotype
        het_variants = [v for v in genotype_data.variants if v.zygosity == "HET"]
        hom_variants = [v for v in genotype_data.variants if v.zygosity == "HOM_ALT"]

        # Case 1: All variants are homozygous for single allele
        if not het_variants and len(sorted_candidates) == 1:
            allele = sorted_candidates[0][0]
            return (f"{allele}/{allele}", 0.95, "Homozygous variant allele", IndeterminateReason.NONE)

        # Case 2: Single allele with high score (likely homozygous)
        if len(sorted_candidates) == 1 or sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
            allele = sorted_candidates[0][0]
            if sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
                return (f"{allele}/{allele}", 0.9, "Likely homozygous", IndeterminateReason.NONE)
            else:
                return (f"*1/{allele}", 0.85, "Heterozygous with wildtype", IndeterminateReason.NONE)

        # Case 3: Two strong candidates (compound heterozygote)
        if len(sorted_candidates) >= 2:
            allele1 = sorted_candidates[0][0]
            allele2 = sorted_candidates[1][0]

            # Check if scores suggest compound het
            score1 = sorted_candidates[0][1]
            score2 = sorted_candidates[1][1]

            min_score = self.diplotype_config.compound_het_min_score

            if score1 >= min_score and score2 >= min_score:
                # Both alleles have at least one copy
                diplotype = self.loader.normalize_diplotype(f"{allele1}/{allele2}")

                # Confidence depends on phasing
                if has_phasing:
                    confidence = 0.9  # High confidence with phasing
                    notes = "Compound heterozygote (phased)"
                    indet = IndeterminateReason.NONE
                else:
                    # Use population frequency priors to inform confidence
                    # Check if trans configuration is more likely than cis
                    most_likely_dip, pop_prob, phase = self.pop_freq.get_most_likely_phase(
                        gene, allele1, allele2, self.population
                    )

                    base_confidence = min(0.8, (score1 + score2) / 4.0)

                    if phase == "trans":
                        # Trans is consistent with our assumption - boost confidence slightly
                        confidence = base_confidence * 1.05 * self.penalties.unphased_heterozygote
                        notes = f"Compound heterozygote (trans assumed, pop. freq. consistent)"
                        indet = IndeterminateReason.NONE if confidence >= 0.6 else IndeterminateReason.AMBIGUOUS
                    else:
                        # Cis might be more likely - note this but keep trans assumption per CPIC
                        confidence = base_confidence * self.penalties.unphased_heterozygote
                        notes = f"Compound heterozygote (trans assumed, cis may be more common)"
                        indet = IndeterminateReason.AMBIGUOUS if confidence < 0.6 else IndeterminateReason.NONE

                return (diplotype, min(confidence, 0.9), notes, indet)
            else:
                # One dominant allele
                return (f"*1/{allele1}", 0.8, "Heterozygous with wildtype", IndeterminateReason.NONE)

        # Default: heterozygous with wildtype
        allele = sorted_candidates[0][0]
        return (f"*1/{allele}", 0.7, "Default heterozygous call", IndeterminateReason.PARTIAL_MATCH)

    def _map_phenotype(self, gene: str, diplotype: str) -> str:
        """Map diplotype to phenotype using CPIC phenotype map."""
        if diplotype in ["Unknown", "Indeterminate"]:
            return "Indeterminate"

        # Try direct lookup
        phenotype = self.loader.lookup_phenotype(gene, diplotype)
        if phenotype:
            return phenotype

        # Try activity score-based mapping
        try:
            activity_score = self.loader.calculate_total_activity_score(gene, diplotype)
            return self._activity_score_to_phenotype(gene, activity_score)
        except:
            return "Indeterminate"

    def _activity_score_to_phenotype(self, gene: str, total_score: float) -> str:
        """
        Map total activity score to phenotype.
        Based on CPIC activity score ranges.
        """
        # Standard CPIC phenotype mappings by activity score
        if gene in ["CYP2D6", "CYP2C19", "CYP2C9"]:
            if total_score == 0:
                return "PM"  # Poor Metabolizer
            elif 0 < total_score < 1.0:
                return "PM"
            elif 1.0 <= total_score < 1.5:
                return "IM"  # Intermediate Metabolizer
            elif 1.5 <= total_score < 2.5:
                return "NM"  # Normal Metabolizer
            elif total_score >= 2.5:
                return "UM"  # Ultra-rapid Metabolizer
            else:
                return "RM"  # Rapid Metabolizer

        # Default mapping
        if total_score < 1.0:
            return "PM"
        elif total_score < 1.5:
            return "IM"
        elif total_score < 2.5:
            return "NM"
        else:
            return "UM"

    def _adjust_confidence_for_coverage(
        self, gene: str, base_confidence: float, genotype_data: GenotypeData
    ) -> Tuple[float, bool]:
        """
        Adjust confidence score based on coverage information.
        Penalizes missing coverage at key positions.
        Returns (adjusted_confidence, has_coverage_issues)
        """
        confidence = base_confidence
        has_coverage_issues = False

        # Get key positions for this gene
        key_positions = set(self.loader.get_key_positions(gene))

        if not key_positions:
            return (confidence, False)  # No position data available

        # Check coverage
        covered_positions = set(genotype_data.covered_positions)

        if not covered_positions:
            # No coverage data provided
            confidence *= self.penalties.no_coverage_data
            return (max(0.0, min(1.0, confidence)), False)

        # Calculate missing positions
        missing_positions = key_positions - covered_positions

        if missing_positions:
            # Penalize for each missing key position
            penalty = self.penalties.missing_key_position ** len(missing_positions)
            confidence *= penalty
            has_coverage_issues = len(missing_positions) > 2  # Flag if >2 positions missing

        return (max(0.0, min(1.0, confidence)), has_coverage_issues)


class PhenotypeMapper:
    """High-level interface for phenotype mapping."""

    def __init__(self, population: str = Population.EUR):
        self.population = population
        self.resolver = DiplotypeResolver(population=population)

    def process_genotype(self, genotype_data: GenotypeData) -> DiplotypeResult:
        """Process genotype data and return diplotype result."""
        return self.resolver.resolve_diplotype(genotype_data)

    def process_multiple_genes(
        self, genotypes: List[GenotypeData]
    ) -> Dict[str, DiplotypeResult]:
        """Process multiple genes and return results keyed by gene."""
        results = {}
        for genotype in genotypes:
            result = self.process_genotype(genotype)
            results[genotype.gene_symbol] = result
        return results
