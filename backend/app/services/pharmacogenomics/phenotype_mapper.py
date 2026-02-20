"""
Phenotype Mapper - Diplotype resolution and phenotype determination.
Implements star allele calling algorithm based on variant data.

All confidence scoring is deterministic and honest.
No confidence boosting. No learned priors.
"""

from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import logging
import re

logger = logging.getLogger(__name__)

from .models import GenotypeData, DiplotypeResult, VariantCall, IndeterminateReason
from .cpic_loader import get_cpic_loader
from .config import get_config, get_confidence_penalties, get_diplotype_config
from .confidence import ConfidenceBreakdown, ConfidenceCalculator

# Short code ↔ CPIC long name mapping (bidirectional)
PHENOTYPE_SHORT_TO_LONG = {
    "PM": "Poor Metabolizer",
    "IM": "Intermediate Metabolizer",
    "NM": "Normal Metabolizer",
    "RM": "Rapid Metabolizer",
    "UM": "Ultrarapid Metabolizer",
}
PHENOTYPE_LONG_TO_SHORT = {v: k for k, v in PHENOTYPE_SHORT_TO_LONG.items()}


class DiplotypeResolver:
    """Resolves diplotypes from variant calls using star allele calling logic."""

    def __init__(self):
        self.loader = get_cpic_loader()
        self.config = get_config()
        self.penalties = get_confidence_penalties()
        self.diplotype_config = get_diplotype_config()
        self.confidence_calc = ConfidenceCalculator()

    def resolve_diplotype(self, genotype_data: GenotypeData) -> DiplotypeResult:
        """
        Main entry point for diplotype resolution.
        Returns a DiplotypeResult with diplotype, phenotype, and confidence.
        """
        gene = genotype_data.gene_symbol
        bd = ConfidenceBreakdown()

        # Check if gene is supported
        if not self.loader.is_gene_supported(gene):
            self.confidence_calc.apply_diplotype_penalties(bd, "Unknown")
            self.confidence_calc.apply_cpic_penalties(bd, has_cpic_rule=False, phenotype_is_indeterminate=True)
            return DiplotypeResult(
                gene=gene,
                diplotype="Unknown",
                phenotype="Unknown",
                confidence=bd.final,
                is_indeterminate=True,
                indeterminate_reason=IndeterminateReason.UNSUPPORTED_GENE,
                notes=f"Gene {gene} not supported",
                confidence_breakdown=bd.to_dict(),
            )

        # Apply CNV penalty for genes that require CNV evaluation
        self.confidence_calc.apply_cnv_penalties(bd, gene, cnv_evaluated=False)

        # If no variants, assume wildtype — but with honest confidence
        if not genotype_data.variants:
            return self._resolve_wildtype(gene, genotype_data, bd)

        # Get variant keys from genotype data
        variant_keys = genotype_data.get_variant_keys()

        # Identify candidate alleles
        candidate_alleles = self._identify_candidate_alleles(gene, genotype_data)

        # Check for phasing information
        has_phasing = any(v.phased for v in genotype_data.variants)

        # Resolve diplotype from candidates
        diplotype, dip_confidence, notes, indet_reason, is_partial = self._select_best_diplotype(
            gene, candidate_alleles, genotype_data, has_phasing
        )

        # Apply diplotype determinism penalties
        self.confidence_calc.apply_diplotype_penalties(
            bd, diplotype,
            is_partial_match=is_partial,
            is_wildtype_unverified=False,
        )

        # Set diplotype_determinism from the resolution method's own confidence
        bd.diplotype_determinism = min(bd.diplotype_determinism, dip_confidence)

        # Map diplotype to phenotype and calculate activity score
        phenotype, total_as, allele_as = self._map_phenotype_with_scores(gene, diplotype)

        # Apply allele coverage penalties
        key_positions = self.loader.get_key_positions(gene)
        has_coverage = bool(genotype_data.covered_positions)
        self.confidence_calc.apply_allele_coverage_penalties(
            bd, gene,
            observed_positions=genotype_data.covered_positions,
            key_positions=key_positions,
            has_coverage_data=has_coverage,
        )

        # Apply phase penalties
        is_compound_het = (
            "/" in diplotype
            and diplotype not in ("Indeterminate", "Unknown")
            and diplotype.split("/")[0] != diplotype.split("/")[1]
            and diplotype.split("/")[0] != "*1"
        )
        self.confidence_calc.apply_phase_penalties(bd, is_compound_het, has_phasing)

        # Determine final indeterminate status
        final_confidence = bd.final
        is_indeterminate = diplotype == "Indeterminate" or final_confidence < 0.5

        # Use most specific indeterminate reason
        if is_indeterminate and indet_reason == IndeterminateReason.NONE:
            if not has_coverage:
                indet_reason = IndeterminateReason.NO_COVERAGE
            else:
                indet_reason = IndeterminateReason.LOW_QUALITY

        return DiplotypeResult(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            confidence=final_confidence,
            is_indeterminate=is_indeterminate,
            indeterminate_reason=indet_reason,
            notes=notes,
            phased=has_phasing,
            confidence_breakdown=bd.to_dict(),
            activity_score=total_as,
            allele_scores=allele_as,
        )

    def _resolve_wildtype(
        self, gene: str, genotype_data: GenotypeData, bd: ConfidenceBreakdown
    ) -> DiplotypeResult:
        """
        Handle case with no variants (wildtype).
        
        Honest confidence:
        - 0.70 if no coverage data (cannot confirm all positions are REF)
        - 0.85 if coverage data provided but incomplete
        - 0.95 if all key positions confirmed covered
        """
        diplotype = "*1/*1"
        phenotype, total_as, allele_as = self._map_phenotype_with_scores(gene, diplotype)

        key_positions = set(self.loader.get_key_positions(gene))
        covered_positions = set(genotype_data.covered_positions)

        if not covered_positions:
            # No coverage data — we have no proof positions are REF
            self.confidence_calc.apply_diplotype_penalties(
                bd, diplotype, is_wildtype_unverified=True
            )
            self.confidence_calc.apply_allele_coverage_penalties(
                bd, gene,
                observed_positions=[],
                key_positions=list(key_positions),
                has_coverage_data=False,
            )
            notes = "No variants detected; wildtype assumed but key positions not verified"
        elif key_positions and not key_positions.issubset(covered_positions):
            missing = key_positions - covered_positions
            self.confidence_calc.apply_allele_coverage_penalties(
                bd, gene,
                observed_positions=list(covered_positions),
                key_positions=list(key_positions),
                has_coverage_data=True,
            )
            notes = f"Wildtype assumed; {len(missing)} key position(s) not covered"
        else:
            notes = "No variants detected; all key positions covered — wildtype confirmed"

        return DiplotypeResult(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            confidence=bd.final,
            is_indeterminate=False,
            notes=notes,
            confidence_breakdown=bd.to_dict(),
            activity_score=total_as,
            allele_scores=allele_as,
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
        rsid_to_cpic_pos = self.loader.get_rsid_map(gene)

        # Count matches for each allele
        allele_scores: Dict[str, float] = defaultdict(float)

        # Track which variant keys (in CPIC coordinates) we've observed
        observed_variants: set = set()

        for v in genotype_data.variants:
            # Try direct position-based key first
            vk = v.variant_key()  # "POS:REF:ALT"

            # If no match in variant_to_allele, try rsID-based translation
            if vk not in variant_to_allele and v.rsid and v.rsid in rsid_to_cpic_pos:
                cpic_pos = rsid_to_cpic_pos[v.rsid]

                # Try direct base match first
                translated_vk = f"{cpic_pos}:{v.ref}:{v.alt}"
                if translated_vk in variant_to_allele:
                    logger.info(
                        f"rsID translation: {v.rsid} {vk} -> {translated_vk} (CPIC canonical)"
                    )
                    vk = translated_vk
                else:
                    # Try strand-flip combinations:
                    # VCF may report opposite strand or swapped REF/ALT
                    comp = str.maketrans("ACGTacgt", "TGCAtgca")
                    comp_ref = v.ref.translate(comp)
                    comp_alt = v.alt.translate(comp)
                    candidates = [
                        f"{cpic_pos}:{comp_ref}:{comp_alt}",   # complement
                        f"{cpic_pos}:{comp_alt}:{comp_ref}",   # reverse complement
                        f"{cpic_pos}:{v.alt}:{v.ref}",         # swapped REF/ALT
                    ]
                    for candidate in candidates:
                        if candidate in variant_to_allele:
                            logger.info(
                                f"rsID translation (strand): {v.rsid} {vk} -> {candidate}"
                            )
                            vk = candidate
                            break

            observed_variants.add(vk)

            # Find alleles that contain this variant
            matching_alleles = variant_to_allele.get(vk, [])

            for allele in matching_alleles:
                if v.zygosity == "HOM_ALT":
                    allele_scores[allele] += 2.0
                elif v.zygosity == "HET":
                    allele_scores[allele] += 1.0

        # Normalize scores by allele definition size
        normalized_scores = {}
        for allele, score in allele_scores.items():
            allele_variants = set(self.loader.get_allele_variants(gene, allele))
            if allele_variants:
                num_variants = len(allele_variants)

                # Check if all defining variants are present
                all_present = allele_variants.issubset(observed_variants)

                if all_present:
                    # Complete match - normalize by definition size
                    normalized_scores[allele] = score / num_variants
                else:
                    # Partial match - penalize
                    completeness = len(allele_variants.intersection(observed_variants)) / num_variants
                    base_norm_score = (score / num_variants)
                    normalized_scores[allele] = base_norm_score * completeness * 0.7

        # Apply VCF-annotated star allele boost
        # This gives precedence to star alleles explicitly called in the VCF (e.g. by a lab panel)
        for v in genotype_data.variants:
            star = v.star_allele
            if star and star in allele_definitions:
                if v.zygosity == "HOM_ALT":
                    normalized_scores[star] = max(normalized_scores.get(star, 0), 2.0 * 1.5)
                elif v.zygosity == "HET":
                    normalized_scores[star] = max(normalized_scores.get(star, 0), 1.0 * 1.5)
                logger.info(
                    f"VCF-annotated allele call boost applied: {gene} {star} ({v.zygosity})"
                )

        return normalized_scores

    def _select_best_diplotype(
        self, gene: str, candidate_alleles: Dict[str, float], genotype_data: GenotypeData,
        has_phasing: bool = False
    ) -> Tuple[str, float, Optional[str], IndeterminateReason, bool]:
        """
        Select the best diplotype from candidate alleles.
        Returns (diplotype, confidence, notes, indeterminate_reason, is_partial_match)
        
        NO confidence boosting. NO learned priors.
        """
        if not candidate_alleles:
            # No matching alleles found
            if genotype_data.variants:
                return ("Indeterminate", 0.5, "Variants present but no matching alleles",
                        IndeterminateReason.NOVEL_VARIANTS, False)
            else:
                return ("*1/*1", 0.95, "No variants", IndeterminateReason.NONE, False)

        # Prepare observed variants for sorting and completeness checks
        observed_variants = set(genotype_data.get_variant_keys())

        # Sort candidates by:
        # 1. Score (Primary) - Represents how well the variants match
        # 2. Completeness (Secondary) - Prefer alleles with more variants matching
        # 3. Specificity (Tertiary) - Prefer more specific alleles (more variants in definition)
        def sort_key(item):
            allele, score = item
            allele_variants = set(self.loader.get_allele_variants(gene, allele))
            num_vars = len(allele_variants)
            match_count = len(allele_variants.intersection(observed_variants))
            return (score, match_count, num_vars)

        sorted_candidates = sorted(
            candidate_alleles.items(), key=sort_key, reverse=True
        )

        # Analyze variant zygosity to determine diplotype
        het_variants = [v for v in genotype_data.variants if v.zygosity == "HET"]
        hom_variants = [v for v in genotype_data.variants if v.zygosity == "HOM_ALT"]

        # Check for partial match (not all defining variants present for top allele)
        top_allele = sorted_candidates[0][0]
        top_allele_variants = set(self.loader.get_allele_variants(gene, top_allele))
        is_partial = not top_allele_variants.issubset(observed_variants)

        # Case 1: All variants are homozygous for single allele
        if not het_variants and len(sorted_candidates) == 1:
            allele = sorted_candidates[0][0]
            diplotype = f"{allele}/{allele}"
            confidence = 0.90 if not is_partial else 0.75
            notes = "Homozygous variant allele"
            return (diplotype, confidence, notes, IndeterminateReason.NONE, is_partial)

        # Case 2: Single allele with high score (likely homozygous)
        if len(sorted_candidates) == 1 or sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
            allele = sorted_candidates[0][0]
            if sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
                diplotype = f"{allele}/{allele}"
                confidence = 0.85 if not is_partial else 0.70
                notes = f"Homozygous {allele} detected"
            else:
                diplotype = f"*1/{allele}"
                confidence = 0.80 if not is_partial else 0.65
                notes = f"Heterozygous {allele} with wildtype *1"

            # Add Structural Variant warning for CYP2D6
            if gene == "CYP2D6":
                notes += ". Note: Structural variants (deletion/duplication) not assessed; AS assumes CN=2."

            return (diplotype, confidence, notes, IndeterminateReason.NONE, is_partial)

        # Case 3: Two strong candidates (compound heterozygote)
        if len(sorted_candidates) >= 2:
            allele1 = sorted_candidates[0][0]
            allele2 = sorted_candidates[1][0]

            score1 = sorted_candidates[0][1]
            score2 = sorted_candidates[1][1]

            min_score = self.diplotype_config.compound_het_min_score

            if score1 >= min_score and score2 >= min_score:
                # Both alleles have at least one copy
                diplotype = self.loader.normalize_diplotype(f"{allele1}/{allele2}")

                if has_phasing:
                    confidence = 0.85  # High confidence with phasing
                    notes = "Compound heterozygote (phased)"
                    indet = IndeterminateReason.NONE
                else:
                    # Unphased compound het — capped at 0.80, never boosted
                    base_confidence = min(0.80, (score1 + score2) / 4.0)
                    confidence = base_confidence
                    notes = "Compound heterozygote (unphased — phase uncertainty)"
                    indet = IndeterminateReason.AMBIGUOUS if confidence < 0.6 else IndeterminateReason.NONE

                # >2 alleles without phasing → Indeterminate
                if len(sorted_candidates) > 2 and not has_phasing:
                    score3 = sorted_candidates[2][1]
                    if score3 >= min_score:
                        diplotype = "Indeterminate"
                        confidence = 0.3
                        notes = f">2 candidate alleles ({len(sorted_candidates)}) without phasing"
                        indet = IndeterminateReason.AMBIGUOUS
                        return (diplotype, confidence, notes, indet, is_partial)

                return (diplotype, min(confidence, 0.90), notes, indet, is_partial)
            else:
                # One dominant allele
                diplotype = f"*1/{allele1}"
                confidence = 0.75
                return (diplotype, confidence, "Heterozygous with wildtype",
                        IndeterminateReason.NONE, is_partial)
 
        # Default: heterozygous with wildtype
        allele = sorted_candidates[0][0]
        diplotype = f"*1/{allele}"
        confidence = 0.65
        return (diplotype, confidence, "Default heterozygous call",
                IndeterminateReason.PARTIAL_MATCH, True)

    def _map_phenotype_with_scores(self, gene: str, diplotype: str) -> Tuple[str, Optional[float], Optional[Dict[str, float]]]:
        """Map diplotype to phenotype and return scores."""
        if diplotype in ["Unknown", "Indeterminate"]:
            return "Indeterminate", None, None

        total_as = None
        allele_as = None

        try:
            # 1. Calculate scores
            if '/' in diplotype:
                alleles = diplotype.split('/')
                a1_score = self.loader.get_activity_score(gene, alleles[0])
                a2_score = self.loader.get_activity_score(gene, alleles[1])
                total_as = a1_score + a2_score
                allele_as = {alleles[0]: a1_score, alleles[1]: a2_score}
            else:
                total_as = self.loader.get_activity_score(gene, diplotype)
                allele_as = {diplotype: total_as}
        except Exception:
            pass

        # 2. Lookup phenotype
        phenotype = self.loader.lookup_phenotype(gene, diplotype)
        if not phenotype and total_as is not None:
            phenotype = self._activity_score_to_phenotype(gene, total_as)

        return phenotype or "Indeterminate", total_as, allele_as

    def _map_phenotype(self, gene: str, diplotype: str) -> str:
        """Deprecated: Use _map_phenotype_with_scores."""
        p, _, _ = self._map_phenotype_with_scores(gene, diplotype)
        return p

    def _activity_score_to_phenotype(self, gene: str, total_score: float) -> str:
        """
        Map total activity score to phenotype.
        Uses gene-specific CPIC cutoffs and long-form names to match
        drug recommendation tables.
        """
        # CYP2D6: CPIC codeine guideline cutoffs
        if gene == "CYP2D6":
            if total_score <= 0:
                return "Poor Metabolizer"
            elif total_score < 1.0:
                return "Intermediate Metabolizer"
            elif total_score <= 2.25:
                return "Normal Metabolizer"
            else:
                return "Ultrarapid Metabolizer"

        # CYP2C19: CPIC clopidogrel guideline cutoffs
        if gene == "CYP2C19":
            if total_score <= 0:
                return "Poor Metabolizer"
            elif total_score < 2.0:
                return "Intermediate Metabolizer"
            elif total_score == 2.0:
                return "Normal Metabolizer"
            elif total_score <= 2.5:
                return "Rapid Metabolizer"
            else:
                return "Ultrarapid Metabolizer"

        # CYP2C9: CPIC warfarin guideline cutoffs
        if gene == "CYP2C9":
            if total_score <= 0.5:
                return "Poor Metabolizer"
            elif total_score < 2.0:
                return "Intermediate Metabolizer"
            else:
                return "Normal Metabolizer"

        # SLCO1B1: uses Function-based phenotypes (not Metabolizer)
        if gene == "SLCO1B1":
            if total_score <= 0:
                return "Poor Function"
            elif total_score <= 1.0:
                return "Decreased Function"
            else:
                return "Normal Function"

        # DPYD: CPIC fluoropyrimidine guideline cutoffs
        if gene == "DPYD":
            if total_score <= 0:
                return "Poor Metabolizer"
            elif total_score <= 1.0:
                return "Intermediate Metabolizer"
            else:
                return "Normal Metabolizer"

        # TPMT: CPIC thiopurine guideline cutoffs
        if gene == "TPMT":
            if total_score <= 0:
                return "Poor Metabolizer"
            elif total_score < 2.0:
                return "Intermediate Metabolizer"
            else:
                return "Normal Metabolizer"

        # Default mapping for unknown genes
        if total_score < 1.0:
            return "Poor Metabolizer"
        elif total_score < 1.5:
            return "Intermediate Metabolizer"
        elif total_score < 2.5:
            return "Normal Metabolizer"
        else:
            return "Ultrarapid Metabolizer"


class PhenotypeMapper:
    """High-level interface for phenotype mapping."""

    def __init__(self):
        self.resolver = DiplotypeResolver()

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
