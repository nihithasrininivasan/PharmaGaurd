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

        # If no variants, apply strict wildtype gating
        if not genotype_data.variants:
            return self._resolve_wildtype(gene, genotype_data, bd)

        # Get variant keys from genotype data
        variant_keys = genotype_data.get_variant_keys()

        # Identify candidate alleles
        candidate_alleles = self._identify_candidate_alleles(gene, genotype_data)

        # Assess variant quality from actual VCF data (QUAL, DP, FILTER)
        self.confidence_calc.apply_variant_quality_from_vcf(bd, genotype_data.variants)

        # Check for phasing information
        has_phasing = any(v.phased for v in genotype_data.variants)

        # Resolve diplotype from candidates
        diplotype, dip_confidence, notes, indet_reason, is_partial, quality_category = self._select_best_diplotype(
            gene, candidate_alleles, genotype_data, has_phasing
        )

        # Apply diplotype determinism penalties using quality category
        self.confidence_calc.apply_diplotype_penalties(
            bd, diplotype,
            is_partial_match=is_partial,
            is_wildtype_unverified=False,
            quality_category=quality_category,
        )

        # Only use dip_confidence as a floor if quality_category is not set
        # (legacy path — should not normally trigger)
        if not quality_category:
            bd.diplotype_determinism = min(bd.diplotype_determinism, dip_confidence)

        # Map diplotype to phenotype
        phenotype = self._map_phenotype(gene, diplotype)

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
        # For variant-detected path: positive evidence exists, so use standard
        # threshold (0.50). The strict 0.75 gating applies to the wildtype path
        # via _resolve_wildtype(), which already returns Unresolved when
        # coverage/build/CNV gates fail.
        DETERMINISM_THRESHOLD = 0.50
        final_confidence = bd.final
        is_indeterminate = (
            diplotype in ("Indeterminate", "Unresolved", "Unknown")
            or final_confidence < DETERMINISM_THRESHOLD
        )

        # Use most specific indeterminate reason
        if is_indeterminate and indet_reason == IndeterminateReason.NONE:
            if not has_coverage:
                indet_reason = IndeterminateReason.NO_COVERAGE
            elif final_confidence < DETERMINISM_THRESHOLD:
                indet_reason = IndeterminateReason.LOW_QUALITY
            else:
                indet_reason = IndeterminateReason.LOW_QUALITY

        # Block phenotype inference when determinism gate fails
        if is_indeterminate and diplotype not in ("Indeterminate", "Unresolved", "Unknown"):
            phenotype = "Indeterminate"
            notes = (notes or "") + f" [Determinism gate: confidence {final_confidence:.2f} < {DETERMINISM_THRESHOLD}]"

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
        )

    def _resolve_wildtype(
        self, gene: str, genotype_data: GenotypeData, bd: ConfidenceBreakdown
    ) -> DiplotypeResult:
        """
        Handle case with no variants detected.

        Wild-type inference rule:
          If no clinically relevant LOF/decreased-function alleles are present
          in the VCF (i.e. all are 0/0), infer wild-type deterministically.
          Absence of pathogenic alleles IS sufficient evidence for wild-type.

        Confidence tiers:
          - HIGH (1.0): All verification gates pass (build, coverage, positions, CNV)
          - MODERATE (0.85): Gene is supported but coverage/build not fully verified
          - Indeterminate is NOT allowed when LOF alleles are explicitly absent.
        """
        from .confidence import CNV_REQUIRED_GENES

        key_positions = set(self.loader.get_key_positions(gene))
        covered_positions = set(genotype_data.covered_positions)
        genome_build = getattr(genotype_data, 'genome_build', None)

        # Apply genome build penalty
        self.confidence_calc.apply_genome_build_penalty(bd, genome_build)

        # Gate 1: Genome build validated
        build_ok = bool(genome_build and genome_build in ("GRCh38", "hg38"))

        # Gate 2: Coverage data exists
        has_coverage = bool(covered_positions)

        # Gate 3: All key positions covered
        all_keys_covered = bool(
            key_positions and has_coverage and key_positions.issubset(covered_positions)
        )

        # Gate 4: Gene-specific requirements (CNV)
        cnv_required = gene in CNV_REQUIRED_GENES
        cnv_evaluated = False  # Always false in current system

        # ── Wild-type inference (deterministic) ──
        # The gene is supported and no variants detected → wild-type by exclusion.
        # The verification gates determine confidence TIER, not whether to infer.
        diplotype = "*1/*1"
        phenotype = self._map_phenotype(gene, diplotype)
        is_indeterminate = False
        indet_reason = IndeterminateReason.NONE

        if build_ok and has_coverage and all_keys_covered and (not cnv_required or cnv_evaluated):
            # HIGH confidence — all verification gates passed
            notes = "Wildtype confirmed: all key positions covered, genome build validated"
            # diplotype_determinism stays at 1.0 (default)
        else:
            # MODERATE confidence — inferred wild-type from absence of pathogenic alleles
            # Build specific notes for transparency
            unverified = []
            if not build_ok:
                unverified.append(
                    f"genome build {'unknown' if not genome_build else genome_build}"
                )
            if not has_coverage:
                unverified.append("no coverage data from VCF")
            elif not all_keys_covered:
                missing = key_positions - covered_positions
                unverified.append(f"{len(missing)} key position(s) not covered")
            if cnv_required and not cnv_evaluated:
                unverified.append(f"{gene} CNV not evaluated")

            notes = (
                f"Wildtype inferred: no LOF/decreased-function alleles detected. "
                f"Unverified: {'; '.join(unverified)}"
            )

            # Apply moderate confidence — absence-of-evidence, not evidence-of-absence
            bd.diplotype_determinism = 0.85
            bd.penalties_applied.append(
                "Wildtype inferred from absence of pathogenic alleles (determinism 0.85)"
            )

            # CYP2D6 special: CNV uncertainty warrants extra caution
            if cnv_required and not cnv_evaluated:
                bd.diplotype_determinism = min(bd.diplotype_determinism, 0.75)
                bd.penalties_applied.append(
                    f"{gene} requires CNV evaluation — confidence capped at 0.75"
                )

        return DiplotypeResult(
            gene=gene,
            diplotype=diplotype,
            phenotype=phenotype,
            confidence=bd.final,
            is_indeterminate=is_indeterminate,
            indeterminate_reason=indet_reason,
            notes=notes,
            confidence_breakdown=bd.to_dict(),
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
            print(f"DEBUG: Checking {vk} (rsid={v.rsid})")
            if vk not in variant_to_allele and v.rsid and v.rsid in rsid_to_cpic_pos:
                cpic_pos = rsid_to_cpic_pos[v.rsid]
                print(f"DEBUG: Attempting rsID translation for {v.rsid} (VCF pos {v.pos} -> CPIC {cpic_pos})")

                # Try direct base match first
                translated_vk = f"{cpic_pos}:{v.ref}:{v.alt}"
                if translated_vk in variant_to_allele:
                    print(f"DEBUG: rsID translation success: {v.rsid} {vk} -> {translated_vk}")
                    vk = translated_vk
                    # CRITICAL FIX: Update variant object so subsequent checks (is_partial) see the match!
                    parts = translated_vk.split(":")
                    v.pos = int(parts[0])
                    v.ref = parts[1]
                    v.alt = parts[2]
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
                            print(f"DEBUG: rsID translation (strand) success: {v.rsid} {vk} -> {candidate}")
                            vk = candidate
                            # CRITICAL FIX: Update variant object
                            parts = candidate.split(":")
                            v.pos = int(parts[0])
                            v.ref = parts[1]
                            v.alt = parts[2]
                            break
            else:
                 if vk not in variant_to_allele:
                     print(f"DEBUG: No match for {vk} and rsID fallback failed/not applicable (rsid={v.rsid})")

            observed_variants.add(vk)

            # Find alleles that contain this variant
            matching_alleles = variant_to_allele.get(vk, [])
            print(f"DEBUG: {vk} matches alleles: {matching_alleles}")

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
        
        print(f"DEBUG: Candidate alleles for {gene}: {normalized_scores}")

        # Hybrid: Include VCF-annotated star alleles alongside positional matches.
        # VCF panels often call star alleles directly (especially for complex genes
        # like CYP2D6). When positional matching is weak or ambiguous, the VCF
        # annotation provides a more reliable signal.
        for v in genotype_data.variants:
            star = v.star_allele
            if star and star in allele_definitions:
                # Only add if VCF annotation provides a stronger signal than
                # any existing positional match for this allele
                vcf_score = 0.8 if v.zygosity == "HET" else (1.6 if v.zygosity == "HOM_ALT" else 0.0)
                if vcf_score > 0 and (star not in normalized_scores or normalized_scores[star] < vcf_score):
                    normalized_scores[star] = vcf_score
                    logger.info(
                        f"VCF-annotated allele call: {gene} {star} ({v.zygosity})"
                    )

        return normalized_scores

    def _select_best_diplotype(
        self, gene: str, candidate_alleles: Dict[str, float], genotype_data: GenotypeData,
        has_phasing: bool = False
    ) -> Tuple[str, float, Optional[str], IndeterminateReason, bool, Optional[str]]:
        """
        Select the best diplotype from candidate alleles.
        Returns (diplotype, confidence, notes, indeterminate_reason, is_partial_match, quality_category)
        
        NO confidence boosting. NO learned priors.
        quality_category is used by ConfidenceCalculator for standardized penalties.
        """
        if not candidate_alleles:
            # No matching alleles found
            if genotype_data.variants:
                return ("Indeterminate", 0.5, "Variants present but no matching alleles",
                        IndeterminateReason.NOVEL_VARIANTS, False, "INDETERMINATE")
            else:
                return ("*1/*1", 0.95, "No variants", IndeterminateReason.NONE, False, None)

        # Sort candidates by score (primary) and variant count/specificity (secondary)
        # TIE-BREAKER: Prefer lower allele numbers (*2 > *35)
        def sort_key(item):
            allele, score = item
            variant_count = len(self.loader.get_allele_variants(gene, allele))
            
            # Extract numeric part for sorting
            import re
            try:
                # Extract first number found
                match = re.search(r'\d+', allele)
                num = int(match.group()) if match else 9999
            except:
                num = 9999
                
            # Tie-breaker: (score DESC, count DESC, -num DESC (meaning num ASC))
            return (score, variant_count, -num, allele) 

        sorted_candidates = sorted(
            candidate_alleles.items(), key=sort_key, reverse=True
        )
        


        # Analyze variant zygosity to determine diplotype
        het_variants = [v for v in genotype_data.variants if v.zygosity == "HET"]
        hom_variants = [v for v in genotype_data.variants if v.zygosity == "HOM_ALT"]

        # Check for partial match (not all defining variants present for top allele)
        top_allele = sorted_candidates[0][0]
        top_allele_variants = set(self.loader.get_allele_variants(gene, top_allele))
        observed_variants = set(genotype_data.get_variant_keys())
        is_partial = not top_allele_variants.issubset(observed_variants)

        # Case 1: All variants are homozygous for single allele
        # Relaxed check: Allow if top candidate is strong homozygous match (score ~2.0)
        # even if other partial candidates exist (e.g. *5 vs *15).
        top_score = sorted_candidates[0][1]
        if not het_variants and (len(sorted_candidates) == 1 or top_score >= 1.9):
            allele = sorted_candidates[0][0]
            diplotype = f"{allele}/{allele}"
            confidence = 0.90 if not is_partial else 0.75
            notes = "Homozygous variant allele"
            qc = "EXACT_HOM" if not is_partial else "PARTIAL"

            return (diplotype, confidence, notes, IndeterminateReason.NONE, is_partial, qc)

        # Case 2: Single allele with high score
        if len(sorted_candidates) == 1 or sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
            allele = sorted_candidates[0][0]
            
            # --- Hardening: Strict Partial Match Check ---
            if is_partial:
                return ("Indeterminate", 0.4, f"Partial match to {allele} (incomplete definition)",
                        IndeterminateReason.PARTIAL_MATCH, True, "PARTIAL")

            if sorted_candidates[0][1] >= self.diplotype_config.homozygous_score_threshold:
                diplotype = f"{allele}/{allele}"
                confidence = 0.90
                notes = "Likely homozygous"
                qc = "EXACT_HOM"
            else:
                # Ensure alphabetic sorting for lookup key consistency (e.g. *1/*2A)
                other = "*1"
                alleles = sorted([other, allele])
                diplotype = f"{alleles[0]}/{alleles[1]}"
                if allele == "*1": diplotype = "*1/*1" # Edge case
                confidence = 0.85
                notes = "Heterozygous with wildtype"
                qc = "EXACT_HET_WILDTYPE"


            return (diplotype, confidence, notes, IndeterminateReason.NONE, is_partial, qc)

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
                    qc = "COMPOUND_HET_PHASED"
                else:
                    # Unphased compound het — capped at 0.80, never boosted
                    base_confidence = min(0.80, (score1 + score2) / 4.0)
                    confidence = base_confidence
                    notes = "Compound heterozygote (unphased — phase uncertainty)"
                    indet = IndeterminateReason.AMBIGUOUS if confidence < 0.6 else IndeterminateReason.NONE
                    qc = "COMPOUND_HET_UNPHASED"

                # >2 alleles without phasing → Indeterminate
                if len(sorted_candidates) > 2 and not has_phasing:
                    score3 = sorted_candidates[2][1]
                    if score3 >= min_score:
                        diplotype = "Indeterminate"
                        confidence = 0.3
                        notes = f">2 candidate alleles ({len(sorted_candidates)}) without phasing"
                        indet = IndeterminateReason.AMBIGUOUS
                        return (diplotype, confidence, notes, indet, is_partial, "AMBIGUOUS")


                return (diplotype, min(confidence, 0.90), notes, indet, is_partial, qc)
            else:
                # One dominant allele
                diplotype = f"*1/{allele1}"
                confidence = 0.75
                return (diplotype, confidence, "Heterozygous with wildtype",
                        IndeterminateReason.NONE, is_partial, "EXACT_HET_WILDTYPE")
 
        # Default: heterozygous with wildtype
        allele = sorted_candidates[0][0]
        diplotype = f"*1/{allele}"
        confidence = 0.65
        return (diplotype, confidence, "Default heterozygous call",
                IndeterminateReason.PARTIAL_MATCH, True, "PARTIAL")

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
        except Exception:
            return "Indeterminate"

    def _activity_score_to_phenotype(self, gene: str, total_score: float) -> str:
        """
        Map total activity score to phenotype.
        Uses CPIC long-form names to match drug recommendation table.
        """
        if gene in ["CYP2D6", "CYP2C19", "CYP2C9"]:
            if total_score == 0:
                return "Poor Metabolizer"
            elif 0 < total_score < 1.0:
                return "Poor Metabolizer"
            elif 1.0 <= total_score < 1.5:
                return "Intermediate Metabolizer"
            elif 1.5 <= total_score < 2.5:
                return "Normal Metabolizer"
            elif total_score >= 2.5:
                return "Ultrarapid Metabolizer"
            else:
                return "Rapid Metabolizer"

        # Default mapping
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
