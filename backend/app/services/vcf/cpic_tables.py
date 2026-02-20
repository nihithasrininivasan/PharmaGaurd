"""
cpic_tables.py
==============
CPIC-based star-allele activity scores and phenotype lookup for all 6 PGx genes.

Activity score model (used by CYP2D6, CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD):
  Each allele carries an activity value. The diplotype activity score = sum of
  both alleles. The score is then mapped to a phenotype per gene.

Sources: CPIC guidelines (https://cpicpgx.org/genes-drugs/)
"""
from __future__ import annotations
from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Activity values per star allele per gene
# ---------------------------------------------------------------------------
# Values: 0 = no function, 0.5 = decreased function, 1 = normal function,
#         1.5 = increased function (rare), 2 = CYP2D6 duplication normal
# "?" = unknown/indeterminate

_NA = None   # unknown

STAR_ACTIVITY: Dict[str, Dict[str, Optional[float]]] = {

    # ── CYP2D6 ──────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/
    "CYP2D6": {
        "*1":   1.0,   # Normal function (reference)
        "*1x2": 2.0,   # Gene duplication – increased function
        "*1xN": 2.0,   # Gene duplication N copies
        "*2":   1.0,   # Normal function
        "*2x2": 2.0,   # Duplication
        "*3":   0.0,   # No function (frameshift)
        "*4":   0.0,   # No function (splice defect)  ← rs3892097
        "*5":   0.0,   # No function (gene deletion)
        "*6":   0.0,   # No function (frameshift)     ← rs5030655
        "*7":   0.0,   # No function
        "*8":   0.0,   # No function
        "*9":   0.5,   # Decreased function
        "*10":  0.5,   # Decreased function
        "*14":  0.0,   # No function
        "*17":  0.5,   # Decreased function
        "*29":  0.5,   # Decreased function
        "*35":  1.0,   # Normal function
        "*36":  0.0,   # No function
        "*41":  0.5,   # Decreased function
        "*2A":  1.0,   # Alias for *2
    },

    # ── CYP2C19 ─────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/
    "CYP2C19": {
        "*1":   1.0,   # Normal function (reference)
        "*2":   0.0,   # No function (splice defect)   ← rs4244285
        "*3":   0.0,   # No function (premature stop)  ← rs4986893
        "*4":   0.0,   # No function
        "*5":   0.0,   # No function
        "*6":   0.0,   # No function
        "*7":   0.0,   # No function
        "*8":   0.0,   # No function
        "*9":   0.5,   # Decreased function
        "*10":  0.5,   # Decreased function
        "*17":  1.5,   # Increased function            ← rs12248560
        "*27":  0.5,   # Decreased function
        "*35":  0.5,   # Decreased function
    },

    # ── CYP2C9 ──────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/
    "CYP2C9": {
        "*1":   1.0,   # Normal function (reference)
        "*2":   0.5,   # Decreased function            ← rs1799853
        "*3":   0.0,   # No function                   ← rs1057910
        "*4":   0.0,   # No function
        "*5":   0.0,   # No function
        "*6":   0.0,   # No function
        "*8":   0.5,   # Decreased function
        "*11":  0.5,   # Decreased function
        "*12":  0.5,   # Decreased function
        "*13":  0.0,   # No function
        "*14":  0.5,   # Decreased function
    },

    # ── SLCO1B1 ─────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/
    # Activity model: normal / decreased / poor (not strictly numeric like CYP)
    # We encode: 1 = normal, 0.5 = decreased, 0 = poor
    "SLCO1B1": {
        "*1a":  1.0,   # Normal function (reference)
        "*1b":  1.0,   # Normal function
        "*1B":  1.0,   # Alias for *1b
        "*2":   0.5,   # Decreased
        "*3":   0.5,   # Decreased
        "*4":   0.5,   # Decreased
        "*5":   0.0,   # Poor function                 ← rs4149056
        "*6":   0.0,   # Poor function
        "*9":   0.5,   # Decreased
        "*10":  0.5,   # Decreased
        "*14":  0.5,   # Decreased
        "*15":  0.0,   # Poor function
        "*16":  0.0,   # Poor function
        "*17":  0.5,   # Decreased
        "*19":  0.0,   # Poor function
        "*20":  0.5,   # Decreased
        "*22":  0.5,   # Decreased
        "*24":  0.0,   # Poor function
        "*25":  0.5,   # Decreased
        "*28":  0.5,   # Decreased
        "*31":  0.5,   # Decreased
        "*37":  0.0,   # Poor function
        "*38":  0.5,   # Decreased
        "*45":  0.5,   # Decreased
    },

    # ── TPMT ────────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-azathioprine-and-tpmt/
    "TPMT": {
        "*1":   1.0,   # Normal (reference)
        "*2":   0.0,   # No function                   ← rs1800462
        "*3A":  0.0,   # No function (combines *3B+*3C)
        "*3B":  0.0,   # No function                   ← rs1800460
        "*3C":  0.0,   # No function                   ← rs1142345
        "*3D":  0.0,   # No function
        "*4":   0.0,   # No function
        "*5":   0.0,   # No function
        "*6":   0.0,   # No function
        "*7":   0.0,   # No function
        "*8":   0.5,   # Decreased function
        "*9":   0.5,   # Decreased function
        "*10":  0.5,   # Decreased function
        "*11":  0.5,   # Decreased function
        "*12":  0.5,   # Decreased function
        "*18":  0.0,   # No function
        "*19":  0.5,   # Decreased function
        "*20":  0.0,   # No function
        "*21":  0.5,   # Decreased function
        "*22":  0.5,   # Decreased function
        "*23":  0.5,   # Decreased function
        "*24":  0.5,   # Decreased function
        "*25":  0.5,   # Decreased function
        "*26":  0.0,   # No function
        "*28":  0.0,   # No function
        "*29":  0.0,   # No function
        "*38":  0.5,   # Decreased function
    },

    # ── DPYD ────────────────────────────────────────────────────────────────
    # CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/
    "DPYD": {
        "*1":   1.0,   # Normal (reference)
        "*2A":  0.0,   # No function (splice defect)   ← rs3918290
        "*3":   0.5,   # Decreased function
        "*4":   0.5,   # Decreased function
        "*5":   0.5,   # Decreased function
        "*6":   0.5,   # Decreased function
        "*7":   0.0,   # No function
        "*8":   0.5,   # Decreased function
        "*9A":  0.5,   # Decreased function
        "*9B":  0.5,   # Decreased function            ← rs55886062
        "*10":  0.5,   # Decreased function
        "*11":  0.0,   # No function
        "*12":  0.0,   # No function
        "*13":  0.0,   # No function                   ← rs67376798
    },
}

# Default: unknown alleles assume normal function (conservative for safety)
_DEFAULT_ACTIVITY = 1.0

# ---------------------------------------------------------------------------
# Activity score → Phenotype per gene
# ---------------------------------------------------------------------------
# Score ranges follow CPIC published activity score cutoffs

def _cyp2d6_phenotype(score: float) -> str:
    if score == 0:        return "Poor Metabolizer"
    if score <= 0.5:      return "Poor Metabolizer"
    if score <= 1.0:      return "Intermediate Metabolizer"
    if score <= 2.0:      return "Normal Metabolizer"
    return                        "Ultrarapid Metabolizer"

def _cyp2c19_phenotype(score: float) -> str:
    if score == 0:        return "Poor Metabolizer"
    if score <= 0.5:      return "Intermediate Metabolizer"
    if score <= 1.0:      return "Normal Metabolizer"
    if score <= 1.5:      return "Rapid Metabolizer"
    return                        "Ultrarapid Metabolizer"

def _cyp2c9_phenotype(score: float) -> str:
    if score == 0:        return "Poor Metabolizer"
    if score <= 0.5:      return "Poor Metabolizer"
    if score <= 1.0:      return "Intermediate Metabolizer"
    return                        "Normal Metabolizer"

def _slco1b1_phenotype(score: float) -> str:
    if score == 0:        return "Poor Function"
    if score <= 0.5:      return "Decreased Function"
    if score <= 1.5:      return "Normal Function"
    return                        "Increased Function"

def _tpmt_phenotype(score: float) -> str:
    if score == 0:        return "Poor Metabolizer"
    if score <= 0.5:      return "Intermediate Metabolizer"
    return                        "Normal Metabolizer"

def _dpyd_phenotype(score: float) -> str:
    if score == 0:        return "Poor Metabolizer"
    if score <= 0.5:      return "Intermediate Metabolizer"
    if score <= 1.5:      return "Normal Metabolizer"
    return                        "Normal Metabolizer"

_PHENOTYPE_FN = {
    "CYP2D6":  _cyp2d6_phenotype,
    "CYP2C19": _cyp2c19_phenotype,
    "CYP2C9":  _cyp2c9_phenotype,
    "SLCO1B1": _slco1b1_phenotype,
    "TPMT":    _tpmt_phenotype,
    "DPYD":    _dpyd_phenotype,
}


def _allele_activity(gene: str, star: Optional[str]) -> float:
    """Return activity value for a single allele. Unknown alleles → 1.0."""
    if not star:
        return _DEFAULT_ACTIVITY
    gene_table = STAR_ACTIVITY.get(gene, {})
    val = gene_table.get(star)
    if val is None:
        return _DEFAULT_ACTIVITY
    return val


def score_to_phenotype(gene: str, score: float) -> str:
    fn = _PHENOTYPE_FN.get(gene)
    if fn is None:
        return "Indeterminate"
    return fn(round(score, 2))


# ---------------------------------------------------------------------------
# Main entry point: variants → (diplotype, phenotype)
# ---------------------------------------------------------------------------

def infer_phenotype_from_variants(
    gene: str,
    variants: Sequence,          # Sequence[ExtractedVariant] — avoid circular import
) -> Tuple[str, str]:
    """
    Derive diplotype string and phenotype from a list of ExtractedVariant objects.

    Strategy:
      1. Collect unique star alleles observed.
      2. For each allele:
           - Hom-Alt  → that allele appears on BOTH chromosomes
           - Het      → one copy of the alt allele, one reference (*1)
           - Hom-Ref  → wildtype (*1) on both chromosomes
           - Unknown  → conservative: treat as *1 (normal)
      3. Build a two-element allele list covering both chromosomes.
      4. Sum activity scores → phenotype.

    Returns (diplotype_str, phenotype_str).
    """
    if not variants:
        # No variants detected → assume *1/*1 (wildtype)
        score = 2.0
        phenotype = score_to_phenotype(gene, score)
        return "*1/*1", phenotype

    # Collect (star, zygosity) pairs
    allele_observations: List[Tuple[str, str]] = []
    for v in variants:
        star = v.star or "*1"          # if no star tag, treat as wildtype
        zygosity = v.zygosity or "Unknown"
        allele_observations.append((star, zygosity))

    # Build the two chromosomal alleles
    # Simple model: use the first clearly non-wildtype allele found.
    # If multiple distinct stars seen, pick the two lowest-activity ones
    # (conservative / worst-case for patient safety).

    star_activity_pairs: List[Tuple[str, float]] = []
    for star, zyg in allele_observations:
        act = _allele_activity(gene, star)
        if zyg == "Hom-Alt":
            # Both chromosomes carry this allele
            star_activity_pairs.append((star, act))
            star_activity_pairs.append((star, act))
        elif zyg in ("Het", "Unknown"):
            # One chromosome carries the alt allele; other is *1 (reference)
            star_activity_pairs.append((star, act))
        # Hom-Ref → both chromosomes are *1, no alt allele → will default to *1/*1

    if not star_activity_pairs:
        # All variants were Hom-Ref → pure wildtype
        return "*1/*1", score_to_phenotype(gene, 2.0)

    # Sort by activity (lowest first = most impactful)
    star_activity_pairs.sort(key=lambda x: x[1])

    if len(star_activity_pairs) == 1:
        # Only one alt allele copy; other is assumed *1
        star1, act1 = star_activity_pairs[0]
        act2 = _allele_activity(gene, "*1")
        total = act1 + act2
        diplotype = f"{star1}/*1"
    else:
        star1, act1 = star_activity_pairs[0]
        star2, act2 = star_activity_pairs[1]
        total = act1 + act2
        diplotype = f"{star1}/{star2}"

    phenotype = score_to_phenotype(gene, total)
    return diplotype, phenotype
