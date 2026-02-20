"""
PharmGKB Data Loader — Deterministic lookup service for pharmacogenomic datasets.

Parses 4 PharmGKB TSV datasets at startup and provides lookup APIs for:
- Gene-drug pair confirmation (relationships.tsv)
- Variant-level annotation (clinicalVariants.tsv)
- Evidence level scoring (clinicalVariants.tsv)
- Clinical annotation linking (relationships.tsv)
- Drug metadata (drugs.tsv)

STRICT CONSTRAINTS:
- No hallucinated data — only returns what exists in the datasets
- No fabricated CPIC logic — data-driven lookups only
- Deterministic and auditable
"""

import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

# PharmGKB TSVs have very large fields (SMILES/InChI in drugs.tsv)
csv.field_size_limit(sys.maxsize)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Evidence level → confidence weight mapping
# ---------------------------------------------------------------------------
EVIDENCE_LEVEL_WEIGHTS: Dict[str, float] = {
    "1A": 1.00,
    "1B": 1.00,
    "2A": 0.85,
    "2B": 0.80,
    "3":  0.65,
    "4":  0.50,
}

# Levels that allow fully automated clinical recommendations
AUTOMATED_RECOMMENDATION_LEVELS = frozenset({"1A", "1B", "2A", "2B"})

# ---------------------------------------------------------------------------
# Data directory
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data"


class PharmGKBLoader:
    """Singleton-style loader for PharmGKB TSV datasets."""

    _instance: Optional["PharmGKBLoader"] = None

    def __init__(self, data_dir: Optional[Path] = None):
        self._data_dir = data_dir or DATA_DIR
        self._relationships: List[Dict[str, str]] = []
        self._clinical_variants: List[Dict[str, str]] = []
        self._drugs: Dict[str, Dict[str, str]] = {}  # keyed by lowercase drug name
        self._genes: Dict[str, Dict[str, str]] = {}   # keyed by gene symbol

        # Pre-built indexes for fast lookup
        self._gene_drug_index: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
        self._variant_gene_index: Dict[str, List[Dict]] = defaultdict(list)
        self._gene_drug_evidence_index: Dict[Tuple[str, str], str] = {}

        self._load_all()

    @classmethod
    def get_instance(cls, data_dir: Optional[Path] = None) -> "PharmGKBLoader":
        if cls._instance is None:
            cls._instance = cls(data_dir)
        return cls._instance

    # ===== Loading =====

    def _load_all(self) -> None:
        self._load_relationships()
        self._load_clinical_variants()
        self._load_drugs()
        self._load_genes()
        logger.info(
            "PharmGKB Loader initialized: %d relationships, %d clinical variants, "
            "%d drugs, %d genes",
            len(self._relationships),
            len(self._clinical_variants),
            len(self._drugs),
            len(self._genes),
        )

    def _load_relationships(self) -> None:
        """Load relationships.tsv and build gene→drug index."""
        path = self._data_dir / "relationships.tsv"
        if not path.exists():
            logger.warning("relationships.tsv not found at %s", path)
            return

        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                self._relationships.append(row)

                # Index Gene→Chemical relationships
                e1_type = row.get("Entity1_type", "")
                e2_type = row.get("Entity2_type", "")
                if e1_type == "Gene" and e2_type == "Chemical":
                    gene = row.get("Entity1_name", "").strip()
                    drug = row.get("Entity2_name", "").strip().lower()
                    if gene and drug:
                        self._gene_drug_index[(gene, drug)].append(row)
                # Also handle reverse: Chemical→Gene
                elif e1_type == "Chemical" and e2_type == "Gene":
                    gene = row.get("Entity2_name", "").strip()
                    drug = row.get("Entity1_name", "").strip().lower()
                    if gene and drug:
                        self._gene_drug_index[(gene, drug)].append(row)

    def _load_clinical_variants(self) -> None:
        """Load clinicalVariants.tsv and build variant→gene index."""
        path = self._data_dir / "clinicalVariants.tsv"
        if not path.exists():
            logger.warning("clinicalVariants.tsv not found at %s", path)
            return

        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                self._clinical_variants.append(row)
                gene = row.get("gene", "").strip()
                variant = row.get("variant", "").strip()
                if gene and variant:
                    self._variant_gene_index[gene].append(row)

        # Build gene-drug evidence level index (top evidence per pair)
        for row in self._clinical_variants:
            gene = row.get("gene", "").strip()
            level = row.get("level of evidence", "").strip()
            chemicals = row.get("chemicals", "").strip()
            if not gene or not chemicals or not level:
                continue
            for chem in chemicals.split(","):
                chem = chem.strip().lower()
                if not chem:
                    continue
                key = (gene, chem)
                existing = self._gene_drug_evidence_index.get(key)
                if existing is None or _level_rank(level) < _level_rank(existing):
                    self._gene_drug_evidence_index[key] = level

    def _load_drugs(self) -> None:
        """Load drugs.tsv."""
        path = self._data_dir / "drugs.tsv"
        if not path.exists():
            logger.warning("drugs.tsv not found at %s", path)
            return

        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = row.get("Name", "").strip()
                if name:
                    self._drugs[name.lower()] = row
                # Also index generic names
                for gn in (row.get("Generic Names") or "").split(","):
                    gn = gn.strip().strip('"').lower()
                    if gn and gn not in self._drugs:
                        self._drugs[gn] = row

    def _load_genes(self) -> None:
        """Load genes.tsv."""
        path = self._data_dir / "genes.tsv"
        if not path.exists():
            logger.warning("genes.tsv not found at %s", path)
            return

        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                symbol = row.get("Symbol", "").strip()
                if symbol:
                    self._genes[symbol] = row

    # ===== Public APIs =====

    def confirm_gene_drug_pair(
        self, gene: str, drug: str
    ) -> Dict:
        """
        Deterministically confirm whether a gene-drug pair exists in the
        PharmGKB relationships dataset.

        Returns dict with:
          - confirmed: bool
          - evidence_types: list of evidence types found
          - association: deterministic classification based on evidence quality
          - source_pmids: list of PMIDs
          - evidence_level: highest evidence level found (1A, 1B, etc.)
        """
        drug_key = drug.strip().lower()
        gene_key = gene.strip()
        rows = self._gene_drug_index.get((gene_key, drug_key), [])

        if not rows:
            return {
                "gene": gene_key,
                "drug": drug_key,
                "confirmed": False,
                "evidence_types": [],
                "association": "not found",
                "source_pmids": [],
                "evidence_level": "none",
            }

        evidence_types = set()
        associations = set()
        pmids = set()

        for row in rows:
            ev = row.get("Evidence", "").strip()
            if ev:
                for e in ev.split(","):
                    evidence_types.add(e.strip())
            assoc = row.get("Association", "").strip()
            if assoc:
                associations.add(assoc)
            pids = row.get("PMIDs", "").strip()
            if pids:
                for p in pids.split(";"):
                    p = p.strip()
                    if p:
                        pmids.add(p)

        # Get evidence level for this gene-drug pair
        evidence_info = self.get_evidence_level(gene_key, drug_key)
        evidence_level = evidence_info.get("level", "none")

        # Classify association using deterministic rule set
        association = _classify_association(
            confirmed=True,
            evidence_level=evidence_level,
            evidence_types=list(evidence_types),
            raw_associations=associations,
        )

        return {
            "gene": gene_key,
            "drug": drug_key,
            "confirmed": True,
            "evidence_types": sorted(evidence_types),
            "association": association,
            "source_pmids": sorted(pmids)[:20],  # Cap for readability
            "evidence_level": evidence_level,
        }

    def get_variant_annotations(
        self, gene: str, variant_ids: List[str]
    ) -> List[Dict]:
        """
        Pull functional annotations for specific variants in a gene.

        Args:
            gene: Gene symbol (e.g., "CYP2D6")
            variant_ids: List of rsIDs or star allele identifiers

        Returns list of annotation dicts with:
          - variant_id, gene, annotation_type, evidence_level,
            associated_chemicals, associated_phenotypes
        """
        if not variant_ids:
            return []

        gene_key = gene.strip()
        gene_rows = self._variant_gene_index.get(gene_key, [])
        if not gene_rows:
            return []

        # Normalize search IDs for matching
        search_ids = set()
        for vid in variant_ids:
            search_ids.add(vid.strip())
            # Also try without * prefix for star alleles
            if vid.startswith("*"):
                search_ids.add(f"{gene_key}{vid}")

        results = []
        for row in gene_rows:
            variant_field = row.get("variant", "").strip()
            # Match: exact rsID or star allele in the comma-separated variant field
            var_parts = [v.strip() for v in variant_field.split(",")]
            matched = False
            for part in var_parts:
                if part in search_ids:
                    matched = True
                    break

            if not matched:
                continue

            chemicals = [c.strip() for c in row.get("chemicals", "").split(",") if c.strip()]
            phenotypes_str = row.get("phenotypes", "").strip()
            phenotypes = [p.strip() for p in phenotypes_str.split(",") if p.strip()] if phenotypes_str else []

            results.append({
                "variant_id": variant_field,
                "gene": gene_key,
                "annotation_type": row.get("type", "").strip(),
                "evidence_level": row.get("level of evidence", "").strip(),
                "associated_chemicals": chemicals,
                "associated_phenotypes": phenotypes,
            })

        return results

    def get_evidence_level(
        self, gene: str, drug: str
    ) -> Dict:
        """
        Get the top (strongest) evidence level for a gene-drug pair.

        Returns dict with:
          - level: str (e.g., "1A")
          - confidence_weight: float
          - allows_automated_recommendation: bool
        """
        gene_key = gene.strip()
        drug_key = drug.strip().lower()
        level = self._gene_drug_evidence_index.get((gene_key, drug_key))

        if not level:
            return {
                "level": "none",
                "confidence_weight": 0.50,
                "allows_automated_recommendation": False,
            }

        weight = EVIDENCE_LEVEL_WEIGHTS.get(level, 0.50)
        allows_auto = level in AUTOMATED_RECOMMENDATION_LEVELS

        return {
            "level": level,
            "confidence_weight": weight,
            "allows_automated_recommendation": allows_auto,
        }

    def get_clinical_annotations(
        self, gene: str, drug: str
    ) -> List[Dict]:
        """
        Get clinical annotation links for a gene-drug pair from relationships.tsv.

        Returns list of dicts with:
          - annotation_id, gene, drug, evidence_type, association, pmids
        Deduplicated by (annotation_id, evidence_type) to prevent duplicates.
        """
        drug_key = drug.strip().lower()
        gene_key = gene.strip()
        rows = self._gene_drug_index.get((gene_key, drug_key), [])

        results = []
        seen_keys = set()  # Deduplication by (annotation_id, evidence_type)
        for row in rows:
            ev = row.get("Evidence", "").strip()
            if not ev:
                continue

            # Use the entity ID as annotation ID
            e1_id = row.get("Entity1_id", "").strip()
            e2_id = row.get("Entity2_id", "").strip()
            # Use whichever is the haplotype/variant entity
            e1_type = row.get("Entity1_type", "")
            if e1_type == "Gene":
                ann_id = e1_id
            else:
                ann_id = e2_id

            # Deduplicate by (annotation_id, evidence_type)
            dedup_key = (ann_id, ev)
            if dedup_key in seen_keys:
                continue
            seen_keys.add(dedup_key)

            pids = row.get("PMIDs", "").strip()
            pmid_list = [p.strip() for p in pids.split(";") if p.strip()] if pids else []

            results.append({
                "annotation_id": ann_id,
                "gene": gene_key,
                "drug": drug_key,
                "evidence_type": ev,
                "association": row.get("Association", "").strip(),
                "pmids": pmid_list[:10],  # Cap for readability
            })

        return results

    def get_drug_info(
        self, drug: str
    ) -> Optional[Dict]:
        """
        Get drug metadata from drugs.tsv.

        Returns dict with:
          - name, cpic_pairs_level, clinical_annotation_level,
            dosing_guideline_sources, has_dosing_info
        """
        drug_key = drug.strip().lower()
        row = self._drugs.get(drug_key)
        if not row:
            return None

        return {
            "name": row.get("Name", "").strip(),
            "cpic_pairs_level": row.get("Top CPIC Pairs Level", "").strip() or None,
            "clinical_annotation_level": row.get("Top Clinical Annotation Level", "").strip() or None,
            "dosing_guideline_sources": row.get("Dosing Guideline Sources", "").strip() or None,
            "has_dosing_info": row.get("Label Has Dosing Info", "").strip().lower() == "yes",
        }

    def has_cpic_guideline(self, gene: str) -> bool:
        """Check if gene has a CPIC dosing guideline from genes.tsv."""
        row = self._genes.get(gene.strip())
        if not row:
            return False
        return row.get("Has CPIC Dosing Guideline", "").strip().lower() == "yes"

    def is_drug_supported(self, drug: str) -> bool:
        """
        Check if a drug has Level 1A or 1B CPIC guidance.

        A drug is considered "supported" if ANY of the following are true:
        1. Clinical variants dataset has 1A/1B annotations for this drug
        2. Gene-drug relationships dataset confirms this drug with 1A/1B evidence

        This aggregates across ALL genes, so warfarin is supported if it has
        1A/1B evidence for ANY of CYP2C9, VKORC1, CYP4F2, etc.

        Args:
            drug: Drug name (case-insensitive)

        Returns:
            True if drug has Level 1A or 1B evidence in PharmGKB
        """
        drug_normalized = drug.strip().lower()

        # Check 1: Clinical variants dataset
        # Look for ANY variant annotation with evidence_level 1A/1B for this drug
        for variant_row in self._clinical_variants:
            chemicals = variant_row.get("chemicals", "").strip()
            if not chemicals:
                continue

            # Check if this drug appears in the chemicals list (case-insensitive)
            chemical_list = [c.strip().lower() for c in chemicals.split(",")]
            if drug_normalized in chemical_list:
                # Check evidence level
                evidence_level = variant_row.get("level of evidence", "").strip()
                if evidence_level in ("1A", "1B"):
                    return True

        # Check 2: Gene-drug evidence index
        # Check if ANY gene has 1A/1B evidence for this drug
        for (gene, drug_key), evidence_level in self._gene_drug_evidence_index.items():
            if drug_key == drug_normalized and evidence_level in ("1A", "1B"):
                return True

        # No Level 1A or 1B evidence found
        return False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _level_rank(level: str) -> int:
    """Lower rank = stronger evidence."""
    ranks = {"1A": 1, "1B": 2, "2A": 3, "2B": 4, "3": 5, "4": 6}
    return ranks.get(level.strip(), 99)


def _classify_association(
    confirmed: bool,
    evidence_level: str,
    evidence_types: List[str],
    raw_associations: set,
) -> str:
    """
    Deterministic association classification based on evidence quality.

    Decision tree:
    1. If not confirmed → "unconfirmed"
    2. If conflicting evidence in raw data → "conflicting"
    3. If level 1A/1B + has guideline → "established"
    4. If level 2A/2B → "moderate"
    5. If level 3 and multiple evidence types → "emerging"
    6. Otherwise → "limited"

    Args:
        confirmed: Whether gene-drug pair is confirmed in dataset
        evidence_level: Highest evidence level (1A, 1B, 2A, 2B, 3, 4)
        evidence_types: List of evidence type strings
        raw_associations: Set of raw association values from dataset

    Returns:
        Classification label: established/moderate/emerging/limited/conflicting/unconfirmed
    """
    # Priority 1: Unconfirmed
    if not confirmed:
        return "unconfirmed"

    # Priority 2: Conflicting evidence (if raw data contains contradictions)
    # Check if we have both "associated" and "not associated"
    has_associated = "associated" in raw_associations
    has_not_associated = "not associated" in raw_associations
    if has_associated and has_not_associated:
        return "conflicting"

    # Priority 3: High-certainty established association
    # Level 1A/1B + guideline annotation = established
    has_guideline = any("Guideline" in et for et in evidence_types)
    if evidence_level in ("1A", "1B") and has_guideline:
        return "established"

    # Priority 4: Moderate evidence
    if evidence_level in ("2A", "2B"):
        return "moderate"

    # Priority 5: Emerging evidence (level 3 with multiple sources)
    if evidence_level == "3" and len(evidence_types) >= 3:
        return "emerging"

    # Default: Limited evidence
    return "limited"


def harmonize_annotation_associations(
    clinical_annotations: List[Dict],
    top_level_association: str,
) -> List[Dict]:
    """
    Harmonize clinical annotation associations with top-level classification.

    Enforces hierarchical consistency:
    - If top-level = "established" → annotations should be "supporting" or "established"
    - If top-level = "moderate" → annotations should be "supporting" or "moderate"
    - If top-level = "conflicting" → keep raw values (conflict is expected)

    Decision tree:
    1. If top_level in ["established", "moderate", "emerging", "limited"]:
       - Map "associated" → "supporting"
       - Map "ambiguous" → "supporting"
       - Keep "not associated" as-is (genuine non-association)
    2. If top_level = "conflicting":
       - Keep raw values (conflict evidence is expected)
    3. If top_level = "unconfirmed" or "not found":
       - Keep raw values

    Args:
        clinical_annotations: List of clinical annotation dicts with "association" field
        top_level_association: Top-level classification (established/moderate/etc.)

    Returns:
        Harmonized list of clinical annotations (creates new dicts, does not modify input)
    """
    if not clinical_annotations:
        return []

    # Associations that require harmonization
    harmonized_associations = {"established", "moderate", "emerging", "limited"}

    if top_level_association not in harmonized_associations:
        # No harmonization needed for conflicting/unconfirmed/not found
        return clinical_annotations

    # Harmonize: normalize raw associations to align with top-level
    harmonized = []
    for ann in clinical_annotations:
        # Create a copy to avoid modifying input
        harmonized_ann = dict(ann)

        raw_assoc = ann.get("association", "").lower()

        # Harmonization mapping
        if raw_assoc in ("associated", "ambiguous"):
            # These support the top-level classification
            harmonized_ann["association"] = "supporting"
        elif raw_assoc == "not associated":
            # Keep as-is (genuine non-association)
            harmonized_ann["association"] = "not associated"
        else:
            # Unknown/empty - mark as supporting if we have evidence type
            if ann.get("evidence_type"):
                harmonized_ann["association"] = "supporting"

        harmonized.append(harmonized_ann)

    return harmonized


def get_pharmgkb_loader(data_dir: Optional[Path] = None) -> PharmGKBLoader:
    """Factory / singleton accessor."""
    return PharmGKBLoader.get_instance(data_dir)
