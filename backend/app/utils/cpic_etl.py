"""
CPIC ETL (Extract, Transform, Load) script.
Processes CPIC Excel files and generates optimized JSON lookup structures.
"""

import pandas as pd
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any
import re


class CPICDataProcessor:
    """Processes CPIC allele definition tables into structured lookup data."""

    def __init__(self, cpic_data_dir: str, output_file: str):
        self.cpic_data_dir = Path(cpic_data_dir)
        self.output_file = Path(output_file)
        self.data = {
            "genes": {},
            "drugs": {},
            "metadata": {
                "version": "1.0",
                "source": "CPIC"
            }
        }

    def process_all_genes(self):
        """Process all gene directories in the CPIC data directory."""
        gene_dirs = [d for d in self.cpic_data_dir.iterdir() if d.is_dir()]

        for gene_dir in gene_dirs:
            gene_symbol = gene_dir.name
            print(f"Processing {gene_symbol}...")

            try:
                allele_file = gene_dir / f"{gene_symbol}_allele_definition_table.xlsx"
                if allele_file.exists():
                    self.process_gene(gene_symbol, allele_file)
                else:
                    print(f"  Warning: No allele definition table found for {gene_symbol}")
            except Exception as e:
                print(f"  Error processing {gene_symbol}: {e}")

    def process_gene(self, gene_symbol: str, allele_file: Path):
        """Process a single gene's allele definition table."""
        # Read the Alleles sheet
        df = pd.read_excel(allele_file, sheet_name='Alleles', header=None)

        # Parse header rows to extract position information
        # Row 0: Gene name
        # Row 1: Nucleotide changes (e.g., "4214G>A")
        # Row 3: Genomic positions (e.g., "g.42126578C>T")
        # Row 5: rsIDs
        # Row 6+: Allele names and their definitions

        nucleotide_changes = df.iloc[1, 1:].tolist()
        genomic_positions = df.iloc[3, 1:].tolist()
        rsids = df.iloc[5, 1:].tolist()

        # Build position map: column index -> position info
        position_map = {}
        for col_idx, (nuc_change, gen_pos, rsid) in enumerate(zip(nucleotide_changes, genomic_positions, rsids), start=1):
            if pd.notna(nuc_change):
                position_info = self._parse_position(gen_pos, nuc_change, rsid)
                if position_info:
                    position_map[col_idx] = position_info

        # Process allele definitions (starting from row 7 or 8, depending on structure)
        # Find the row where allele names start (look for "*1")
        allele_start_row = None
        for idx, val in enumerate(df.iloc[:, 0]):
            if str(val).strip() == "*1":
                allele_start_row = idx
                break

        if allele_start_row is None:
            print(f"  Warning: Could not find allele definitions for {gene_symbol}")
            return

        # Extract allele definitions
        allele_definitions = {}
        variant_to_allele = {}

        for row_idx in range(allele_start_row, len(df)):
            allele_name = str(df.iloc[row_idx, 0]).strip()
            if not allele_name or allele_name == 'nan' or not allele_name.startswith('*'):
                continue

            # Get variant calls for this allele
            variants = []
            for col_idx, pos_info in position_map.items():
                cell_value = df.iloc[row_idx, col_idx]

                # Skip if NaN (means reference/wildtype)
                if pd.isna(cell_value):
                    continue

                # Parse the cell value
                variant_key = self._parse_variant_call(cell_value, pos_info, allele_name)
                if variant_key:
                    variants.append(variant_key)

                    # Add to variant_to_allele map
                    if variant_key not in variant_to_allele:
                        variant_to_allele[variant_key] = []
                    if allele_name not in variant_to_allele[variant_key]:
                        variant_to_allele[variant_key].append(allele_name)

            if variants or allele_name == "*1":
                allele_definitions[allele_name] = variants

        # Store processed data
        self.data["genes"][gene_symbol] = {
            "allele_definitions": allele_definitions,
            "variant_to_allele": variant_to_allele,
            "positions": position_map
        }

        print(f"  Processed {len(allele_definitions)} alleles for {gene_symbol}")

    def _parse_position(self, genomic_pos: str, nuc_change: str, rsid: Any) -> Optional[Dict]:
        """Parse genomic position information."""
        if pd.isna(genomic_pos):
            return None

        # Extract position from genomic notation (e.g., "g.42126578C>T")
        match = re.search(r'g\.(\d+)([ACGT])>([ACGT])', str(genomic_pos))
        if match:
            pos = int(match.group(1))
            ref = match.group(2)
            alt = match.group(3)

            return {
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "rsid": str(rsid) if pd.notna(rsid) else None,
                "nuc_change": str(nuc_change) if pd.notna(nuc_change) else None
            }

        return None

    def _parse_variant_call(self, cell_value: Any, pos_info: Dict, allele_name: str) -> Optional[str]:
        """Parse a variant call from a cell value."""
        if pd.isna(cell_value):
            return None

        cell_str = str(cell_value).strip().upper()

        # Handle special cases
        if cell_str in ['NAN', '']:
            return None

        # If it's a single base (the alternate allele)
        if cell_str in ['A', 'C', 'G', 'T']:
            # This is the alternate allele
            if cell_str == pos_info['alt']:
                return f"{pos_info['pos']}:{pos_info['ref']}:{pos_info['alt']}"

        # Handle IUPAC ambiguity codes (for diploid representations in some tables)
        # R = A or G, Y = C or T, S = G or C, W = A or T, K = G or T, M = A or C
        iupac_codes = {
            'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
            'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C']
        }

        if cell_str in iupac_codes:
            # This indicates heterozygosity - the allele contains the alt
            if pos_info['alt'] in iupac_codes[cell_str]:
                return f"{pos_info['pos']}:{pos_info['ref']}:{pos_info['alt']}"

        return None

    def add_phenotype_map(self, gene: str, phenotype_mapping: Dict[str, str]):
        """Add diplotype to phenotype mapping for a gene."""
        if gene in self.data["genes"]:
            self.data["genes"][gene]["phenotype_map"] = phenotype_mapping

    def add_drug_recommendations(self, drug: str, gene: str, recommendations: Dict[str, Dict]):
        """Add drug-specific recommendations."""
        if drug not in self.data["drugs"]:
            self.data["drugs"][drug] = {}

        self.data["drugs"][drug] = {
            "gene": gene,
            "recommendations": recommendations
        }

    def save(self):
        """Save processed data to JSON file."""
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(self.output_file, 'w') as f:
            json.dump(self.data, f, indent=2)

        print(f"\nSaved CPIC cache to {self.output_file}")
        print(f"Processed {len(self.data['genes'])} genes")
        print(f"Configured {len(self.data['drugs'])} drugs")


def add_manual_phenotype_mappings(processor: CPICDataProcessor):
    """
    Add phenotype mappings that are typically found in diplotype-phenotype tables
    or derived from activity scores.
    """
    # CYP2D6 phenotype mappings (based on activity scores)
    cyp2d6_phenotypes = {
        "*1/*1": "NM", "*1/*2": "NM", "*2/*2": "NM",
        "*1/*3": "IM", "*1/*4": "IM", "*1/*5": "IM", "*1/*6": "IM",
        "*2/*3": "IM", "*2/*4": "IM", "*2/*5": "IM", "*2/*6": "IM",
        "*3/*3": "PM", "*3/*4": "PM", "*3/*5": "PM", "*3/*6": "PM",
        "*4/*4": "PM", "*4/*5": "PM", "*4/*6": "PM",
        "*5/*5": "PM", "*5/*6": "PM", "*6/*6": "PM",
        "*1/*41": "IM", "*2/*41": "IM", "*4/*41": "IM",
        "*41/*41": "IM"
    }
    processor.add_phenotype_map("CYP2D6", cyp2d6_phenotypes)

    # CYP2C19 phenotype mappings
    cyp2c19_phenotypes = {
        "*1/*1": "NM", "*1/*2": "IM", "*2/*2": "PM",
        "*1/*3": "IM", "*2/*3": "PM", "*3/*3": "PM",
        "*1/*17": "RM", "*17/*17": "UM"
    }
    processor.add_phenotype_map("CYP2C19", cyp2c19_phenotypes)

    # CYP2C9 phenotype mappings
    cyp2c9_phenotypes = {
        "*1/*1": "NM", "*1/*2": "IM", "*1/*3": "IM",
        "*2/*2": "IM", "*2/*3": "PM", "*3/*3": "PM"
    }
    processor.add_phenotype_map("CYP2C9", cyp2c9_phenotypes)

    # TPMT phenotype mappings
    tpmt_phenotypes = {
        "*1/*1": "NM", "*1/*2": "IM", "*1/*3A": "IM", "*1/*3C": "IM",
        "*2/*2": "PM", "*2/*3A": "PM", "*2/*3C": "PM",
        "*3A/*3A": "PM", "*3A/*3C": "PM", "*3C/*3C": "PM"
    }
    processor.add_phenotype_map("TPMT", tpmt_phenotypes)


def add_drug_guidelines(processor: CPICDataProcessor):
    """Add drug-specific clinical guidelines."""

    # Codeine - CYP2D6
    processor.add_drug_recommendations("codeine", "CYP2D6", {
        "PM": {
            "risk": "Avoid codeine use",
            "severity": "high",
            "implication": "Reduced morphine formation, insufficient analgesia",
            "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        "IM": {
            "risk": "Use alternative analgesic or monitor closely",
            "severity": "moderate",
            "implication": "Reduced morphine formation",
            "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        "NM": {
            "risk": "Standard dosing",
            "severity": "none",
            "implication": "Normal morphine formation",
            "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        },
        "UM": {
            "risk": "Avoid codeine use",
            "severity": "high",
            "implication": "Increased morphine formation, risk of toxicity",
            "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
        }
    })

    # Warfarin - CYP2C9
    processor.add_drug_recommendations("warfarin", "CYP2C9", {
        "PM": {
            "risk": "Reduce initial dose by 50-75%",
            "severity": "high",
            "implication": "Increased bleeding risk, reduced clearance",
            "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/"
        },
        "IM": {
            "risk": "Reduce initial dose by 25-50%",
            "severity": "moderate",
            "implication": "Moderately increased bleeding risk",
            "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/"
        },
        "NM": {
            "risk": "Standard dosing",
            "severity": "none",
            "implication": "Normal warfarin metabolism",
            "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/"
        }
    })

    # Clopidogrel - CYP2C19
    processor.add_drug_recommendations("clopidogrel", "CYP2C19", {
        "PM": {
            "risk": "Use alternative antiplatelet (prasugrel, ticagrelor)",
            "severity": "high",
            "implication": "Reduced active drug formation, reduced platelet inhibition",
            "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        "IM": {
            "risk": "Consider alternative antiplatelet or higher dose",
            "severity": "moderate",
            "implication": "Moderately reduced platelet inhibition",
            "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        "NM": {
            "risk": "Standard dosing",
            "severity": "none",
            "implication": "Normal clopidogrel activation",
            "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        },
        "UM": {
            "risk": "Standard dosing",
            "severity": "low",
            "implication": "Increased platelet inhibition (may be beneficial)",
            "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
        }
    })


def main():
    """Main ETL execution."""
    # Set up paths
    backend_dir = Path(__file__).parent.parent.parent
    cpic_data_dir = backend_dir / "data" / "cpic"
    output_file = backend_dir / "data" / "cpic_cache.json"

    print(f"CPIC ETL Script")
    print(f"Source: {cpic_data_dir}")
    print(f"Output: {output_file}")
    print("-" * 80)

    # Create processor and process all genes
    processor = CPICDataProcessor(str(cpic_data_dir), str(output_file))
    processor.process_all_genes()

    # Add manual phenotype mappings
    print("\nAdding phenotype mappings...")
    add_manual_phenotype_mappings(processor)

    # Add drug guidelines
    print("Adding drug guidelines...")
    add_drug_guidelines(processor)

    # Save results
    processor.save()

    print("\nETL complete!")


if __name__ == "__main__":
    main()
