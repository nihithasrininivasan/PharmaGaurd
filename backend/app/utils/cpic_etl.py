"""
CPIC ETL (Extract, Transform, Load) script.
Processes official CPIC Excel files and generates optimized JSON lookup structures.

Handles 6 data types per gene:
  1. Allele definition tables   → allele_definitions, variant_to_allele, rsid_map
  2. Diplotype-Phenotype tables → phenotype_map
  3. Allele functionality       → activity_scores, allele_functions
  4. Population frequencies     → frequencies
  5. Drug CDS alerts            → drug recommendations
  6. Allele definition metadata → positions (GRCh38)
"""

import pandas as pd
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import re
import logging

logger = logging.getLogger(__name__)


class CPICDataProcessor:
    """Processes CPIC allele definition tables into structured lookup data."""

    def __init__(self, cpic_data_dir: str, output_file: str):
        self.cpic_data_dir = Path(cpic_data_dir)
        self.output_file = Path(output_file)
        self.data = {
            "genes": {},
            "drugs": {},
            "metadata": {
                "version": "2.0",
                "source": "CPIC",
                "genome_build": "GRCh38",
            }
        }

    def process_all_genes(self):
        """Process all gene directories in the CPIC data directory."""
        gene_dirs = sorted(
            d for d in self.cpic_data_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        )

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
                import traceback; traceback.print_exc()

    def process_gene(self, gene_symbol: str, allele_file: Path):
        """Process a single gene's allele definition table."""
        try:
            df = pd.read_excel(allele_file, sheet_name='Alleles', header=None)
        except Exception as e:
            print(f"  Error reading allele table for {gene_symbol}: {e}")
            return

        # ---------------------------------------------------------------
        # Parse header rows
        # Row 0: Gene name
        # Row 1: Nucleotide changes
        # Row 2: Protein effect
        # Row 3: Genomic positions (GRCh38 NC_ coords)
        # Row 4: RefSeqGene positions
        # Row 5: rsIDs
        # Row 6: Label row (e.g. "CYP2C19 Allele")
        # Row 7+: *38 reference row, then allele definitions
        # ---------------------------------------------------------------
        nucleotide_changes = df.iloc[1, 1:].tolist()
        genomic_positions = df.iloc[3, 1:].tolist()
        rsids = df.iloc[5, 1:].tolist()

        # Build position map and rsID index
        position_map = {}
        rsid_map = {}  # rsid → position

        for col_idx, (nuc_change, gen_pos, rsid) in enumerate(
            zip(nucleotide_changes, genomic_positions, rsids), start=1
        ):
            if pd.notna(nuc_change):
                position_info = self._parse_position(gen_pos, nuc_change, rsid)
                if position_info:
                    position_map[col_idx] = position_info
                    # Build rsID → position mapping
                    if position_info.get("rsid"):
                        rs = position_info["rsid"]
                        if rs.startswith("rs"):
                            rsid_map[rs] = position_info["pos"]

        # Find allele start row
        # Most genes use "*1" but DPYD uses "Reference" and "c.XXX" notation
        allele_start_row = None
        uses_star_notation = True

        for idx, val in enumerate(df.iloc[:, 0]):
            val_str = str(val).strip()
            if val_str == "*1" or val_str.startswith("*1 "):
                allele_start_row = idx
                break

        # Try *<number> pattern (e.g. *38 reference rows)
        if allele_start_row is None:
            for idx, val in enumerate(df.iloc[:, 0]):
                val_str = str(val).strip()
                if re.match(r'^\*\d+', val_str):
                    allele_start_row = idx
                    break

        # DPYD fallback: look for "Reference" row
        if allele_start_row is None:
            for idx, val in enumerate(df.iloc[:, 0]):
                val_str = str(val).strip()
                if val_str.lower() == "reference":
                    allele_start_row = idx
                    uses_star_notation = False
                    break

        if allele_start_row is None:
            print(f"  Warning: Could not find allele definitions for {gene_symbol}")
            return

        # Extract allele definitions
        allele_definitions = {}
        variant_to_allele = {}

        for row_idx in range(allele_start_row, len(df)):
            allele_name = str(df.iloc[row_idx, 0]).strip()
            if not allele_name or allele_name == 'nan':
                continue

            # Filter valid allele names
            if uses_star_notation:
                if not allele_name.startswith('*'):
                    continue
                allele_name = allele_name.split()[0]
            else:
                # DPYD: accept "Reference", "c.XXX" patterns (may include " (*NX)" suffixes)
                if allele_name.lower() == "reference":
                    allele_name = "Reference"
                elif not (allele_name.startswith('c.') or allele_name.startswith('*')):
                    continue
                # Extract star allele alias if present: "c.85T>C (*9A)" → store both
                star_match = re.search(r'\(\*(\S+)\)', allele_name)
                if star_match:
                    star_alias = "*" + star_match.group(1)
                    # Use the star alias as primary name for consistency
                    allele_name = star_alias

            variants = []
            for col_idx, pos_info in position_map.items():
                if col_idx >= len(df.columns):
                    continue
                cell_value = df.iloc[row_idx, col_idx]

                if pd.isna(cell_value):
                    continue

                variant_key = self._parse_variant_call(cell_value, pos_info, allele_name)
                if variant_key:
                    variants.append(variant_key)

                    if variant_key not in variant_to_allele:
                        variant_to_allele[variant_key] = []
                    if allele_name not in variant_to_allele[variant_key]:
                        variant_to_allele[variant_key].append(allele_name)

            if variants or allele_name in ("*1", "Reference"):
                allele_definitions[allele_name] = variants

        # Store processed data
        self.data["genes"][gene_symbol] = {
            "allele_definitions": allele_definitions,
            "variant_to_allele": variant_to_allele,
            "positions": {str(k): v for k, v in position_map.items()},
            "rsid_map": rsid_map,
            "phenotype_map": {},
            "activity_scores": {},
            "allele_functions": {},
            "frequencies": {},
        }

        print(f"  Processed {len(allele_definitions)} alleles, {len(rsid_map)} rsIDs")

        gene_dir = allele_file.parent

        # ---- Process Diplotype-Phenotype table ----
        phenotype_file = gene_dir / f"{gene_symbol}_Diplotype_Phenotype_Table.xlsx"
        if phenotype_file.exists():
            self._process_phenotypes(gene_symbol, phenotype_file)
        else:
            print(f"  Warning: No phenotype table for {gene_symbol}")

        # ---- Process Allele Functionality ----
        func_file = gene_dir / f"{gene_symbol}_allele_functionality_reference.xlsx"
        if func_file.exists():
            self._process_allele_functionality(gene_symbol, func_file)
        else:
            print(f"  Warning: No functionality table for {gene_symbol}")

        # ---- Process Population Frequencies ----
        freq_file = gene_dir / f"{gene_symbol}_frequency_table.xlsx"
        if freq_file.exists():
            self._process_frequencies(gene_symbol, freq_file)

        # ---- Process Drug Alerts (CDS) ----
        for alert_file in gene_dir.glob("*_Pre_and_Post_Test_Alerts.xlsx"):
            self._process_drug_alerts(gene_symbol, alert_file)

    # ===== Diplotype-Phenotype Parsing =====

    def _process_phenotypes(self, gene_symbol: str, file_path: Path):
        """
        Parse Diplotype-Phenotype table.

        Columns (consistent across all genes):
          0: {GENE} Diplotype
          1: Activity Score
          2: Coded Diplotype/Phenotype Summary
          3: EHR Priority Notation
        """
        try:
            df = pd.read_excel(file_path)
            cols = df.columns.tolist()

            # Find columns by pattern matching
            dip_col = next((c for c in cols if "Diplotype" in c), None)
            pheno_col = next((c for c in cols if "Phenotype Summary" in c or "Phenotype" in c), None)
            score_col = next((c for c in cols if "Activity Score" in c), None)

            if not dip_col or not pheno_col:
                # Fallback: use positional
                dip_col = cols[0]
                pheno_col = cols[2] if len(cols) > 2 else cols[1]

            phenotype_map = {}
            activity_score_map = {}

            for _, row in df.iterrows():
                dip_raw = str(row[dip_col]).strip()
                pheno_raw = str(row[pheno_col]).strip()

                if dip_raw in ('nan', '') or pheno_raw in ('nan', ''):
                    continue

                # Extract phenotype from "CYP2C19 Normal Metabolizer" → "Normal Metabolizer"
                phenotype = pheno_raw
                if phenotype.startswith(gene_symbol + " "):
                    phenotype = phenotype[len(gene_symbol) + 1:]

                phenotype_map[dip_raw] = phenotype

                # Store activity score if available
                if score_col and pd.notna(row.get(score_col)):
                    try:
                        score = float(row[score_col])
                        activity_score_map[dip_raw] = score
                    except (ValueError, TypeError):
                        pass  # "n/a" or other non-numeric

            self.data["genes"][gene_symbol]["phenotype_map"] = phenotype_map
            if activity_score_map:
                self.data["genes"][gene_symbol]["diplotype_activity_scores"] = activity_score_map

            print(f"  Loaded {len(phenotype_map)} phenotype mappings")

        except Exception as e:
            print(f"  Error processing phenotypes for {gene_symbol}: {e}")
            import traceback; traceback.print_exc()

    # ===== Allele Functionality Parsing =====

    def _process_allele_functionality(self, gene_symbol: str, file_path: Path):
        """
        Parse Allele Functionality Reference table.

        Row 0: Gene label (skip)
        Row 1: Headers:
          - Allele/cDNA/rsID
          - Activity Value (Optional)
          - Allele Clinical Functional Status (Required)
          - ...
        """
        try:
            # Skip first row (gene label), use second row as header
            df = pd.read_excel(file_path, header=1)
            cols = df.columns.tolist()

            # Find columns by pattern
            allele_col = next((c for c in cols if c.startswith("Allele")), None)
            activity_col = next((c for c in cols if "Activity Value" in c), None)
            func_col = next((c for c in cols if "Clinical Functional Status" in c), None)

            if not allele_col:
                print(f"  Error: No allele column found in {file_path.name}")
                return

            activity_scores = {}
            allele_functions = {}

            for _, row in df.iterrows():
                allele_raw = str(row[allele_col]).strip()
                if allele_raw in ('nan', '') or not allele_raw.startswith('*'):
                    continue

                # Clean allele name
                allele = allele_raw.split()[0] if ' ' in allele_raw else allele_raw

                # Activity score
                if activity_col and pd.notna(row.get(activity_col)):
                    try:
                        score = float(row[activity_col])
                        activity_scores[allele] = score
                    except (ValueError, TypeError):
                        pass  # Some are "n/a"

                # Function status
                if func_col and pd.notna(row.get(func_col)):
                    func_status = str(row[func_col]).strip()
                    if func_status and func_status != 'nan':
                        allele_functions[allele] = func_status

            self.data["genes"][gene_symbol]["activity_scores"] = activity_scores
            self.data["genes"][gene_symbol]["allele_functions"] = allele_functions

            print(f"  Loaded {len(activity_scores)} activity scores, {len(allele_functions)} function statuses")

        except Exception as e:
            print(f"  Error processing functionality for {gene_symbol}: {e}")
            import traceback; traceback.print_exc()

    # ===== Population Frequency Parsing =====

    def _process_frequencies(self, gene_symbol: str, file_path: Path):
        """
        Parse Population Frequency table.

        Row 0: Label row (skip)
        Row 1: Headers — allele name in col 0, then population names
        """
        try:
            df = pd.read_excel(file_path, header=1)
            cols = df.columns.tolist()

            if len(cols) < 2:
                return

            allele_col = cols[0]  # e.g. "CYP2C19 allele"
            population_cols = cols[1:]  # e.g. "European", "East Asian", ...

            frequencies = {}

            for pop in population_cols:
                pop_name = str(pop).strip()
                if pop_name in ('nan', ''):
                    continue
                frequencies[pop_name] = {}

            for _, row in df.iterrows():
                allele_raw = str(row[allele_col]).strip()
                if allele_raw in ('nan', '') or not allele_raw.startswith('*'):
                    continue

                allele = allele_raw.split()[0]

                for pop in population_cols:
                    pop_name = str(pop).strip()
                    if pop_name in ('nan', '') or pop_name not in frequencies:
                        continue
                    val = row.get(pop)
                    if pd.notna(val):
                        try:
                            freq_val = float(val)
                            if 0 <= freq_val <= 1:
                                frequencies[pop_name][allele] = round(freq_val, 6)
                        except (ValueError, TypeError):
                            pass

            # Remove empty populations
            frequencies = {k: v for k, v in frequencies.items() if v}

            self.data["genes"][gene_symbol]["frequencies"] = frequencies
            total_pop = len(frequencies)
            total_alleles = sum(len(v) for v in frequencies.values())
            print(f"  Loaded frequencies: {total_pop} populations, {total_alleles} total entries")

        except Exception as e:
            print(f"  Error processing frequencies for {gene_symbol}: {e}")
            import traceback; traceback.print_exc()

    # ===== Drug Alert Parsing =====

    def _process_drug_alerts(self, gene_symbol: str, file_path: Path):
        """
        Parse Pre and Post Test Alerts for drug-specific recommendations.

        Columns vary by gene but follow patterns:
          - Drug Ordered
          - {GENE} Phenotype  OR  {GENE} Activity Score  (col 1)
          - CDS Context, Relative to Genetic Testing
          - CDS Alert Text
        """
        try:
            df = pd.read_excel(file_path)
            cols = df.columns.tolist()

            # Find columns
            drug_col = next((c for c in cols if "Drug" in c), cols[0])
            alert_col = next((c for c in cols if "Alert Text" in c), None)
            context_col = next((c for c in cols if "Context" in c), None)

            # Phenotype/score column is the second one
            pheno_col = cols[1] if len(cols) > 1 else None

            if not alert_col:
                print(f"  Warning: No Alert Text column in {file_path.name}")
                return

            # Detect whether this uses phenotype names or activity scores
            uses_activity_scores = pheno_col and "Activity Score" in str(pheno_col)

            for _, row in df.iterrows():
                drug = str(row[drug_col]).lower().strip()
                alert_text = str(row.get(alert_col, '')).strip()
                context = str(row.get(context_col, '')).strip() if context_col else ''

                if drug in ('nan', '') or alert_text in ('nan', ''):
                    continue

                # Get phenotype/score label
                pheno_raw = str(row.get(pheno_col, '')).strip() if pheno_col else ''

                # Skip rows without useful phenotype
                if pheno_raw in ('nan', '', 'No Result'):
                    continue

                # For multi-gene alerts (e.g., TPMT with NUDT15), handle separately
                # Check if there's also a TPMT-specific phenotype column
                tpmt_pheno = None
                if gene_symbol == "TPMT" and len(cols) >= 4:
                    tpmt_col = next((c for c in cols if "TPMT Phenotype" in c), None)
                    if tpmt_col:
                        tpmt_pheno = str(row.get(tpmt_col, '')).strip()

                # Determine phenotype label for the recommendation key
                if tpmt_pheno and tpmt_pheno not in ('nan', ''):
                    phenotype_label = tpmt_pheno
                else:
                    phenotype_label = pheno_raw

                # Strip gene prefix: "CYP2C19 Normal Metabolizer" → "Normal Metabolizer"
                if phenotype_label.startswith(gene_symbol + " "):
                    phenotype_label = phenotype_label[len(gene_symbol) + 1:]

                # Determine severity from CDS context
                if context == "No CDS" or alert_text in ('n/a', 'N/A'):
                    severity = "none"
                    alert_text = f"Standard dosing of {drug}. No clinical intervention needed."
                elif context == "Post-test" and any(
                    kw in alert_text.lower()
                    for kw in ("avoid", "contraindicated", "do not")
                ):
                    severity = "critical"
                elif context == "Post-test":
                    severity = "high"
                elif context == "Pre-test":
                    severity = "moderate"
                else:
                    severity = "moderate"

                # Initialize drug entry
                if drug not in self.data["drugs"]:
                    self.data["drugs"][drug] = {
                        "gene": gene_symbol,
                        "recommendations": {}
                    }

                self.data["drugs"][drug]["recommendations"][phenotype_label] = {
                    "risk": alert_text[:50] + "..." if len(alert_text) > 50 else alert_text,
                    "severity": severity,
                    "implication": alert_text,
                    "alert_context": context,
                    "url": "https://cpicpgx.org/guidelines/"
                }

            drug_name = file_path.stem.split("_")[0]
            recs = self.data["drugs"].get(drug_name, {}).get("recommendations", {})
            print(f"  Loaded {len(recs)} drug alert(s) from {file_path.name}")

        except Exception as e:
            print(f"  Error processing drug alerts from {file_path.name}: {e}")
            import traceback; traceback.print_exc()

    # ===== Position Parsing =====

    def _parse_position(self, genomic_pos: str, nuc_change: str, rsid: Any) -> Optional[Dict]:
        """Parse genomic position information from allele definition header."""
        if pd.isna(genomic_pos):
            return None

        # Extract position from notation like "g.42126578C>T"
        match = re.search(r'g\.(\d+)([ACGT])>([ACGT])', str(genomic_pos))
        if match:
            pos = int(match.group(1))
            ref = match.group(2)
            alt = match.group(3)

            return {
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "rsid": str(rsid).strip() if pd.notna(rsid) else None,
                "nuc_change": str(nuc_change).strip() if pd.notna(nuc_change) else None
            }

        return None

    def _parse_variant_call(self, cell_value: Any, pos_info: Dict, allele_name: str) -> Optional[str]:
        """Parse a variant call from a cell value."""
        if pd.isna(cell_value):
            return None

        cell_str = str(cell_value).strip().upper()

        if cell_str in ['NAN', '', 'NONE']:
            return None

        # If it's a single base that matches the alt allele
        if cell_str in ['A', 'C', 'G', 'T']:
            if cell_str == pos_info['alt']:
                return f"{pos_info['pos']}:{pos_info['ref']}:{pos_info['alt']}"
            # Could be the ref allele (wildtype marker) — skip
            if cell_str == pos_info['ref']:
                return None

        # Handle IUPAC ambiguity codes
        iupac_codes = {
            'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
            'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C']
        }

        if cell_str in iupac_codes:
            if pos_info['alt'] in iupac_codes[cell_str]:
                return f"{pos_info['pos']}:{pos_info['ref']}:{pos_info['alt']}"

        # Handle deletion notation (e.g., "delXXX" or structural variants)
        if 'DEL' in cell_str:
            return f"{pos_info['pos']}:{pos_info['ref']}:DEL"

        # Handle insertion
        if 'INS' in cell_str:
            return f"{pos_info['pos']}:{pos_info['ref']}:INS"

        return None

    def save(self):
        """Save processed data to JSON file."""
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(self.output_file, 'w') as f:
            json.dump(self.data, f, indent=2)

        # Summary statistics
        total_alleles = sum(
            len(g.get("allele_definitions", {}))
            for g in self.data["genes"].values()
        )
        total_phenotypes = sum(
            len(g.get("phenotype_map", {}))
            for g in self.data["genes"].values()
        )
        total_activity = sum(
            len(g.get("activity_scores", {}))
            for g in self.data["genes"].values()
        )
        total_freq = sum(
            sum(len(p) for p in g.get("frequencies", {}).values())
            for g in self.data["genes"].values()
        )

        print(f"\nSaved CPIC cache to {self.output_file}")
        print(f"  Genes:      {len(self.data['genes'])}")
        print(f"  Drugs:      {len(self.data['drugs'])}")
        print(f"  Alleles:    {total_alleles}")
        print(f"  Phenotypes: {total_phenotypes}")
        print(f"  Activity:   {total_activity}")
        print(f"  Frequencies: {total_freq}")


def main():
    """Main ETL execution."""
    backend_dir = Path(__file__).parent.parent.parent
    cpic_data_dir = backend_dir / "data" / "cpic"
    output_file = backend_dir / "data" / "cpic_cache.json"

    print(f"CPIC ETL Script v2.0")
    print(f"Source: {cpic_data_dir}")
    print(f"Output: {output_file}")
    print("-" * 80)

    processor = CPICDataProcessor(str(cpic_data_dir), str(output_file))
    processor.process_all_genes()
    processor.save()

    print("\nETL complete!")


if __name__ == "__main__":
    main()
