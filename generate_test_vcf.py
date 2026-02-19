"""
generate_test_vcf.py
====================
Development utility to generate synthetic VCF test files.

Usage:
    python generate_test_vcf.py                  # writes test_data/ directory
    python generate_test_vcf.py --out my_dir     # custom output directory

NOT a deployed project asset - this is for local testing only.
The generated .vcf files are user-uploaded inputs for the VCF parser.
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

# ---------------------------------------------------------------------------
# Known pharmacogenomic variants (hg38 coords)
# ---------------------------------------------------------------------------
GENE_VARIANTS = {
    "CYP2D6": [
        ("chr22", 42522613, "rs3892097", "A", "G", "GENE=CYP2D6;STAR=*4"),
        ("chr22", 42523943, "rs5030655", "delA", ".", "GENE=CYP2D6;STAR=*6"),
        ("chr22", 42524175, "rs16947",   "G", "A", "GENE=CYP2D6;STAR=*2"),
    ],
    "CYP2C9": [
        ("chr10", 96741053, "rs1799853", "C", "T", "GENE=CYP2C9;STAR=*2"),
        ("chr10", 96748850, "rs1057910", "A", "C", "GENE=CYP2C9;STAR=*3"),
    ],
    "CYP2C19": [
        ("chr10", 96522463, "rs4244285", "G", "A", "GENE=CYP2C19;STAR=*2"),
        ("chr10", 96535173, "rs4986893", "G", "A", "GENE=CYP2C19;STAR=*3"),
        ("chr10", 96541616, "rs12248560","C", "T", "GENE=CYP2C19;STAR=*17"),
    ],
    "SLCO1B1": [
        ("chr12", 21331549, "rs4149056", "T", "C", "GENE=SLCO1B1;STAR=*5"),
        ("chr12", 21239145, "rs2306283", "A", "G", "GENE=SLCO1B1;STAR=*1B"),
    ],
    "TPMT": [
        ("chr6",  18143955, "rs1142345", "T", "C", "GENE=TPMT;STAR=*3C"),
        ("chr6",  18130918, "rs1800460", "C", "A", "GENE=TPMT;STAR=*3B"),
        ("chr6",  18128556, "rs1800462", "G", "C", "GENE=TPMT;STAR=*2"),
    ],
    "DPYD": [
        ("chr1",  97981395, "rs3918290", "C", "T", "GENE=DPYD;STAR=*2A"),
        ("chr1",  97915614, "rs55886062","A", "T", "GENE=DPYD"),
        ("chr1",  97450058, "rs67376798","T", "A", "GENE=DPYD"),
    ],
}

ZYGOSITIES = ["0/1", "1/1", "0/0", "1/0", "0|1", "1|1", "0|0"]

NOISE_VARIANTS = [
    ("chr1",  925952,  "rs2799066", "G", "A", "DP=30"),
    ("chr3",  12393541,"rs145536",  "C", "T", "DP=20"),
    ("chr5",  88888888,"rs123456",  "A", "G", "DP=45"),
    ("chr7",  42000000,"rs987654",  "T", "C", "DP=12"),
    ("chr11", 65400000,"rs111111",  "G", "T", "DP=60"),
]


def _header(patient_id: str, vcf_version: str = "VCFv4.2") -> str:
    return "\n".join([
        f"##fileformat={vcf_version}",
        f"##SAMPLE=<ID={patient_id}>",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">',
        '##INFO=<ID=STAR,Number=1,Type=String,Description="Star allele">',
        '##INFO=<ID=RS,Number=1,Type=String,Description="dbSNP RS ID">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{patient_id}",
    ])


def _row(chrom: str, pos: int, rsid: str, ref: str, alt: str,
         info: str, gt: str, qual: float = 100.0, flt: str = "PASS") -> str:
    return f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\tGT\t{gt}"


# ---------------------------------------------------------------------------
# VCF file generators
# ---------------------------------------------------------------------------

def make_comprehensive_80(patient_id: str = "PATIENT_COMP_80") -> str:
    """80-variant VCF covering all 6 pharmacogenes (10 variants each) + noise."""
    random.seed(42)
    rows = [_header(patient_id)]
    for gene, variants in GENE_VARIANTS.items():
        base_variants = variants * 4   # repeat to fill ~10 rows per gene
        for i in range(10):
            chrom, pos, rsid, ref, alt, info = base_variants[i % len(base_variants)]
            # Jitter position slightly to make rows unique
            pos_jitter = pos + i * 13
            gt = random.choice(ZYGOSITIES)
            rows.append(_row(chrom, pos_jitter, rsid, ref, alt, info, gt))
    # 20 noise rows (non-pharmacogene regions)
    for i in range(20):
        chrom, pos, rsid, ref, alt, info = NOISE_VARIANTS[i % len(NOISE_VARIANTS)]
        pos_jitter = pos + i * 7
        rows.append(_row(chrom, pos_jitter, rsid, ref, alt, info, random.choice(ZYGOSITIES)))
    return "\n".join(rows) + "\n"


def make_manual_test(patient_id: str = "PATIENT_MANUAL_TEST") -> str:
    """Small hand-crafted VCF used for quick manual inspection."""
    rows = [_header(patient_id)]
    # CYP2C19 het
    rows.append(_row("chr22", 96541616, ".", "G", "A", "GENE=CYP2C19", "0/1", qual=100.0, flt="PASS"))
    # CYP2D6 hom-alt
    rows.append(_row("chr22", 42522613, "rs3892097", "A", "G", "GENE=CYP2D6;STAR=*4", "1/1", qual=99.0, flt="PASS"))
    # CYP2C9 het phased
    rows.append(_row("chr10", 96741053, "rs1799853", "C", "T", "GENE=CYP2C9;STAR=*2", "0|1", qual=95.5, flt="PASS"))
    # DPYD missing GT
    rows.append(_row("chr1",  97981395, "rs3918290", "C", "T", "GENE=DPYD;STAR=*2A", "./.", qual=30.0, flt="."))
    # Noise
    rows.append(_row("chr1", 925952, "rs2799066", "G", "A", "DP=30", "0/1"))
    return "\n".join(rows) + "\n"


def make_single_gene(gene: str, patient_id: str = "PATIENT_SINGLE",
                     gt: str = "0/1") -> str:
    """One-variant VCF for a specific gene - used for unit-level tests."""
    rows = [_header(patient_id)]
    chrom, pos, rsid, ref, alt, info = GENE_VARIANTS[gene][0]
    rows.append(_row(chrom, pos, rsid, ref, alt, info, gt))
    return "\n".join(rows) + "\n"


def make_empty_body(patient_id: str = "PATIENT_EMPTY") -> str:
    """VCF with no data rows (edge case)."""
    return _header(patient_id) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic pharmacogenomics VCF test files.")
    parser.add_argument("--out", default="test_data", help="Output directory (default: test_data)")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    files: dict[str, str] = {
        "comprehensive_test_80.vcf": make_comprehensive_80(),
        "manual_test.vcf": make_manual_test(),
        "empty_body.vcf": make_empty_body(),
    }
    # One file per gene per zygosity
    for gene in GENE_VARIANTS:
        for label, gt in [("het", "0/1"), ("hom_alt", "1/1"), ("hom_ref", "0/0")]:
            fname = f"single_{gene.lower()}_{label}.vcf"
            files[fname] = make_single_gene(gene, patient_id=f"PT_{gene}_{label.upper()}", gt=gt)

    for fname, content in files.items():
        path = out_dir / fname
        path.write_text(content, encoding="utf-8")
        print(f"  wrote {path}  ({len(content.splitlines())} lines)")

    print(f"\nDone — {len(files)} VCF files written to '{out_dir}/'")


if __name__ == "__main__":
    main()
