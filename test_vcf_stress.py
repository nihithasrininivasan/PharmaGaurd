"""
Stress test: run 50 synthetic VCF datasets through analyze_vcf_for_drugs
and print a summary table to the terminal.

Run with:
    python test_vcf_stress.py
"""
from __future__ import annotations
import json
import sys
import traceback
from typing import List, Tuple

from services.vcf.pharmaguard_adapter import analyze_vcf_for_drugs

# ── Known PGx variants per gene ───────────────────────────────────────────────
GENE_VARIANTS = {
    "CYP2D6": [
        ("22", 42522613, "rs3892097", "A", "G", "GENE=CYP2D6;STAR=*4"),
        ("22", 42523943, "rs5030655", "delA", ".", "GENE=CYP2D6;STAR=*6"),
        ("22", 42524175, "rs16947",   "G", "A", "GENE=CYP2D6;STAR=*2"),
    ],
    "CYP2C9": [
        ("10", 96741053, "rs1799853", "C", "T", "GENE=CYP2C9;STAR=*2"),
        ("10", 96748850, "rs1057910", "A", "C", "GENE=CYP2C9;STAR=*3"),
    ],
    "CYP2C19": [
        ("10", 96522463, "rs4244285", "G", "A", "GENE=CYP2C19;STAR=*2"),
        ("10", 96535173, "rs4986893", "G", "A", "GENE=CYP2C19;STAR=*3"),
    ],
    "SLCO1B1": [
        ("12", 21331549, "rs4149056", "T", "C", "GENE=SLCO1B1;STAR=*5"),
    ],
    "TPMT": [
        ("6",  18143955, "rs1142345", "T", "C", "GENE=TPMT;STAR=*3C"),
        ("6",  18130918, "rs1800460", "C", "A", "GENE=TPMT;STAR=*3B"),
    ],
    "DPYD": [
        ("1",  97981395, "rs3918290", "C", "T", "GENE=DPYD;STAR=*2A"),
        ("1",  97915614, "rs55886062","A", "T", "GENE=DPYD"),
    ],
}

DRUG_GENE_MAP = {
    "WARFARIN":    ["CYP2C9"],
    "CODEINE":     ["CYP2D6"],
    "CLOPIDOGREL": ["CYP2C19"],
    "SIMVASTATIN": ["SLCO1B1"],
    "AZATHIOPRINE":["TPMT"],
    "FLUOROURACIL":["DPYD"],
}

ZYGOSITIES = ["0/1", "1/1", "0/0", "1/0", "0|1", "1|1", "./.", "."]
DRUGS_ALL  = list(DRUG_GENE_MAP.keys())


def make_vcf(patient_id: str, variants: list, vcf_version: str = "VCFv4.2") -> bytes:
    lines = [
        f"##fileformat={vcf_version}",
        f"##SAMPLE=<ID={patient_id}>",
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">",
        "##INFO=<ID=STAR,Number=1,Type=String,Description=\"Star allele\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{patient_id}",
    ]
    for chrom, pos, rsid, ref, alt, info, gt in variants:
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{gt}")
    return "\n".join(lines).encode("utf-8")


def build_datasets() -> List[Tuple[str, bytes, List[str]]]:
    """Build 50 varied test datasets."""
    import itertools, random
    random.seed(42)
    datasets = []

    # 1-6: one variant per gene, het
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 1):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "0/1")])
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_single_het_{gene}", vcf, [drug]))

    # 7-12: hom-alt
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 7):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "1/1")])
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_single_homalt_{gene}", vcf, [drug]))

    # 13-18: hom-ref (no variants → Safe)
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 13):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "0/0")])
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_homref_{gene}", vcf, [drug]))

    # 19-24: missing GT (.)/(..)
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 19):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "./.")])
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_missing_gt_{gene}", vcf, [drug]))

    # 25-30: phased genotype (pipe notation)
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 25):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "0|1")])
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_phased_{gene}", vcf, [drug]))

    # 31-36: multiple variants per patient
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 31):
        rows = [(c, p, r, re, al, inf, "0/1") for c, p, r, re, al, inf in var_list]
        vcf = make_vcf(f"PT{i:03d}", rows)
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_multi_variant_{gene}", vcf, [drug]))

    # 37-42: all drugs queried at once
    for i, (gene, var_list) in enumerate(GENE_VARIANTS.items(), 37):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "0/1")])
        datasets.append((f"DS{i:02d}_alldrugs_{gene}", vcf, DRUGS_ALL))

    # 43-46: VCF 4.1 version tag
    for i, (gene, var_list) in enumerate(list(GENE_VARIANTS.items())[:4], 43):
        chrom, pos, rsid, ref, alt, info = var_list[0]
        vcf = make_vcf(f"PT{i:03d}", [(chrom, pos, rsid, ref, alt, info, "1/1")], vcf_version="VCFv4.1")
        drug = [d for d, gs in DRUG_GENE_MAP.items() if gene in gs][0]
        datasets.append((f"DS{i:02d}_vcf41_{gene}", vcf, [drug]))

    # 47-50: no variants at all (empty body)
    for i in range(47, 51):
        vcf = make_vcf(f"PT{i:03d}", [])
        datasets.append((f"DS{i:02d}_empty_body", vcf, DRUGS_ALL))

    return datasets[:50]


def main():
    datasets = build_datasets()
    passed = 0
    failed = 0
    rows = []

    lines_out = []
    lines_out.append(f"\n{'='*90}")
    lines_out.append(f"  VCF PARSER STRESS TEST  -  {len(datasets)} datasets")
    lines_out.append(f"{'='*90}")
    lines_out.append(f"{'#':<4} {'Dataset':<35} {'Drugs':<3} {'Status':<8} {'Variants':<10} {'Risk Sample'}")
    lines_out.append(f"{'-'*90}")

    for idx, (name, vcf_bytes, drugs) in enumerate(datasets, 1):
        try:
            results = analyze_vcf_for_drugs(vcf_bytes, drugs)
            total_variants = sum(len(r["pharmacogenomic_profile"]["detected_variants"]) for r in results)
            risk_sample = results[0]["risk_assessment"]["risk_label"] if results else "N/A"
            lines_out.append(f"{idx:<4} {name:<35} {len(drugs):<3} {'PASS':<8} {total_variants:<10} {risk_sample}")
            passed += 1
        except Exception as e:
            lines_out.append(f"{idx:<4} {name:<35} {len(drugs):<3} {'FAIL':<8} {'-':<10} {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1

    lines_out.append(f"{'='*90}")
    lines_out.append(f"  RESULTS: {passed} PASSED  |  {failed} FAILED  |  {len(datasets)} TOTAL")
    lines_out.append(f"{'='*90}")

    output = "\n".join(lines_out)
    print(output)

    # Write to file for easy inspection
    import pathlib
    out_path = pathlib.Path(__file__).parent / "stress_test_results.txt"
    out_path.write_text(output, encoding="utf-8")
    print(f"\nResults also saved to: {out_path}")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
