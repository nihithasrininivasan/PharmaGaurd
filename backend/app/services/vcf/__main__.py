from __future__ import annotations

import json
import sys
from pathlib import Path

from .parser import parse_vcf
from .variant_extractor import extract_pharmacogenes, extract_pharmacogenes_from_vcf_path


def main(argv: list[str]) -> int:
    if len(argv) < 2 or "--help" in argv:
        print("Usage: python -m services.vcf <path-to.vcf> [--stream] [--max-variants N]")
        return 0
    
    path = Path(argv[1])
    if not path.exists():
        print(f"File not found: {path}")
        return 2

    stream = "--stream" in argv
    max_variants = None
    if "--max-variants" in argv:
        try:
            idx = argv.index("--max-variants")
            max_variants = int(argv[idx + 1])
        except (ValueError, IndexError):
            print("Error: --max-variants requires an integer argument")
            return 2

    if stream:
        header, by_gene = extract_pharmacogenes_from_vcf_path(
            path, max_variants_kept=max_variants, max_scanned_variants=max_variants
        )
        patient_id = header.patient_id
        vcf_version = header.vcf_version
        quality_metrics = {
            "vcf_parsing_success": True,
            "variant_count_kept": sum(len(vs) for vs in by_gene.values()),
            "sample_count": len(header.samples),
            "has_column_header": True,
            "vcf_version": vcf_version or "",
            "streaming_mode": True,
        }
    else:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            parsed = parse_vcf(f, max_variants=max_variants)
        patient_id = parsed.patient_id
        vcf_version = parsed.vcf_version
        quality_metrics = parsed.quality_metrics
        by_gene = extract_pharmacogenes(parsed.variants)

    payload = {
        "patient_id": patient_id,
        "vcf_version": vcf_version,
        "quality_metrics": quality_metrics,
        "pharmacogene_variant_counts": {g: len(vs) for g, vs in by_gene.items()},
        "pharmacogene_variants": {
            g: [
                {
                    "rsid": v.rsid,
                    "star": v.star,
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "ref": v.ref,
                    "alt": v.alt,
                    "gt": v.gt,
                    "zygosity": v.zygosity,
                    "gene": v.gene,
                }
                for v in vs
            ]
            for g, vs in by_gene.items()
        },
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

