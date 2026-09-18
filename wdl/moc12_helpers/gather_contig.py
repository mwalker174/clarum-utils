#!/usr/bin/env python3
"""gather_contig.py -- WDL helper: per-shard join outputs -> one VCF and one summary per contig.

The join writes one annotated VCF plus one summary JSON per input shard. This concatenates them in
filename order (shard URIs arrive coordinate-ordered from group_by_contig) and folds the per-shard
tallies into a single auditable summary. Header comes from the first part, which is the same header
every part carries because they all come from the same script run.
"""
import argparse
import glob
import json
import os
import sys

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True, help="directory holding the per-shard join outputs")
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--dataset", default="joint")
    ap.add_argument("--maf", type=float, default=0.01)
    args = ap.parse_args()

    parts = sorted(glob.glob(os.path.join(args.out_dir, "*.g41.vcf.gz")))
    if not parts:
        sys.stderr.write(f"FATAL no per-shard annotated VCFs in {args.out_dir}\n")
        return 3
    # NB the `*` must be allowed to match nothing: the join script writes the fixed name
    # `gnomad_join_summary.json`, which a pattern requiring a leading dot would silently miss -- and a
    # missing-summary set reads as "0 PASS / 0 FAIL", which is how an empty tally becomes a fake result.
    summaries = sorted(glob.glob(os.path.join(args.out_dir, "*gnomad_join_summary.json")))
    if not summaries:
        sys.stderr.write(f"FATAL no join summaries in {args.out_dir}; refusing to report zero verdicts\n")
        return 6

    agg = {"records": 0, "found_in_gnomad": 0, "absent_from_gnomad": 0, "unmatched_allele": 0,
           "PASS": 0, "FAIL": 0, "NOT MEASURED": 0, "vep_cache_empty_gnomad_af_records": 0,
           "of_those_now_joined": 0}
    shard_rows = []
    for p in summaries:
        with open(p) as fh:
            doc = json.load(fh)
        for inp in doc.get("inputs", []):
            shard_rows.append({"input": os.path.basename(inp.get("input", "")),
                               "records": inp.get("records"),
                               "found_in_gnomad": inp.get("found_in_gnomad"),
                               "absent_from_gnomad": inp.get("absent_from_gnomad"),
                               "unmatched_allele": inp.get("unmatched_allele"),
                               "verdict": inp.get("verdict"),
                               "vep_cache_empty_gnomad_af_records": inp.get("vep_cache_empty_gnomad_af_records"),
                               "of_those_now_joined": inp.get("of_those_now_joined")})
            for k in ("records", "found_in_gnomad", "absent_from_gnomad", "unmatched_allele",
                      "vep_cache_empty_gnomad_af_records", "of_those_now_joined"):
                agg[k] += inp.get(k) or 0
            for k in ("PASS", "FAIL", "NOT MEASURED"):
                agg[k] += (inp.get("verdict") or {}).get(k, 0)

    first = pysam.VariantFile(parts[0])
    out = pysam.VariantFile(args.out_vcf, "wz", header=first.header)
    n = 0
    for p in parts:
        for rec in pysam.VariantFile(p):
            out.write(rec)
            n += 1
    out.close()

    if agg["records"] and n != agg["records"]:  # counts disagree across artifacts -> say so loudly
        sys.stderr.write(f"WARNING record count mismatch: concat wrote {n}, summaries said {agg['records']}\n")

    doc = {"contig": args.contig, "dataset": args.dataset, "maf_threshold": args.maf,
           "parts": [os.path.basename(p) for p in parts], "records_written": n,
           "summaries_found": len(summaries),
           "shard_summaries": shard_rows, "aggregate": agg}
    with open(args.out_summary, "w") as fh:
        json.dump(doc, fh, indent=2)
    print(f"{args.contig}: wrote {n} records from {len(parts)} parts; "
          f"{agg['PASS']} PASS / {agg['FAIL']} FAIL / {agg['NOT MEASURED']} NOT MEASURED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
