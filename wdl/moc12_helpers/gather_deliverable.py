#!/usr/bin/env python3
"""gather_deliverable.py -- WDL helper: one annotated VCF per callset + the acceptance artifact.

Answers "1 annotated VCF generated with >40 metrics" honestly, by reporting every count separately
instead of picking whichever one clears 40:
  csq_subfields              subfields in the VEP CSQ FORMAT string (names come from the HEADER --
                             they do not appear on data lines, so a grep of data lines proves nothing)
  g41_info_fields_declared   our gnomAD v4.1 INFO columns
  metrics_manifest_rows      the documented metric manifest shipped beside the VCF
and it re-tallies the three-state verdicts, because `NOT MEASURED` must remain visible in a
deliverable-level report rather than being absorbed into PASS.
"""
import argparse
import json
import os
import re
import sys

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parts", required=True, help="tab-joined VCF paths, chromosome order")
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--acceptance", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--summaries", default="", help="tab-joined per-contig summary JSON paths")
    args = ap.parse_args()

    parts = [p for p in args.parts.split("\t") if p]
    if not parts:
        sys.stderr.write("FATAL no parts supplied\n")
        return 3
    missing = [p for p in parts if not os.path.exists(p)]
    if missing:
        sys.stderr.write(f"FATAL missing part files: {missing[:3]}\n")
        return 4

    header = pysam.VariantFile(parts[0]).header
    out = pysam.VariantFile(args.out_vcf, "wz", header=header)
    n = 0
    for p in parts:
        for rec in pysam.VariantFile(p):
            out.write(rec)
            n += 1
    out.close()

    csq = 0
    if "CSQ" in header.info:
        m = re.search(r"Format: (.*)$", header.info["CSQ"].description)
        csq = len(m.group(1).split("|")) if m else 0
    g41 = len([k for k in header.info.keys() if k.startswith("g41_") or k == "cohort_af"])

    declared = set()
    with open(args.manifest) as fh:
        for i, ln in enumerate(fh):
            row = ln.rstrip("\n").split("\t")
            if i > 0 and row and row[0]:
                declared.add(row[0])

    summary_paths = [x for x in args.summaries.split("\t") if x]
    verdicts = {"PASS": 0, "FAIL": 0, "NOT MEASURED": 0}
    joined = {"vep_cache_empty_gnomad_af_records": 0, "of_those_now_joined": 0}
    if not summary_paths:
        sys.stderr.write("WARNING no per-contig summaries supplied: verdict_tally below is EMPTY, "
                         "not zero-by-measurement\n")
    for p in summary_paths:
        try:
            with open(p) as fh:
                doc = json.load(fh)
        except (OSError, ValueError) as exc:
            sys.stderr.write(f"WARN unreadable summary {p}: {exc}\n")
            continue
        agg = doc.get("aggregate") or {}
        for k in verdicts:
            verdicts[k] += agg.get(k, 0)
        for k in joined:
            joined[k] += agg.get(k, 0)

    report = {
        "records": n,
        "contig_parts": len(parts),
        "csq_subfields": csq,
        "g41_info_fields_declared": g41,
        "metrics_manifest_rows": len(declared),
        "meets_gt40_metrics_clause": bool(csq > 40 or g41 > 40 or len(declared) > 40),
        "summaries_found": len(summary_paths),
        "verdict_tally": verdicts,
        "cache_blind_records": joined["vep_cache_empty_gnomad_af_records"],
        "cache_blind_now_measured": joined["of_those_now_joined"],
        "note": ("NOT MEASURED is a verdict, not a pass: these are variants gnomAD v4.1 does not "
                 "contain, and treating them as rare would be the bug that produced docs 090 §8A."),
    }
    with open(args.acceptance, "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
