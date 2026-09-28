#!/usr/bin/env python3
"""gather_callset.py -- WDL helper: concatenate annotated parts into one callset VCF + acceptance JSON.

Replaces gather_deliverable.py's record-by-record pysam rewrite for full callsets: a 620-sample WGS
callset is ~47.7M records, and re-encoding every genotype through pysam single-threaded is hours of
work that `bcftools concat --threads` does in a fraction of the time. The concat itself is bcftools;
this script decides the ORDER, runs it, and then refuses to certify anything it cannot count.

Order is derived from the data, never from the input array. Cromwell's scatter-of-glob order is not
a genomic order, and a mis-ordered concat either fails or -- with --allow-overlaps -- silently
produces an unsorted file. Each part's first and last record is read, parts are sorted by
(contig index in the header, first position), and any overlap between neighbours is fatal.

Certification, all of which must hold or the exit code is non-zero:
  * output record count == sum of part record counts (from the index, not a text scan);
  * none of the --dropped INFO tags survives in the output header;
  * every per-part summary JSON parses.

Counts are reported SEPARATELY (csq_subfields, g41_info_fields, clinvar_info_fields, ...), never
as whichever one clears 40 -- same rule as gather_deliverable.py.
"""

import argparse
import json
import os
import re
import subprocess
import sys

import pysam


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write(f"FATAL {' '.join(cmd)}\n{r.stderr}\n")
        raise SystemExit(5)
    return r.stdout


def n_records(path, nonempty=True):
    """Record count from the .tbi. htslib's `tabix` writes per-contig counts; a GATK-written .tbi
    does NOT, and `bcftools index -n` then answers 0 with rc 0 (measured on the JointGenotyping
    shards, 2026-09-25). So a 0 from a file known to hold records is fatal, never a count."""
    n = int(run(["bcftools", "index", "-n", path]).strip())
    if nonempty and n == 0:
        sys.stderr.write(f"FATAL {path}: index reports 0 records for a non-empty file -- the .tbi "
                         "carries no counts; re-index with tabix\n")
        raise SystemExit(5)
    return n


def span(path):
    vf = pysam.VariantFile(path)
    contigs = list(vf.header.contigs)
    first = last = None
    for rec in vf:
        first = (contigs.index(rec.chrom), rec.pos, rec.chrom)
        break
    if first is None:
        return None
    # last record: the index lists which contigs carry records; walk only the final one
    with_data = [ln.split("\t")[0] for ln in run(["bcftools", "index", "-s", path]).splitlines() if ln]
    with_data.sort(key=contigs.index)
    tail = None
    for rec in vf.fetch(with_data[-1]):
        tail = rec
    last = (contigs.index(tail.chrom), tail.pos, tail.chrom)
    return first, last


def lines_of(path):
    with open(path) as fh:
        return [ln for ln in fh.read().split("\n") if ln]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", required=True, help="file listing part VCF paths, one per line (.tbi alongside)")
    ap.add_argument("--gnomad-summaries", required=True, help="file listing per-part gnomad_join_summary.json")
    ap.add_argument("--clinvar-summaries", required=True, help="file listing per-part clinvar_join_summary.json")
    ap.add_argument("--g41-manifest", required=True)
    ap.add_argument("--clinvar-manifest", required=True)
    ap.add_argument("--dropped", default="", help="comma-separated INFO tags that must be absent")
    ap.add_argument("--blanked-csq", default="",
                    help="comma-separated CSQ subfields that must be empty on every record")
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--acceptance", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    parts = lines_of(args.parts)
    if not parts:
        sys.stderr.write("FATAL no parts supplied\n")
        return 3

    spans, empty = [], []
    for p in parts:
        s = span(p)
        (empty if s is None else spans).append(p if s is None else (s, p))
    spans.sort(key=lambda x: (x[0][0][0], x[0][0][1]))
    for (a, pa), (b, pb) in zip(spans, spans[1:]):
        if (b[0][0], b[0][1]) <= (a[1][0], a[1][1]):
            sys.stderr.write(f"FATAL parts overlap: {pa} ends {a[1][2]}:{a[1][1]}, "
                             f"{pb} starts {b[0][2]}:{b[0][1]}\n")
            return 4
    ordered = [p for _s, p in spans]
    with open("parts.ordered.list", "w") as fh:
        fh.write("\n".join(ordered) + "\n")

    expected = sum(n_records(p) for p in ordered)
    run(["bcftools", "concat", "--threads", str(args.threads), "-f", "parts.ordered.list",
         "-Oz", "-o", args.out_vcf])
    run(["tabix", "-f", "-p", "vcf", args.out_vcf])
    got = n_records(args.out_vcf)

    header = pysam.VariantFile(args.out_vcf).header
    info = list(header.info.keys())
    dropped = [t for t in args.dropped.split(",") if t]
    survived = [t for t in dropped if t in info]

    csq_names = []
    if "CSQ" in header.info:
        m = re.search(r"Format: (.*)$", header.info["CSQ"].description)
        csq_names = m.group(1).split("|") if m else []
    csq = len(csq_names)

    # Blanked CSQ subfields: the header marker says the step ran, and a full-file scan proves it
    # left nothing behind. `bcftools query` streams CSQ only, so this is one pass at I/O speed.
    blanked = [t for t in args.blanked_csq.split(",") if t]
    blank_missing = [t for t in blanked if t not in csq_names]
    marker = any(ln.startswith("##clarumCsqBlanked=") for ln in str(header).splitlines())
    blank_survivors = 0
    if blanked and not blank_missing:
        idx = [csq_names.index(t) for t in blanked]
        q = subprocess.Popen(["bcftools", "query", "-f", "%INFO/CSQ\n", args.out_vcf],
                             stdout=subprocess.PIPE, text=True)
        for line in q.stdout:
            for entry in line.rstrip("\n").split(","):
                parts = entry.split("|")
                blank_survivors += sum(1 for i in idx if i < len(parts) and parts[i])
        if q.wait():
            sys.stderr.write("FATAL bcftools query over CSQ failed\n")
            return 5

    def manifest_rows(path):
        return len({ln.split("\t")[0] for ln in lines_of(path)[1:] if ln.split("\t")[0]})

    g41_verdict = {"PASS": 0, "FAIL": 0, "NOT MEASURED": 0}
    g41 = {"found_in_gnomad": 0, "absent_from_gnomad": 0, "unmatched_allele": 0,
           "vep_cache_empty_gnomad_af_records": 0, "of_those_now_joined": 0}
    unpublished = set()
    cv = {"found": 0, "absent_from_clinvar": 0, "unmatched_allele": 0,
          "records_with_somatic_classification": 0}
    cv_verdicts, cv_stars = {}, {}
    bad = []
    for p in lines_of(args.gnomad_summaries):
        try:
            doc = json.load(open(p))
        except (OSError, ValueError) as exc:
            bad.append(f"{p}: {exc}")
            continue
        for inp in doc.get("inputs", []):
            for k in g41_verdict:
                g41_verdict[k] += inp.get("verdict", {}).get(k, 0)
            for k in g41:
                g41[k] += inp.get(k, 0)
            for cs in (inp.get("contigs_not_published_by_gnomad") or {}).values():
                unpublished.update(cs)
    for p in lines_of(args.clinvar_summaries):
        try:
            doc = json.load(open(p))
        except (OSError, ValueError) as exc:
            bad.append(f"{p}: {exc}")
            continue
        for inp in doc.get("inputs", []):
            for k in cv:
                cv[k] += inp.get(k, 0)
            for k, v in inp.get("germline_verdicts", {}).items():
                cv_verdicts[k] = cv_verdicts.get(k, 0) + v
            for k, v in inp.get("review_stars", {}).items():
                cv_stars[k] = cv_stars.get(k, 0) + v

    report = {
        "records": got,
        "records_expected_from_parts": expected,
        "parts": len(ordered),
        "empty_parts_skipped": empty,
        "samples": len(header.samples),
        "csq_subfields": csq,
        "g41_info_fields": len([k for k in info if k.startswith("g41_")]),
        "clinvar_info_fields": len([k for k in info if k.startswith("clinvar_")]),
        "other_info_fields": len([k for k in info if not k.startswith(("g41_", "clinvar_")) and k != "CSQ"]),
        "g41_manifest_rows": manifest_rows(args.g41_manifest),
        "clinvar_manifest_rows": manifest_rows(args.clinvar_manifest),
        "dropped_info_tags": dropped,
        "dropped_info_tags_surviving": survived,
        "blanked_csq_subfields": blanked,
        "blanked_csq_marker_in_header": marker,
        "blanked_csq_nonempty_values": blank_survivors,
        "gnomad": {**g41, "verdict_maf": g41_verdict,
                   "contigs_not_published_by_gnomad": sorted(unpublished)},
        "clinvar": {**cv, "germline_verdicts": cv_verdicts, "review_stars": cv_stars},
        "unreadable_summaries": bad,
        "note": ("NOT MEASURED is a verdict, not a pass: absent from gnomAD is not rare, and "
                 "absent from ClinVar is not benign."),
    }
    with open(args.acceptance, "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))

    failures = []
    if got != expected:
        failures.append(f"output has {got} records, parts sum to {expected}")
    if survived:
        failures.append(f"dropped tags still in header: {survived}")
    if bad:
        failures.append(f"{len(bad)} unreadable summaries")
    if blank_missing:
        failures.append(f"blanked CSQ subfields not in the CSQ Format: {blank_missing}")
    if blanked and not marker:
        failures.append("no ##clarumCsqBlanked header line: the blanking step did not run")
    if blank_survivors:
        failures.append(f"{blank_survivors} non-empty values remain in blanked CSQ subfields")
    for f in failures:
        sys.stderr.write(f"REFUSING TO CERTIFY: {f}\n")
    return 6 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
