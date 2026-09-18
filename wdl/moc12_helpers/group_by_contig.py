#!/usr/bin/env python3
"""group_by_contig.py -- WDL helper: decide which callset shards share one task.

Reads a newline-delimited file of shard URIs and writes a TSV to stdout with exactly three columns:

    contig <TAB> group_index <TAB> shard_uris_joined_with_pipe

Two rules that are not arbitrary:
  * The third column is pipe-joined, never tab-joined. This script's own output is read by Cromwell's
    read_tsv(), so a tab inside a field would silently create a fourth column and mis-position every
    downstream input.
  * group_index is emitted as a small integer string and is consumed as a String by the workflow,
    because read_tsv hands back Strings and WDL will not coerce String->Int at a call site.

Only the FIRST data record of each shard is read, and iteration (not a tabix query) is used, so no
index file is required and nothing is localized whole-file-wise.
"""
import sys

import pysam

# Order chromosomes so gather is deterministic and chromosome-ordered rather than discovery-ordered.
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
SEP = "|"


def chrom_key(contig):
    try:
        return (0, CHROM_ORDER.index(contig))
    except ValueError:
        return (1, CHROM_ORDER.__len__())  # unknown contigs last, stably


def main(argv):
    if len(argv) != 3:
        sys.stderr.write(__doc__)
        return 2
    shards_file, max_per = argv[1], int(argv[2])
    if max_per < 1:
        sys.stderr.write("max_shards_per_group must be >= 1\n")
        return 2

    groups, empty, failed = {}, [], []
    with open(shards_file) as fh:
        uris = [ln.strip() for ln in fh if ln.strip()]
    if not uris:
        sys.stderr.write("no shard URIs supplied\n")
        return 3

    for u in uris:
        try:
            vf = pysam.VariantFile(u)
        except Exception as exc:  # noqa: BLE001 - report which URI and why; do not guess a contig
            failed.append((u, f"{type(exc).__name__}: {exc}"))
            continue
        contig = None
        for rec in vf:
            contig = rec.chrom
            break
        if contig is None:
            empty.append(u)
            continue
        groups.setdefault(contig, []).append(u)

    if failed:
        for u, why in failed:
            sys.stderr.write(f"UNREADABLE shard {u}: {why}\n")
        sys.stderr.write(f"{len(failed)} of {len(uris)} shards could not be opened; refusing to group\n")
        return 4
    for u in empty:
        sys.stderr.write(f"WARN empty shard skipped: {u}\n")
    if not groups:
        sys.stderr.write("every shard was empty; nothing to annotate\n")
        return 5

    for contig in sorted(groups, key=chrom_key):
        uris_g = groups[contig]
        for i in range(0, len(uris_g), max_per):
            chunk = uris_g[i:i + max_per]
            print(f"{contig}\t{i // max_per}\t{SEP.join(chunk)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
