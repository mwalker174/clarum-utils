#!/usr/bin/env python3
"""Render one task's `command <<< >>>` block from a WDL with Cromwell-style localized
paths, so the bash can be replayed outside Cromwell.

Not a WDL engine: it substitutes only the placeholders the task command uses, and refuses
to guess if it meets one it does not know. The point is to reproduce the *cloud* layout -
inputs materialized at <call_root>/<bucket>/<prefix>/<name> while the task runs in
<call_root> - because a plain local run of this WDL hides exactly that difference (see
test/replay_cromwell_layout.sh and docs/progress/072 of the clarum workspace).
"""
import argparse
import re

PLACEHOLDER = re.compile(r"~\{([^{}]*)\}")


def command_block(wdl_text, task):
    m = re.search(
        r"task %s \{.*?command <<<\n(.*?)\n\s*>>>" % re.escape(task), wdl_text, re.S
    )
    if not m:
        raise SystemExit("task %s not found, or has no command <<< >>> block" % task)
    return m.group(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wdl")
    ap.add_argument("out")
    ap.add_argument("--task", default="ConvertBamToCram")
    ap.add_argument("--root", required=True, help="emulated Cromwell call root")
    ap.add_argument("--ref-md5", required=True)
    args = ap.parse_args()

    cr = args.root.rstrip("/")
    values = {
        "input_bam": "%s/datarepo-fake/uuid1/A.bam" % cr,
        "output_basename": "A",
        "cpu": "2",
        'default="" input_bam_index': "%s/datarepo-fake/uuid2/A.bam.bai" % cr,
        'default="" input_bam_md5': "%s/datarepo-fake/uuid3/A.bam.md5sum" % cr,
        "ref_fasta": "%s/gcp-public-data--broad-references/hg38/v0/mini.fasta" % cr,
        "ref_fasta_index": "%s/gcp-public-data--broad-references/hg38/v0/mini.fasta.fai" % cr,
        "ref_fasta_md5": args.ref_md5,
        'default="" reference_m5_table': "",
        "run_roundtrip": "true",
        "strict_roundtrip": "true",
        "roundtrip_window_bp": "500",
        "roundtrip_windows_per_contig": "2",
        "roundtrip_contigs": "12",
        # Revisions before the fix addressed the optional inputs through workflow-level
        # Strings, which hold the raw cloud URI (that is the bug, not the harness). Mapping
        # them keeps historical revisions replayable: git show HEAD~1:bam_to_cram.wdl.
        "bam_index_path": "gs://datarepo-fake/uuid2/A.bam.bai",
        "bam_md5_path": "gs://datarepo-fake/uuid3/A.bam.md5sum",
    }

    text = open(args.wdl).read()
    body = command_block(text, args.task)

    def sub(mo):
        key = mo.group(1).strip()
        if key not in values:
            raise SystemExit("unmapped placeholder in %s command: %r" % (args.task, key))
        return values[key]

    open(args.out, "w").write(PLACEHOLDER.sub(sub, body))


if __name__ == "__main__":
    main()
