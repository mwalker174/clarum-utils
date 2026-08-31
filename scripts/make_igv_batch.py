#!/usr/bin/env python3
"""Generate an IGV 2.4.14 batch script + snapshot manifest for manual review
of selected Mutect2 mosaic candidate calls (point SNV/indel), one CRAM per
sample.

Point-variant analog of talkowski-lab/gatk-sv-internal's
dockerfiles/igv/makeigvsplit_cram.py (which does the SV/breakpoint-pair
version of the same idea, called by wdl/IGVGeneratePlotsWholeGenome.wdl).
Batch-command reference: https://github.com/igvteam/igv/wiki/Batch-commands

Input is any TSV with sample/chrom/pos/ref/alt columns (e.g.
tables/igv_review_sample.tsv, or a hand-curated subset of it — this script
doesn't care where the rows came from, only that they're for one sample).
No pandas dependency, so it runs unmodified inside a bare IGV/Xvfb container.

Usage (local dry run, no Terra):
    python make_igv_batch.py --tsv tables/igv_review_sample.tsv --sample 007 \
        --cram data/007.cram --crai data/007.cram.crai \
        --fasta ref/Homo_sapiens_assembly38.fasta --out-dir igv_batch/007
    bash IGV_2.4.14/igv.sh -b igv_batch/007/igv_batch.txt

Companion WDL: IGVSnapshotMutectVariants.wdl (per-sample scatter).
"""
import argparse
import csv
import os


def load_rows(tsv_path, sample):
    rows = []
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["sample"] == sample:
                rows.append(row)
    return rows


def build_batch(rows, sample, cram, fasta, out_dir, buff=75, genome_id=None, max_panel_height=1200,
                pon_vcf=None, bamout=None):
    screenshot_dir = os.path.join(out_dir, "screenshots")
    os.makedirs(screenshot_dir, exist_ok=True)
    batch_path = os.path.join(out_dir, "igv_batch.txt")
    manifest_path = os.path.join(out_dir, "snapshot_manifest.tsv")
    manifest_rows = []

    with open(batch_path, "w") as f:
        f.write("new\n")
        f.write(f"genome {genome_id or fasta}\n")
        f.write(f"snapshotDirectory {screenshot_dir}\n")
        f.write(f"maxPanelHeight {max_panel_height}\n")
        # Off by default in IGV 2.4.14 -- without this, insertion/deletion
        # support that the aligner soft-clipped instead of representing as an
        # indel op is invisible in the screenshot even when reads carry it
        # (see IGV_REVIEW_SOP.md).
        f.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
        # Draws a vertical line at the exact screenshot-center coordinate so a
        # reviewer can tell which column is the variant without counting bp
        # ticks -- goto centers the locus pixel-exact, but nothing marked it.
        f.write("preference SAM.SHOW_CENTER_LINE true\n")
        f.write(f"load {cram}\n")
        if pon_vcf:
            f.write(f"load {pon_vcf}\n")  # PoN-recurrence check, IGV_REVIEW_SOP.md
        for r in rows:
            chrom, pos, ref, alt = r["chrom"], int(r["pos"]), r["ref"], r["alt"]
            start, end = max(1, pos - buff), pos + buff
            locus = f"{chrom}:{start}-{end}"
            png = f"{sample}__{chrom}-{pos}-{ref}-{alt}.png"
            f.write(f"goto {locus}\n")
            f.write(f"sort base {chrom}:{pos}\n")
            # "colorBy" isn't a real batch command in IGV 2.4.14 (added in 2.19.2) --
            # silently logged as an unknown command and ignored. "group strand" IS
            # supported here (AlignmentTrack.GroupOption.STRAND) and gives the same
            # visual fwd/rev separation the strand-balance check needs.
            f.write("group strand\n")
            f.write("viewaspairs\n")
            f.write("squish\n")
            f.write(f"snapshot {png}\n")
            row = {**r, "locus": locus, "png": png}
            if bamout:
                row["png_bamout"] = f"{sample}__{chrom}-{pos}-{ref}-{alt}__bamout.png"
            manifest_rows.append(row)

        # Second pass over the same sites with Mutect2's bamout loaded.
        #
        # IGV shows BWA's alignment, but Mutect2 doesn't call from that -- it
        # reassembles the active region into haplotypes and realigns reads to
        # them with PairHMM, and AD/TLOD come from that realignment. So a site
        # that looks unreadable here can still be a confident call, and the
        # only way to see the caller's own view is bamout (GATK --bam-output).
        # Reads carry an HC tag naming the haplotype they were assigned to;
        # grouping on it splits the pileup the way the caller split it.
        #
        # Separate screenshots rather than one combined image because IGV's
        # "group" applies to every alignment track at once -- there is no
        # per-track form -- so the CRAM cannot stay grouped by strand while
        # the bamout is grouped by HC.
        if bamout:
            f.write(f"load {bamout}\n")
            # GroupOption.TAG exists in 2.4.14 and getAlignmentGroupOption
            # falls through to valueOf(), so "group TAG HC" resolves even
            # though only strand/sample/library/read_group are named
            # explicitly in that method.
            f.write("group TAG HC\n")
            for r in rows:
                chrom, pos, ref, alt = r["chrom"], int(r["pos"]), r["ref"], r["alt"]
                start, end = max(1, pos - buff), pos + buff
                locus = f"{chrom}:{start}-{end}"
                png = f"{sample}__{chrom}-{pos}-{ref}-{alt}__bamout.png"
                f.write(f"goto {locus}\n")
                f.write(f"sort base {chrom}:{pos}\n")
                # No viewaspairs: the artificial haplotype records bamout adds
                # aren't paired, and pairing the real reads hides the
                # read-to-haplotype assignment this pass exists to show.
                f.write("squish\n")
                f.write(f"snapshot {png}\n")

        f.write("exit\n")

    with open(manifest_path, "w", newline="") as f:
        fieldnames = list(manifest_rows[0].keys()) if manifest_rows else ["sample", "chrom", "pos", "ref", "alt", "locus", "png"]
        if bamout and "png_bamout" not in fieldnames:
            fieldnames.append("png_bamout")
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    return batch_path, manifest_path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv", required=True, help="variant sheet with sample/chrom/pos/ref/alt columns")
    p.add_argument("--sample", required=True, help="restrict to this sample column value")
    p.add_argument("--cram", required=True)
    p.add_argument("--crai", help="unused directly by IGV (must sit alongside --cram); kept so WDL localizes it")
    p.add_argument("--fasta", required=True, help="reference fasta — must match the CRAM's contig naming (chr-prefixed b38)")
    p.add_argument("--fasta-idx", help="unused directly; kept so WDL localizes it")
    p.add_argument("--pon-vcf", default=None, help="cross-cohort PoN VCF, loaded as a second track for PoN-recurrence check")
    p.add_argument("--pon-vcf-idx", default=None, help="unused directly (Tribble .idx alongside --pon-vcf); kept so WDL localizes it")
    p.add_argument("--bamout", help="Mutect2 --bam-output BAM (reads realigned to assembled haplotypes); "
                                    "adds a second screenshot per site grouped by the HC haplotype tag")
    p.add_argument("--bamout-index", help="unused directly by IGV (must sit alongside --bamout); kept so the "
                                          "WDL localizes it")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--buff", type=int, default=75, help="flanking bp each side of the variant")
    p.add_argument("--max-panel-height", type=int, default=1200)
    p.add_argument("--genome-id", default=None, help="IGV registered genome ID to use instead of --fasta (e.g. hg38)")
    a = p.parse_args()

    rows = load_rows(a.tsv, a.sample)
    if not rows:
        raise SystemExit(f"no rows for sample {a.sample!r} in {a.tsv}")
    os.makedirs(a.out_dir, exist_ok=True)
    batch_path, manifest_path = build_batch(
        rows, a.sample, a.cram, a.fasta, a.out_dir,
        buff=a.buff, genome_id=a.genome_id, max_panel_height=a.max_panel_height, pon_vcf=a.pon_vcf,
                                        bamout=a.bamout)
    print(f"[igv_batch] {len(rows)} variants for {a.sample} -> {batch_path}")
    print(f"[igv_batch] snapshot manifest -> {manifest_path}")


if __name__ == "__main__":
    main()
