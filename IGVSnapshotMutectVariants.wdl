version 1.0

## IGVSnapshotMutectVariants
##
## Batch IGV screenshots of selected Mutect2 mosaic candidate calls (point
## SNV/indel), one CRAM per sample, for manual review — CLARUM's MOC10
## mosaic-detection QC pipeline's IGV read-evidence step is fully manual,
## this just automates taking the screenshots.
##
## Point-variant analog of talkowski-lab/gatk-sv-internal's
## wdl/IGVGeneratePlotsWholeGenome.wdl + dockerfiles/igv/makeigvsplit_cram.py,
## which do the SV/breakpoint-pair version of the same idea. Batch-command
## reference: https://github.com/igvteam/igv/wiki/Batch-commands
## Runs IGV 2.4.14: http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip
##
## variant_tsv: any TSV with sample/chrom/pos/ref/alt columns (e.g.
## clarum's tables/igv_review_sample.tsv, or a hand-curated subset — CHIP
## hits, high-VAF-band candidates, whatever the reviewer wants a look at).
## samples/crams/cram_indexes are parallel arrays (same order, one per
## sample); each sample's CRAM is loaded once and every one of its rows in
## variant_tsv gets its own screenshot.
##
## igv_batch_script: scripts/make_igv_batch.py in this repo, staged wherever
## Terra can localize it (e.g. uploaded to the workspace bucket) — kept as a
## single checked-in source of truth rather than duplicated inline here.
##
## Output per sample: a tarball of screenshots + a snapshot_manifest.tsv
## (variant row -> expected .png filename) to reconcile against after a
## human looks at the pictures and fills in igv_verdict/notes.

workflow IGVSnapshotMutectVariants {
  input {
    File variant_tsv
    File igv_batch_script
    Array[String] samples
    Array[File] crams
    Array[File] cram_indexes
    File ref_fasta
    File ref_fasta_index
    Int buff = 25
    Int max_panel_height = 1200
    String docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    Int preemptible = 1
  }

  scatter (i in range(length(samples))) {
    call MakeSnapshots {
      input:
        variant_tsv = variant_tsv,
        igv_batch_script = igv_batch_script,
        sample = samples[i],
        cram = crams[i],
        cram_index = cram_indexes[i],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        buff = buff,
        max_panel_height = max_panel_height,
        docker = docker,
        preemptible = preemptible
    }
  }

  output {
    Array[File] screenshot_tarballs = MakeSnapshots.tar_gz
    Array[File] snapshot_manifests = MakeSnapshots.snapshot_manifest
  }

  meta { description: "Per-sample IGV batch screenshots of selected Mutect2 mosaic candidate calls, for manual review." }
}

task MakeSnapshots {
  input {
    File variant_tsv
    File igv_batch_script
    String sample
    File cram
    File cram_index
    File ref_fasta
    File ref_fasta_index
    Int buff
    Int max_panel_height
    String docker
    Int preemptible
  }

  Int disk_gb = ceil(size(cram, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    apt-get update -qq
    apt-get install -y -qq --no-install-recommends xvfb openjdk-11-jre-headless wget unzip python3 > /dev/null

    wget -q http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip -O igv.zip
    unzip -q igv.zip

    python3 ~{igv_batch_script} \
      --tsv ~{variant_tsv} --sample ~{sample} \
      --cram ~{cram} --crai ~{cram_index} \
      --fasta ~{ref_fasta} --fasta-idx ~{ref_fasta_index} \
      --buff ~{buff} --max-panel-height ~{max_panel_height} \
      --out-dir out

    xvfb-run --auto-servernum --server-args="-screen 0 1920x1080x24" \
      bash IGV_2.4.14/igv.sh -b out/igv_batch.txt

    tar -czf ~{sample}.igv_screenshots.tar.gz -C out screenshots
    cp out/snapshot_manifest.tsv ~{sample}.snapshot_manifest.tsv
  >>>

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File tar_gz = "~{sample}.igv_screenshots.tar.gz"
    File snapshot_manifest = "~{sample}.snapshot_manifest.tsv"
  }
}
