version 1.0

## Fix invalid read-group PL tag in a BAM or CRAM.
##
## Picard ValidateSamFile rejects PL:"NovaSeq X" because PL must come from the
## SAM spec controlled vocabulary (ILLUMINA, PACBIO, ONT, ...). "NovaSeq X" is
## an instrument *model*, which belongs in PM. This workflow rewrites every
## @RG so that:
##   PL:NovaSeq X  ->  PL:ILLUMINA   (and adds PM:NovaSeq X if no PM present)
## It is a header-only rewrite (samtools reheader) -- no realignment, no
## reference required, runs in seconds-to-minutes depending on file size.
##
## Input format (BAM vs CRAM) is detected by content (magic bytes), not the
## filename, so a mislabeled extension is handled correctly. The output and its
## index are written with the matching extension (.bam/.bai or .cram/.crai).
##
## Run once per affected sample (3 of them); sample_name identifies the sample
## and names the outputs.

workflow FixReadGroupPlatform {
  input {
    File   input_reads              # BAM or CRAM (format auto-detected)
    String sample_name              # SM, used for output naming / traceability
    String bad_pl    = "NovaSeq X"  # the offending PL value to replace
    String good_pl   = "ILLUMINA"   # SAM-spec platform to set
    String model_pm  = "NovaSeq X"  # value to preserve in PM
    String samtools_docker = "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    Int    additional_disk = 20
  }

  call ReheaderReads {
    input:
      input_reads     = input_reads,
      sample_name     = sample_name,
      bad_pl          = bad_pl,
      good_pl         = good_pl,
      model_pm        = model_pm,
      samtools_docker = samtools_docker,
      additional_disk = additional_disk
  }

  output {
    File reheadered_reads = ReheaderReads.reheadered_reads
    File reheadered_index = ReheaderReads.reheadered_index
  }
}

task ReheaderReads {
  input {
    File   input_reads
    String sample_name
    String bad_pl
    String good_pl
    String model_pm
    String samtools_docker
    Int    additional_disk
  }

  # reheader copies the whole file -> need room for input + output.
  Int disk_gb = ceil(size(input_reads, "GB") * 2) + additional_disk

  command <<<
    set -euo pipefail

    # Detect format by content, not extension. CRAM files begin with the
    # ASCII bytes "CRAM"; BAM files are bgzip (0x1f 0x8b ...). Anything that
    # is not CRAM is treated as BAM.
    if [ "$(head -c 4 "~{input_reads}")" = "CRAM" ]; then
      ext=cram; idx=crai
    else
      ext=bam;  idx=bai
    fi
    echo "Detected format: ${ext}" >&2

    samtools quickcheck "~{input_reads}" || { echo "ERROR: input failed samtools quickcheck" >&2; exit 1; }

    samtools view -H "~{input_reads}" > header.sam

    # Tab-delimited @RG fields. Replace the exact PL value; append PM:<model>
    # once per @RG line that lacks a PM tag. Non-@RG lines pass through.
    awk -v bad="PL:~{bad_pl}" -v good="PL:~{good_pl}" -v pm="PM:~{model_pm}" '
      BEGIN { FS = OFS = "\t" }
      /^@RG/ {
        haspm = 0
        for (i = 1; i <= NF; i++) if ($i ~ /^PM:/) haspm = 1
        for (i = 1; i <= NF; i++) if ($i == bad) $i = good
        if (!haspm) $0 = $0 OFS pm
        print; next
      }
      { print }
    ' header.sam > header.fixed.sam

    # Fail loudly if the bad value survived anywhere in the header.
    if grep -qF "PL:~{bad_pl}" header.fixed.sam; then
      echo "ERROR: PL:~{bad_pl} still present after fix" >&2
      exit 1
    fi

    out="~{sample_name}.reheadered.${ext}"
    samtools reheader header.fixed.sam "~{input_reads}" > "${out}"
    samtools index "${out}"           # writes ${out}.${idx}
  >>>

  runtime {
    docker: samtools_docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }

  output {
    # Exactly one of each pair exists per run; flatten + index 0 picks it.
    File reheadered_reads = flatten([glob("~{sample_name}.reheadered.bam"),
                                     glob("~{sample_name}.reheadered.cram")])[0]
    File reheadered_index = flatten([glob("~{sample_name}.reheadered.bam.bai"),
                                     glob("~{sample_name}.reheadered.cram.crai")])[0]
  }
}
