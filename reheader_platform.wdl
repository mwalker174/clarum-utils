version 1.0

## Fix invalid read-group PL tag in a CRAM.
##
## Picard ValidateSamFile rejects PL:"NovaSeq X" because PL must come from the
## SAM spec controlled vocabulary (ILLUMINA, PACBIO, ONT, ...). "NovaSeq X" is
## an instrument *model*, which belongs in PM. This workflow rewrites every
## @RG so that:
##   PL:NovaSeq X  ->  PL:ILLUMINA   (and adds PM:NovaSeq X if no PM present)
## It is a header-only rewrite (samtools reheader) -- no realignment, no
## reference required, runs in seconds-to-minutes depending on CRAM size.
##
## Run once per affected sample (3 of them); sample_name identifies the sample
## and names the outputs.

workflow FixReadGroupPlatform {
  input {
    File   input_cram
    String sample_name              # SM, used for output naming / traceability
    String bad_pl    = "NovaSeq X"  # the offending PL value to replace
    String good_pl   = "ILLUMINA"   # SAM-spec platform to set
    String model_pm  = "NovaSeq X"  # value to preserve in PM
    String samtools_docker = "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    Int    additional_disk = 20
  }

  call ReheaderCram {
    input:
      input_cram      = input_cram,
      sample_name     = sample_name,
      bad_pl          = bad_pl,
      good_pl         = good_pl,
      model_pm        = model_pm,
      samtools_docker = samtools_docker,
      additional_disk = additional_disk
  }

  output {
    File reheadered_cram      = ReheaderCram.reheadered_cram
    File reheadered_cram_index = ReheaderCram.reheadered_cram_index
  }
}

task ReheaderCram {
  input {
    File   input_cram
    String sample_name
    String bad_pl
    String good_pl
    String model_pm
    String samtools_docker
    Int    additional_disk
  }

  # reheader copies the whole CRAM -> need room for input + output.
  Int disk_gb = ceil(size(input_cram, "GB") * 2) + additional_disk

  command <<<
    set -euo pipefail

    samtools view -H "~{input_cram}" > header.sam

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

    samtools reheader header.fixed.sam "~{input_cram}" > "~{sample_name}.reheadered.cram"
    samtools index "~{sample_name}.reheadered.cram"
  >>>

  runtime {
    docker: samtools_docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }

  output {
    File reheadered_cram       = "~{sample_name}.reheadered.cram"
    File reheadered_cram_index = "~{sample_name}.reheadered.cram.crai"
  }
}
