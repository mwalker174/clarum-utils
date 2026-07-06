version 1.0

## Strip the @SQ sequence dictionary from an unmapped BAM header.
##
## The wilkinshaug flowcell "unmapped" BAMs carry a full @SQ sequence
## dictionary (3366 contigs) even though every read is unmapped. A proper uBAM
## has no @SQ. That stale dictionary makes WholeGenomeGermlineSingleSample fail
## in SamToFastqAndBwaMemAndMba at Picard MergeBamAlignment:
##   "Do not use this function to merge dictionaries with different sequences...
##    Found [] and [chr1, chr2, ...]"
##
## Fix is header-only: drop the @SQ lines with `samtools reheader`. The reads
## are already unmapped (RNAME=*), so removing the dictionary is safe and lets
## MergeBamAlignment take its dictionary from the reference. No realignment, no
## revert, no reference required -- runs in seconds-to-minutes.
##
## Output: <basename>.nodict.unmapped.bam -> wire into
##   WholeGenomeGermlineSingleSample.sample_and_unmapped_bams.flowcell_unmapped_bams

workflow StripSequenceDictionary {
  input {
    File   input_bam                       # an UNMAPPED BAM whose header carries a stale @SQ
    String output_basename = basename(input_bam, ".bam")
    String samtools_docker = "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    Int    additional_disk_gb = 20
    Int    mem_gb = 4
    Int    preemptible = 1
  }

  call StripSQ {
    input:
      input_bam        = input_bam,
      output_basename  = output_basename,
      docker           = samtools_docker,
      additional_disk_gb = additional_disk_gb,
      mem_gb           = mem_gb,
      preemptible      = preemptible
  }

  output {
    File unmapped_bam = StripSQ.unmapped_bam
  }

  meta { description: "Remove the @SQ sequence dictionary from an unmapped BAM header so WGSingleSample's MergeBamAlignment succeeds." }
}

task StripSQ {
  input {
    File   input_bam
    String output_basename
    String docker
    Int    additional_disk_gb
    Int    mem_gb
    Int    preemptible
  }

  Int disk_gb = ceil(size(input_bam, "GB") * 2 + additional_disk_gb)

  command <<<
    set -euo pipefail

    # Safety: this is a header-only fix and is only valid when the reads are
    # unmapped. Abort if any read is mapped (would leave dangling RNAME refs).
    if [ "$(samtools view -c -F 0x4 ~{input_bam})" -ne 0 ]; then
      echo "ERROR: input has mapped reads; strip-@SQ is unsafe. Revert alignment first." >&2
      exit 1
    fi

    # Drop every @SQ line; keep @HD/@RG/@PG/@CO.
    samtools view -H ~{input_bam} | grep -v '^@SQ' > header.noSQ.sam

    samtools reheader header.noSQ.sam ~{input_bam} > ~{output_basename}.nodict.unmapped.bam

    # Verify: no @SQ remains and the file is intact.
    if samtools view -H ~{output_basename}.nodict.unmapped.bam | grep -q '^@SQ'; then
      echo "ERROR: @SQ still present after reheader" >&2; exit 1
    fi
    samtools quickcheck ~{output_basename}.nodict.unmapped.bam
  >>>

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File unmapped_bam = "~{output_basename}.nodict.unmapped.bam"
  }
}
