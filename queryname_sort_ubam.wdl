version 1.0

## Queryname-sort a flowcell "unmapped" BAM so WholeGenomeGermlineSingleSample's
## MergeBamAlignment succeeds.
##
## The wilkinshaug stillbirth uBAMs were built CRAM -> CramToFastq ->
## PairedFastQsToUnmappedBAM (Picard FastqToSam). FastqToSam preserves the input
## read order (here: the source CRAM's coordinate order, across 4 read groups /
## flowcells) and tags the header `SO:unsorted`. WARP trusts the header and never
## sorts, so SamToFastqAndBwaMemAndMba feeds MergeBamAlignment an iterator that is
## not in Picard lexicographic queryname order. The first merge attempt throws
##   "Underlying iterator is not queryname sorted: H77L7...:2/2 > H77J7...:1/2"
## Picard then retries with forceSort, but ALIGNED_BAM=/dev/stdin (the bwa pipe)
## is already drained, so the aligned dictionary reads back empty and the visible
## failure becomes the misleading
##   "Do not use this function to merge dictionaries... Found [] and [chr1, ...]"
## The [] is the drained aligned stream, NOT the uBAM's @SQ. Stripping @SQ (the
## earlier strip_sequence_dictionary.wdl fix) therefore does nothing.
##
## Real fix: sort the uBAM by queryname with Picard (matching the SO:queryname
## order that CramToUnmappedBams produces for standard WARP inputs). Use Picard,
## not `samtools sort -n`: MergeBamAlignment requires Picard's lexicographic
## queryname order, which samtools' natural name sort does not reproduce.
##
## Output: <basename>.qname_sorted.unmapped.bam (SO:queryname) -> wire into
##   WholeGenomeGermlineSingleSample.sample_and_unmapped_bams.flowcell_unmapped_bams

workflow QuerynameSortUbam {
  input {
    File   input_bam                       # an UNMAPPED BAM that is not queryname-sorted
    String output_basename = basename(input_bam, ".bam")
    String picard_docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    Int    additional_disk_gb = 50
    Int    mem_gb = 16
    Int    preemptible = 0
  }

  call SortQueryname {
    input:
      input_bam          = input_bam,
      output_basename    = output_basename,
      docker             = picard_docker,
      additional_disk_gb = additional_disk_gb,
      mem_gb             = mem_gb,
      preemptible        = preemptible
  }

  output {
    File unmapped_bam = SortQueryname.unmapped_bam
  }

  meta { description: "Picard queryname-sort an unmapped BAM so WGSingleSample's MergeBamAlignment iterator check passes." }
}

task SortQueryname {
  input {
    File   input_bam
    String output_basename
    String docker
    Int    additional_disk_gb
    Int    mem_gb
    Int    preemptible
  }

  # input + output + Picard spill-to-disk temp (~1x records each)
  Int disk_gb = ceil(size(input_bam, "GB") * 3.25 + additional_disk_gb)
  Int jvm_mem_mb = (mem_gb - 2) * 1024

  command <<<
    set -euo pipefail

    java -Xms~{jvm_mem_mb}m -Xmx~{jvm_mem_mb}m -jar /usr/gitc/picard.jar SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_basename}.qname_sorted.unmapped.bam \
      SORT_ORDER=queryname \
      MAX_RECORDS_IN_RAM=3000000 \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
  >>>

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File unmapped_bam = "~{output_basename}.qname_sorted.unmapped.bam"
  }
}
