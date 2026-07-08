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
    Int    additional_disk_gb = 200
    Int    mem_gb = 7
    Int    compression_level = 2
    Int    preemptible = 0
  }

  call SortQueryname {
    input:
      input_bam          = input_bam,
      output_basename    = output_basename,
      docker             = picard_docker,
      additional_disk_gb = additional_disk_gb,
      mem_gb             = mem_gb,
      compression_level  = compression_level,
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
    Int    compression_level
    Int    preemptible
  }

  # SortSam spills to disk heavily (MAX_RECORDS_IN_RAM=300000, uncompressed spill),
  # so size the local disk for input + output + spill. Same 3.25x multiplier WARP
  # uses for its SortSam task (tasks/broad/BamProcessing.wdl).
  Int disk_gb        = ceil(size(input_bam, "GB") * 3.25 + additional_disk_gb)
  Int machine_mem_mb = mem_gb * 1024
  Int java_xmx_mb    = machine_mem_mb - 512
  Int java_xms_mb    = machine_mem_mb - 1024

  command <<<
    set -euo pipefail

    # Spill to the mounted local-disk (cwd), NOT /tmp: java.io.tmpdir defaults to
    # /tmp, which on Terra/Cromwell is the small boot disk. Picard's TMP_DIR falls
    # back to java.io.tmpdir, so set both. Use $PWD (resolves to the execution dir
    # on the local-disk regardless of PAPIv2 /cromwell_root vs Batch mount path).
    mkdir -p "${PWD}/tmp"

    java -Dsamjdk.compression_level=~{compression_level} \
      -Djava.io.tmpdir="${PWD}/tmp" \
      -Xms~{java_xms_mb}m -Xmx~{java_xmx_mb}m \
      -jar /usr/picard/picard.jar SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_basename}.qname_sorted.unmapped.bam \
      SORT_ORDER=queryname \
      MAX_RECORDS_IN_RAM=300000 \
      TMP_DIR="${PWD}/tmp" \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
  >>>

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }

  output {
    File unmapped_bam = "~{output_basename}.qname_sorted.unmapped.bam"
  }
}
