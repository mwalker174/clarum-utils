version 1.0

## Convert an aligned BAM to a CRAM. Deliverables: .cram, .crai, and an md5 of each.
##
## TDD 2.H.3/2.H.4: the prospective somatic ES deliverable requires CRAM generation
## with integrity by file size and MD5. The BI TDR snapshots (81 rows, both
## prospective-WES workspaces) ship BAM only - this is the converter.
##
## The reference is the whole ballgame: CRAM pins the reference by MD5, and the
## delivered DRAGEN BAMs do not say which hg38 they used. Settled in
## docs/progress/070 (BI GS-CRAM M5 table, genome-wide NM audit, full round trip):
##   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
## (file md5 7ff134953dcca8c8997453bbb80b6b5e). The ONLY gate left in the task is
## that md5: a wrong FASTA burns a 3 h conversion to make a CRAM nobody can decode.
## Everything else the task used to do (source-BAM md5, record-count and
## sampled-window round-trip, @SQ M5 join) is redundant evidence - the record count
## and MD5 of the source are TDR's integrity claim, and the round trip was proved
## offline on 616k real records. The CRAM's own @SQ M5 block is checked offline
## against data/qc-inputs/bi_gs_cram_sq_m5.tsv after the fact.
##
## Output contract: the four Files below, deliberately nothing else. All File? so a
## late failure cannot strand the CRAM behind a missing required output
## (submission 0c6ccba5); no derived outputs because every extra expression is
## another thing to evaluate after a two-hour call (cf4cb95f died on size(x,
## "Bytes") - "Bytes" is not a unit Cromwell knows, and miniwdl accepts it).
##
## Cromwell traps this file still obeys (miniwdl check passes all of them):
## - meta { } entries are newline-separated, no commas;
## - File inputs pass into the task as File (casting to String in the workflow
##   yields the gs:// URI and the object is never localized);
## - the command uses the rendered localized paths, never basename();
## - size()/memory units come from B/KB/MB/GB/TB/PB + KiB/MiB/GiB/TiB/PiB.

workflow BamToCram {
  input {
    File    input_bam
    # Keeps the collaborator sample id, drops any .bam suffix.
    String  output_basename = sub(basename(input_bam), "\\.bam$", "")

    # The reference DRAGEN used. Required, deliberately un-defaulted - see the header.
    File    ref_fasta
    File    ref_fasta_index
    # md5 of ref_fasta (broad hg38/v0 FASTA = 7ff134953dcca8c8997453bbb80b6b5e).
    # Set it: it turns "we pointed at a hg38" into "we pointed at THE hg38 we audited".
    String  ref_fasta_md5 = ""

    # samtools 1.11 image WARP's ConvertToCram runs on.
    String  samtools_docker = "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"

    Int     cpu = 8                        # samtools -@ threads
    Int     mem_gb = 12
    Int     preemptible = 3                # a restart re-streams the BAM; raise to 0 if flaky
  }

  call ConvertBamToCram {
    input:
      input_bam        = input_bam,
      output_basename  = output_basename,
      ref_fasta        = ref_fasta,
      ref_fasta_index  = ref_fasta_index,
      ref_fasta_md5    = ref_fasta_md5,
      samtools_docker  = samtools_docker,
      cpu              = cpu,
      mem_gb           = mem_gb,
      preemptible      = preemptible
  }

  output {
    File? output_cram           = ConvertBamToCram.output_cram
    File? output_cram_index     = ConvertBamToCram.output_cram_index
    File? output_cram_md5       = ConvertBamToCram.output_cram_md5
    File? output_cram_index_md5 = ConvertBamToCram.output_cram_index_md5
  }

  ## meta { } is newline-separated, no commas - Terra's Cromwell parser rejects commas
  ## ("Expected rbrace, got ','") while miniwdl accepts them.
  meta {
    description: "BAM -> CRAM for the prospective somatic ES delivery (TDD 2.H.3 / 2.H.4): convert, index, md5. Reference identity gated on file md5."
    summary: "Convert an aligned BAM to CRAM with index and md5s"
    author: "CLARUM / Talkowski lab"
  }
}

task ConvertBamToCram {
  input {
    File    input_bam
    String  output_basename
    File    ref_fasta
    File    ref_fasta_index
    String  ref_fasta_md5
    String  samtools_docker
    Int     cpu
    Int     mem_gb
    Int     preemptible
  }

  # Localized BAM (1x) + CRAM (~0.4x) + reference + headroom.
  Int disk_gb        = ceil(size(input_bam, "GB") + size(ref_fasta, "GB") + 30)
  Int machine_mem_mb = mem_gb * 1024

  command <<<
    set -euo pipefail

    samtools --version | sed -n '1,2p'

    REF='~{ref_fasta}'
    # The .fai must sit beside the FASTA (htslib looks there before building one,
    # and building it costs minutes per call).
    [ -s "${REF}.fai" ] || ln -sf '~{ref_fasta_index}' "${REF}.fai"

    # The one gate: htslib only checks reference md5 when it decodes a slice, so a
    # wrong FASTA can burn the whole conversion before it says anything. Check up
    # front; it is ~20 s of md5sum on a 3.2 GB file.
    REF_MD5=$(md5sum "$REF" | awk '{print $1}')
    echo "reference md5=${REF_MD5}"
    if [ -n '~{ref_fasta_md5}' ] && [ "${REF_MD5}" != '~{ref_fasta_md5}' ]; then
      echo "FATAL: reference md5 ${REF_MD5} != expected ~{ref_fasta_md5}" >&2
      exit 1
    fi

    OUT="~{output_basename}"
    # version=3.0 explicitly: the delivered GS CRAMs are 3.0; 4.0 would be a
    # unilateral change. -T gives the compressor the reference directly.
    samtools view -C -@ ~{cpu} -O cram,version=3.0 -T "$REF" -o "${OUT}.cram" '~{input_bam}'

    # Indexing a CRAM resolves the reference through the path embedded at encode
    # time, which still exists in this container - no REF_CACHE dance needed.
    samtools index -@ ~{cpu} "${OUT}.cram"

    md5sum "${OUT}.cram" | awk '{print $1}' > "${OUT}.cram.md5"
    md5sum "${OUT}.cram.crai" | awk '{print $1}' > "${OUT}.cram.crai.md5"
    echo "cram_md5=$(cat "${OUT}.cram.md5")  crai_md5=$(cat "${OUT}.cram.crai.md5")"
  >>>

  runtime {
    docker: samtools_docker
    memory: machine_mem_mb + " MB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  # The deliverable, and nothing else: four optional Files, completeness implied by
  # set -e (an rc of 0 means every line above, including both md5 writes, succeeded).
  output {
    File? output_cram           = "~{output_basename}.cram"
    File? output_cram_index     = "~{output_basename}.cram.crai"
    File? output_cram_md5       = "~{output_basename}.cram.md5"
    File? output_cram_index_md5 = "~{output_basename}.cram.crai.md5"
  }

  meta {
    description: "samtools view -C + index + md5; reference identity gated on file md5."
  }
}
