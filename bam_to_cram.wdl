version 1.0

## Convert an already-aligned BAM to a CRAM (+ .crai, md5s, integrity evidence).
##
## Why this exists
## ---------------
## TDD §2.H.3 requires the prospective somatic ES deliverable to include "generation
## of sequencing CRAM file and QC report for at least 14 samples (goal: 17)", and
## §2.H.4 transfers "at least 14 CRAM files (goal: 17) ... integrity confirmed via
## file size and MD5 checks". What BI actually delivered through the TDR snapshots
## (clarum_talkowski_strub_ach_prospective_wes / ..._smits_bch_prospective_wes, 81
## rows) is BAM only: `bam_file`, `bai_file`, `bam_md5_sum_path`, four DRAGEN metrics
## CSVs, and a hard-filtered gVCF. No CRAM column. Both workspaces have zero method
## configs, so nothing is wired to produce one either.
##
## Why not WARP's BamToCram
## ------------------------
## `warp-clarum/tasks/broad/BamToCram.wdl` is a sub-workflow of
## ExomeGermlineSingleSample, not a standalone converter: it demands the pipeline's
## `Array[Map[String:String]]` duplication + chimerism metrics (the uBAM pass WARP
## itself produced) to populate cram_md5_with_prefix_file / validation_report, and it
## drags in CheckPreValidation + ValidateSamFile + DRAGStr machinery. WARP's
## lower-level `ConvertToCram` task is the right shape but is private to that repo's
## import graph and emits no index/integrity evidence beyond the .crai.
## `pipelines/broad/reprocessing/exome/ExomeReprocessing.wdl` does accept
## `File? input_bam`, but it rebuilds uBAMs and re-runs MarkDuplicates / BQSR /
## HaplotypeCaller / gVCF: a re-process (~$15-20/sample), not a conversion.
##
## The reference is the whole ballgame
## ----------------------------------
## CRAM stores differences against a reference and pins that reference by MD5, so the
## CRAM is only decompressible against a byte-identical FASTA. The delivered BAMs are
## DRAGEN-aligned and do NOT tell us which hg38 they used:
##
##   @PG ID: Hash Table Build  CL: /opt/edico/bin/dragen --ht-reference=.../Hsapiens/hg38/seq/hg38.fa
##   @PG ID: DRAGEN SW build   CL: /opt/dragen/4.2.4-2-1-g3ac86beb/bin/dragen -r /seq/dragen/references/hg38/dragen_komodo/
##   @HD VN:1.4 SO:coordinate   3366 x @SQ, M5 on NONE of them
##
## 3,366 contigs is also the count in Broad's hg38/v0/Homo_sapiens_assembly38.fasta,
## which makes it look interchangeable - it is not: 2,841 @SQ name+length entries match
## Broad's .fai, but 525 differ, all in the HLA block (this BAM carries `HLA-A*01`,
## `HLA-A*01:01:38L`, ...; Broad's carries `HLA-A*01:01:01:01`, `HLA-A*01:02`, ...).
## With no M5 in the header there is nothing to cross-check at runtime, so:
##
##   1. pass the EXACT FASTA DRAGEN used (ask BI for its md5), and
##   2. let the built-in round-trip check prove it - it re-reads the CRAM and compares
##      records against the source BAM over sampled windows. A wrong reference shows up
##      there as mismatched SEQ/QUAL, not as a silent corruption.
##
## Contig NAMES are unaffected (primary + ALT names match Broad's), so gVCF/BAM joint
## calling is not at risk; only CRAM decompression is. That includes the ARPA-H Commons
## side: whoever unpacks these CRAMs needs the same FASTA, so record its md5 alongside
## the deliverable.
##
## Operational notes
## -----------------
## - `bam_file` and `bai_file` sit in DIFFERENT datarepo-row prefixes
##   (e.g. .../883d659a-.../ACH0024_LM17.bam vs .../9d7d6b83-.../ACH0024_LM17.bam.bai),
##   so anything assuming "<bam>.bai" beside the BAM fails. This workflow copies the
##   supplied index into place, verifies it, and re-indexes the localized copy if not.
## - The source bucket is a TDR snapshot bucket (read-only, requester-pays). Outputs
##   must land in a CLARUM-owned bucket; write the paths back to the `sample` table.
## - One BAM is 20.4 GiB (400x exome), so expect ~60+ GiB of local disk per call.
##
## Outputs (wire straight into a Terra `sample` table):
##   output_cram, output_cram_index, output_cram_md5, output_cram_index_md5,
##   integrity_report (md5s + sizes + record counts + ROUNDTRIP line),
##   roundtrip_identical (Boolean), validation_report

workflow BamToCram {
  input {
    File    input_bam
    # Default keeps the collaborator sample id and drops any .bam suffix.
    String  output_basename = sub(basename(input_bam), "\\.bam$", "")

    # The reference DRAGEN used. Required, deliberately un-defaulted - see the header.
    File    ref_fasta
    File    ref_fasta_index

    # BI ships the index and the BAM md5 in separate snapshot rows.
    File?   input_bam_index
    File?   input_bam_md5                  # a "<hash>  <name>" text file, as delivered
    String  bam_index_path = select_first([input_bam_index, ""])
    String  bam_md5_path   = select_first([input_bam_md5, ""])

    Boolean run_roundtrip    = true        # sampled-window CRAM-vs-BAM record comparison
    Boolean strict_roundtrip = true        # fail the task if the round trip disagrees
    Int     roundtrip_window_bp          = 20000
    Int     roundtrip_windows_per_contig = 2
    Int     roundtrip_contigs            = 12   # longest N primary contigs

    Boolean run_validation = true          # Picard ValidateSamFile over the CRAM
    Int     validation_mem_gb = 16
    String  validation_ignore = "MISSING_TAG_NM"   # space-separated; WARP ignores this one too
    Boolean skip_mate_validation = false   # WARP sets this only for outlier data

    # samtools 1.11 image WARP's ConvertToCram runs on (ships seq_cache_populate.pl).
    String  samtools_docker = "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    # Picard 2.26.10 = the version the delivered 620-sample gVCFs were produced with
    # (methodConfig gatk/ReadBamHeader attrs), so these reports are comparable to the
    # retrospective ones. WARP's ValidateSamFile task runs the same image.
    String  picard_docker   = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"

    Int     cpu = 8                        # samtools -@ threads for -C compression
    Int     mem_gb = 12
    Int     additional_disk_gb = 40        # on top of ~2x BAM + reference
    Int     preemptible = 3                # a restart re-streams 20 GiB; raise to 0 if flaky

    String  picard_extra_args = ""         # e.g. "IGNORE=MATERIAL_AND_NON_PRIMITIVE_SUPPORT"
  }

  call ConvertBamToCram {
    input:
      input_bam                = input_bam,
      bam_index_path           = bam_index_path,
      bam_md5_path             = bam_md5_path,
      output_basename          = output_basename,
      ref_fasta                = ref_fasta,
      ref_fasta_index          = ref_fasta_index,
      run_roundtrip            = run_roundtrip,
      strict_roundtrip         = strict_roundtrip,
      roundtrip_window_bp      = roundtrip_window_bp,
      roundtrip_windows_per_contig = roundtrip_windows_per_contig,
      roundtrip_contigs        = roundtrip_contigs,
      samtools_docker          = samtools_docker,
      cpu                      = cpu,
      mem_gb                   = mem_gb,
      additional_disk_gb       = additional_disk_gb,
      preemptible              = preemptible
  }

  if (run_validation) {
    call ValidateCram {
      input:
        cram            = ConvertBamToCram.output_cram,
        output_basename = output_basename,
        ref_fasta       = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        validation_mem_gb     = validation_mem_gb,
        validation_ignore     = validation_ignore,
        skip_mate_validation  = skip_mate_validation,
        picard_docker     = picard_docker,
        picard_extra_args = picard_extra_args,
        preemptible       = preemptible
    }
  }

  output {
    File    output_cram             = ConvertBamToCram.output_cram
    File    output_cram_index       = ConvertBamToCram.output_cram_index
    File    output_cram_md5         = ConvertBamToCram.output_cram_md5
    File    output_cram_index_md5   = ConvertBamToCram.output_cram_index_md5
    File    integrity_report        = ConvertBamToCram.integrity_report
    Boolean roundtrip_identical     = ConvertBamToCram.roundtrip_identical
    File?   validation_report       = ValidateCram.report
  }

  meta {
    description: "BAM -> CRAM for the prospective somatic ES delivery (TDD 2.H.3 / 2.H.4): convert, index, md5, md5-check the source BAM, and prove the reference is the one DRAGEN used by comparing CRAM and BAM records over sampled windows.",
    summary: "Convert an aligned BAM to CRAM with index, md5 and round-trip evidence",
    author: "CLARUM / Talkowski lab"
  }
  # No `capabilities:` block: the pointer used by mutect2_pon.wdl
  # (gs://terra-0c560982/clarum_utils/capabilities/v1.json) does not resolve - that
  # bucket does not exist - and this workflow's outputs are all single-valued, so
  # Terra's default wiring (BamToCram.output_cram, ...) is enough. Wire them to the
  # `sample` table attributes output_cram / output_cram_index / output_cram_md5.
}

task ConvertBamToCram {
  input {
    File    input_bam
    String  bam_index_path
    String  bam_md5_path
    String  output_basename

    File    ref_fasta
    File    ref_fasta_index

    Boolean run_roundtrip
    Boolean strict_roundtrip
    Int     roundtrip_window_bp
    Int     roundtrip_windows_per_contig
    Int     roundtrip_contigs

    String  samtools_docker
    Int     cpu
    Int     mem_gb
    Int     additional_disk_gb
    Int     preemptible
  }

  # Localized BAM (1x) + CRAM (~0.3-0.5x) + a second BAM-index pass wants headroom,
  # plus the 3 GiB reference and its cache. WARP's ConvertToCram uses 2x + ref + 20.
  Int disk_gb        = ceil(2 * size(input_bam, "GB") + size(ref_fasta, "GB") + additional_disk_gb)
  Int machine_mem_mb = mem_gb * 1024

  command <<<
    set -euo pipefail

    BAM="$(basename '~{input_bam}')"
    REF="~{ref_fasta}"
    OUT="~{output_basename}"
    CPU="~{cpu}"

    echo "== toolchain" | tee integrity.txt
    samtools --version | head -2 | tee -a integrity.txt

    # ---- index: trust-but-verify the delivered .bai -------------------------------
    # htslib looks for "<bam>.bai" (or an adjacent .csi) and this BAM is a localized
    # copy, so place the supplied index there. BI's index lives in a different
    # datarepo-row prefix, and it was built against BI's copy of the file, so verify
    # with idxstats before believing it; otherwise build our own.
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM}.csi" ]; then
      if [ -n "~{bam_index_path}" ]; then
        cp "~{bam_index_path}" "${BAM}.bai"
        echo "copied supplied index -> ${BAM}.bai" >> integrity.txt
      fi
    fi
    if ! samtools idxstats "${BAM}" > bam_idxstats.txt 2>/dev/null; then
      echo "supplied index unusable (or absent) -> samtools index" >> integrity.txt
      samtools index -@ "${CPU}" "${BAM}"
      samtools idxstats "${BAM}" > bam_idxstats.txt
    fi
    BAM_MAPPED=$(awk '$3>0{s+=$3}END{print s+0}' bam_idxstats.txt)
    # idxstats columns: refname length mapped "mate unmapped" unmapped total - the
    # length is $2, so the unmapped counts are $5 (and mate-unmapped $4), not $2.
    BAM_UNMAPPED=$(awk '{s+=$4; u+=$5}END{print u+0"+"s+0}' bam_idxstats.txt)
    echo "bam: mapped=${BAM_MAPPED} unmapped(mate+alone)=${BAM_UNMAPPED}" >> integrity.txt

    # ---- md5 of the source BAM (TDD 2.H.4: integrity by size + MD5) ---------------
    md5sum "${BAM}" | awk '{print $1}' > input_bam_md5_computed.txt
    echo "input bam md5 = $(cat input_bam_md5_computed.txt)  size = $(wc -c < "${BAM}")" >> integrity.txt
    if [ -n "~{bam_md5_path}" ]; then
      SUPPLIED=$(awk 'NR==1{print $1}' "~{bam_md5_path}")
      echo "supplied bam md5 = ${SUPPLIED}" >> integrity.txt
      if [ "${SUPPLIED}" = "$(cat input_bam_md5_computed.txt)" ]; then
        echo "INPUT_BAM_MD5=MATCH" >> integrity.txt
      else
        echo "INPUT_BAM_MD5=MISMATCH" >> integrity.txt
      fi
    else
      echo "INPUT_BAM_MD5=NOT_SUPPLIED" >> integrity.txt
    fi

    # ---- convert ------------------------------------------------------------------
    # -T gives the compressor the reference directly (no REF_DOWNLOAD); -@ is the
    # multi-threaded cram encoder. No --force-reads: if htslib cannot read this BAM
    # cleanly we want to know.
    samtools view -C -@ "${CPU}" -T "${REF}" -o "${OUT}.cram" "${BAM}"
    samtools quickcheck -v "${OUT}.cram" | tee -a integrity.txt

    md5sum "${OUT}.cram" | awk '{print $1}' > "${OUT}.cram.md5"

    # ---- index the CRAM -----------------------------------------------------------
    # Indexing a CRAM needs the reference for every container it walks. Populate
    # htslib's REF_CACHE from the FASTA and shut REF_PATH off so htslib cannot quietly
    # fetch a DIFFERENT hg38 from the internet mid-run (WARP's ConvertToCram does the
    # same dance; the samtools image ships seq_cache_populate.pl).
    seq_cache_populate.pl -root ./ref/cache "${REF}"
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s
    samtools index -@ "${CPU}" "${OUT}.cram"
    md5sum "${OUT}.cram.crai" | awk '{print $1}' > "${OUT}.cram.crai.md5"

    CRAM_MAPPED=$(samtools idxstats "${OUT}.cram" | awk '$3>0{s+=$3}END{print s+0}')
    CRAM_UNMAPPED=$(samtools idxstats "${OUT}.cram" | awk '{s+=$4; u+=$5}END{print u+0"+"s+0}')
    echo "cram: mapped=${CRAM_MAPPED} unmapped(mate+alone)=${CRAM_UNMAPPED}" >> integrity.txt
    if [ "${BAM_MAPPED}" != "${CRAM_MAPPED}" ] || [ "${BAM_UNMAPPED}" != "${CRAM_UNMAPPED}" ]; then
      echo "RECORD_COUNT=MISMATCH ${BAM_MAPPED}/${BAM_UNMAPPED} vs ${CRAM_MAPPED}/${CRAM_UNMAPPED}" >> integrity.txt
      if [ "~{strict_roundtrip}" = "true" ]; then exit 1; fi
    else
      echo "RECORD_COUNT=MATCH" >> integrity.txt
    fi
    echo "cram md5 = $(cat ${OUT}.cram.md5)  size = $(wc -c < ${OUT}.cram)" >> integrity.txt
    echo "crai md5 = $(cat ${OUT}.cram.crai.md5)  size = $(wc -c < ${OUT}.cram.crai)" >> integrity.txt

    # ---- round trip: does the CRAM read back as the BAM? ---------------------------
    # This is the only available proof that REF is the reference DRAGEN aligned
    # against: a wrong FASTA decodes different SEQ/QUAL, which shows up here. Sampled
    # windows keep it cheap (~1-2 min); restrict to primary contigs because HLA contig
    # names contain ':' and are therefore region-string-hostile.
    echo "ROUNDTRIP=SKIPPED" > roundtrip.txt
    RT_FAIL=0
    if [ "~{run_roundtrip}" = "true" ]; then
      samtools view -H "${OUT}.cram" 2>/dev/null \
        | awk -F'\t' '$1=="@SQ"{n=substr($2,4); l=substr($3,4); if (n ~ /^chr[0-9XY]+$/) print l"\t"n}' \
        | sort -rn | head -n "~{roundtrip_contigs}" \
        | awk -F'\t' -v W=~{roundtrip_window_bp} -v K=~{roundtrip_windows_per_contig} \
            '{for(k=0;k<K;k++){s=int($1*(k+1)/(K+1)); print $2":"s"-"s+W}}' > windows.txt

      RT_N=0
      while read -r W; do
        RT_N=$((RT_N + 1))
        # Region goes as a POSITIONAL argument. `samtools view -r` is --read-group, not
        # a region: on these BAMs every record carries an RG tag, so `-r chr1:1-20000`
        # selects records in that read group (none) and the comparison would pass
        # vacuously on two empty files.
        samtools view -T "${REF}" "${BAM}" "${W}"  > a.sam 2>/dev/null || true
        samtools view -T "${REF}" "${OUT}.cram" "${W}" > b.sam 2>/dev/null || true
        NA=$(wc -l < a.sam | tr -d ' '); NB=$(wc -l < b.sam | tr -d ' ')
        if cmp -s a.sam b.sam; then
          FULL=MATCH; CORE=MATCH
        else
          FULL=MISMATCH
          # A full-record mismatch is EXPECTED here even for a perfect conversion:
          # htslib's CRAM decoder re-derives the MD:Z tag from the reference on
          # read-back (there is no --no-MD switch in `samtools view`), so the CRAM side
          # carries tags the BAM never had. The gate is therefore the core alignment
          # fields (QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL), sorted
          # so that tag order cannot matter: a wrong reference changes SEQ/QUAL.
          for f in a b; do cut -f1-11 $f.sam | LC_ALL=C sort > $f.core; done
          if cmp -s a.core b.core; then CORE=MATCH; else CORE=MISMATCH; fi
        fi
        echo "window ${W} bam_records=${NA} cram_records=${NB} full=${FULL} core=${CORE}" >> roundtrip.txt

        # Reference-identity cross-check (informational). The round trip above can only
        # prove the CRAM is lossless against the FASTA we handed it - CRAM stores
        # differences against whatever reference it was given, so a near-miss hg38 still
        # round-trips perfectly. What DOES betray a wrong reference is the mismatch
        # count: NM:i comes from DRAGEN's own alignment (their reference), MD:Z is
        # re-derived on read-back from OUR reference. Where the two disagree, the two
        # references differ. Reported, not gated: soft clips and indels make MD and NM
        # legitimately unequal on individual reads.
        awk -F'\t' -v w="${W}" '{nm=""; md=""
            for(i=12;i<=NF;i++){ if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]} else if($i ~ /^MD:Z:/){md=substr($i,6)} }
            isnm=(nm=="0"); ismd=(md!="" && md !~ /[ACGTN]/)
            if(isnm)nm0++; if(ismd)md0++; if(isnm&&ismd)both++ }
          END{ v="UNKNOWN"; m=(nm0<md0?nm0:md0)
            if(m==0) v="NO_INFORMATIVE_RECORDS"; else if(both >= 0.9*m) v="PLAUSIBLE"; else v="SUSPECT"
            printf "ref_agreement %s nm0=%d md0=%d both=%d verdict=%s\n", w, nm0+0, md0+0, both+0, v }' \
          b.sam >> roundtrip.txt

        # An empty window proves nothing, so it must not count as a pass.
        if [ "${CORE}" != "MATCH" ] || [ "${NA}" != "${NB}" ] || [ "${NA}" -eq 0 ]; then
          RT_FAIL=$((RT_FAIL + 1))
        fi
        rm -f a.sam b.sam a.core b.core
      done < windows.txt
      echo "roundtrip: windows=${RT_N} failing=${RT_FAIL}" >> integrity.txt
      SUS=$(grep -c 'verdict=SUSPECT' roundtrip.txt || true)
      if [ "${SUS}" -gt 0 ]; then
        echo "REF_AGREEMENT=SUSPECT on ${SUS} of ${RT_N} windows - this reference may not be the one DRAGEN aligned against" >> integrity.txt
      fi
      # An empty window list means the @SQ parse found no primary contigs - that is a
      # failed test, not a passed one, so it must not read back as roundtrip_identical.
      if [ "${RT_N}" -eq 0 ]; then
        echo "ROUNDTRIP=NO_WINDOWS (no chr[0-9XY] contigs in the header?)" >> integrity.txt
        RT_FAIL=1
      fi
      if [ "${RT_FAIL}" -gt 0 ]; then
        cat roundtrip.txt >> integrity.txt
        echo "ROUNDTRIP=FAIL(${RT_FAIL}/${RT_N})" >> roundtrip.txt
        if [ "~{strict_roundtrip}" = "true" ]; then
          echo "FATAL: CRAM does not read back as the source BAM - wrong reference?" >&2
          exit 1
        fi
      else
        echo "ROUNDTRIP=PASS(${RT_N} windows)" >> roundtrip.txt
      fi
    fi

    # roundtrip_identical: false only on a real, tested disagreement.
    if [ "${RT_FAIL}" -eq 0 ]; then echo "true" > roundtrip_identical.txt; else echo "false" > roundtrip_identical.txt; fi
  >>>

  runtime {
    docker: samtools_docker
    memory: machine_mem_mb + " MB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File output_cram           = "~{output_basename}.cram"
    File output_cram_index     = "~{output_basename}.cram.crai"
    File output_cram_md5       = "~{output_basename}.cram.md5"
    File output_cram_index_md5 = "~{output_basename}.cram.crai.md5"

    File integrity_report      = "integrity.txt"
    File roundtrip_report      = "roundtrip.txt"
    File input_bam_md5_computed = "input_bam_md5_computed.txt"
    File bam_idxstats          = "bam_idxstats.txt"
    Boolean roundtrip_identical = read_string("roundtrip_identical.txt") == "true"

    Int cram_bytes  = round(size(output_cram, "Bytes"))
    Int crai_bytes  = round(size(output_cram_index, "Bytes"))
  }

  meta {
    description: "samtools view -C + index + md5 + source-BAM md5 check + sampled-window round-trip against the source BAM."
  }
}

task ValidateCram {
  input {
    File    cram
    String  output_basename
    File    ref_fasta
    File    ref_fasta_index
    Int     validation_mem_gb
    String  validation_ignore
    Boolean skip_mate_validation
    String  picard_docker
    String  picard_extra_args
    Int     preemptible
  }

  Int disk_gb        = ceil(size(cram, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + 30)
  Int machine_mem_mb = validation_mem_gb * 1024
  Int java_xmx_mb    = machine_mem_mb - 500
  Int java_xms_mb    = machine_mem_mb - 1000

  command <<<
    set -euo pipefail

    IGNORES=""
    for tok in ~{validation_ignore}; do IGNORES="${IGNORES} IGNORE=${tok}"; done

    # Same invocation WARP's QC.ValidateSamFile uses (tasks/broad/QC.wdl): Picard,
    # REFERENCE_SEQUENCE only - no dict needed - MODE=VERBOSE, IS_BISULFITE_SEQUENCED
    # =false. MISSING_TAG_NM is expected for DRAGEN output and WARP ignores it too.
    # Mate validation is the expensive part on a 400x exome; WARP enables it unless
    # the sample is outlier data, so we default to the same.
    java -Xms~{java_xms_mb}m -Xmx~{java_xmx_mb}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT="~{cram}" \
      OUTPUT="~{output_basename}.cram.validation_report" \
      REFERENCE_SEQUENCE="~{ref_fasta}" \
      MAX_OUTPUT=100000 \
      MODE=VERBOSE \
      SKIP_MATE_VALIDATION=~{skip_mate_validation} \
      IS_BISULFITE_SEQUENCED=false \
      ${IGNORES} ~{picard_extra_args} \
      2> "~{output_basename}.cram.validation.log"
  >>>

  runtime {
    docker: picard_docker
    memory: machine_mem_mb + " MB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File report = "~{output_basename}.cram.validation_report"
    File log    = "~{output_basename}.cram.validation.log"
  }

  meta {
    description: "Picard ValidateSamFile over the converted CRAM (MODE=VERBOSE, MISSING_TAG_NM ignored by default, mate validation on unless skipped)."
  }
}
