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
## The reference is the whole ballgame - and it is settled
## -------------------------------------------------------
## CRAM stores differences against a reference and pins that reference by MD5, so the
## CRAM is only decompressible against a byte-identical FASTA. The delivered ES BAMs are
## DRAGEN-aligned and do NOT tell us which hg38 they used:
##
##   @PG ID: Hash Table Build  CL: /opt/edico/bin/dragen --ht-reference=.../Hsapiens/hg38/seq/hg38.fa
##   @PG ID: DRAGEN SW build   CL: /opt/dragen/4.2.4-2-1-g3ac86beb/bin/dragen -r /seq/dragen/references/hg38/dragen_komodo/
##   @HD VN:1.4 SO:coordinate   3366 x @SQ, M5 on NONE of them
##
## An earlier draft of this header claimed the delivered contig dictionary disagreed
## with Broad's on 525 HLA names. That was wrong (a diff against a missing file).
## Measured against `Homo_sapiens_assembly38.fasta`: ALL 3,366 @SQ name+length entries
## match, in the same order, and every one of those 3,366 contigs carries reads - so the
## full ALT/HLA/decoy build is required, not a no-alt analysis set.
##
## The reference is `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta`
## (file md5 7ff134953dcca8c8997453bbb80b6b5e, per-sequence md5 22f6161a7683386226afe1831a4fafc6).
## Three independent lines of evidence, all in `docs/progress/070`:
##
##   1. The BI-delivered GS CRAMs pin per-contig M5 for 3,366/3,366 contigs, and those
##      values equal this FASTA's 3,366/3,366 - across three cohorts (ACH, BCH, CSG),
##      which also agree with each other. Table: `data/qc-inputs/bi_gs_cram_sq_m5.tsv`,
##      rebuilt with `src/scripts/reference_m5_table.py`.
##   2. Genome-wide NM audit of this very BAM vs this FASTA: 55,531 pure-M reads over
##      3,783 sampled windows on 1,203 contigs (primary, HLA, ALT, decoy, chrUn) -
##      99.888% reproduce DRAGEN's own NM exactly. Artifacts under
##      `artifacts/qc/moc11_prospective_es/cram/`, tool `src/scripts/verify_reference_against_bam.py`.
##   3. A local BAM->CRAM->SAM round trip of 616,215 real records from this BAM decoded
##      byte-for-byte as the source (core fields + full tag set).
##
## What the 0.112% NM disagreement is
## ----------------------------------
## Three windows out of 3,783 disagree systematically, and each disagreement sits on an
## exact dbSNP ALT allele (chr22:20337902 rs28411685 C>A, chr15:20398748 rs1846757 A>G,
## chr16:18076851 rs199761340 C>T) carried by reads with NM=0. The same site in unrelated
## BCH samples charges NM=2..4 for the same base. That is the pop-alt/graph path in the
## HT-build @PG (`--ht-pop-alt-contigs ... 32global.af.05... --ht-pop-alt-liftover ...`):
## a read aligned on a pop-alt contig is re-emitted at primary coordinates with NM scored
## on that contig. It is DRAGEN's NM bookkeeping, not a different FASTA - so **do not gate
## conversion on NM or MD**. An earlier version of this task inferred reference identity
## from NM-vs-MD agreement; that heuristic fires on exactly these legitimate sites, so it
## is replaced by the deterministic `reference_m5_table` comparison below.
##
## The near-miss hazard is real and it is partial, not total
## ---------------------------------------------------------
## Decoding a CRAM built on this FASTA against a copy with ONE chrM base changed gives
## `MD5 checksum reference mismatch at chrM:827-1002` and exit 1 - but only after
## emitting ~20k reads from earlier slices, of which some came back with altered SEQ.
## So: never pipe a decode into `head`, never `|| true` over its exit status, and never
## trust output from a run that did not exit 0. `REF_PATH=:` is exported so htslib cannot
## fetch a *different* hg38 from EBI mid-run instead of failing.
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
    # md5 of ref_fasta itself (broad hg38/v0 FASTA = 7ff134953dcca8c8997453bbb80b6b5e).
    # Set it: it turns "we pointed at a hg38" into "we pointed at THE hg38 we audited".
    String  ref_fasta_md5 = ""
    # Per-contig M5 the delivery expects (contig/length/m5 TSV, see reference_m5_table.py).
    # Our CRAM's own @SQ M5 must equal it, else the CRAM is pinned to a reference the
    # Commons cannot reproduce. Optional only because the file has to be staged in a
    # CLARUM bucket first.
    File?   reference_m5_table

    # BI ships the index and the BAM md5 in separate snapshot rows.
    #
    # Keep these File? ALL THE WAY INTO THE TASK. They used to be aliased
    #   String bam_index_path = select_first([input_bam_index, ""])
    # which is a File->String cast evaluated in the WORKFLOW, where no localization
    # exists yet: it yields the original gs:// URI, and because nothing File-typed
    # reaches the task, Cromwell does not localize those two objects at all. Both
    # halves were proved on submission 58698ceb-782a-46a2-9384-89f063edb192 - the
    # rendered task script contained `cp "gs://datarepo-e6f844bb-bucket/…bai"` (no
    # gsutil in the samtools image) and its gcs_localization.sh listed only the BAM
    # and the reference. Passing the File? through unchanged is what makes Cromwell
    # localize it and interpolate an in-container path.
    File?   input_bam_index
    File?   input_bam_md5                  # a "<hash>  <name>" text file, as delivered

    Boolean run_roundtrip    = true        # sampled-window CRAM-vs-BAM record comparison
    Boolean strict_roundtrip = true        # fail the task if the round trip disagrees
    Int     roundtrip_window_bp          = 20000
    Int     roundtrip_windows_per_contig = 2
    Int     roundtrip_contigs            = 12   # longest N primary contigs
    # How many sampled windows must actually contain reads for the comparison to mean
    # anything. Off-target windows on capture data are skipped, not failed (see the
    # round-trip block), so this is what fails a run that never got to compare anything.
    Int     min_roundtrip_windows        = 8

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
      input_bam_index          = input_bam_index,
      input_bam_md5            = input_bam_md5,
      output_basename          = output_basename,
      ref_fasta                = ref_fasta,
      ref_fasta_index          = ref_fasta_index,
      ref_fasta_md5            = ref_fasta_md5,
      reference_m5_table       = reference_m5_table,
      run_roundtrip            = run_roundtrip,
      strict_roundtrip         = strict_roundtrip,
      roundtrip_window_bp      = roundtrip_window_bp,
      roundtrip_windows_per_contig = roundtrip_windows_per_contig,
      roundtrip_contigs        = roundtrip_contigs,
      min_roundtrip_windows    = min_roundtrip_windows,
      samtools_docker          = samtools_docker,
      cpu                      = cpu,
      mem_gb                   = mem_gb,
      additional_disk_gb       = additional_disk_gb,
      preemptible              = preemptible
  }

  # `defined(...)` because ConvertBamToCram's outputs are optional (see its completion
  # manifest): a call that died before converting must not launch Picard on a null CRAM.
  if (run_validation && defined(ConvertBamToCram.output_cram)) {
    call ValidateCram {
      input:
        cram            = select_first([ConvertBamToCram.output_cram]),
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

  # Everything is optional at this boundary on purpose. Cromwell aborts a call's
  # delocalization list at the first missing REQUIRED file and skips everything behind it
  # in an order the WDL does not control; submission 0c6ccba5-1045-4479-be9e-0b75548801f7
  # therefore lost a CRAM that had already been converted and written, because one report
  # the strict gate had not reached sat earlier in the list. Completeness is enforced
  # inside the task (completion manifest) instead, where a missing artifact is a loud
  # failure rather than a silent stranding. Nulls here mean "that evidence was not
  # produced", which a caller must still check.
  output {
    File?    output_cram             = ConvertBamToCram.output_cram
    File?    output_cram_index       = ConvertBamToCram.output_cram_index
    File?    output_cram_md5         = ConvertBamToCram.output_cram_md5
    File?    output_cram_index_md5   = ConvertBamToCram.output_cram_index_md5
    File?    integrity_report        = ConvertBamToCram.integrity_report
    File?    roundtrip_report        = ConvertBamToCram.roundtrip_report
    File?    reference_m5_check      = ConvertBamToCram.reference_m5_check
    File?    reference_sq_table      = ConvertBamToCram.reference_sq_table
    String   roundtrip_verdict       = ConvertBamToCram.roundtrip_verdict
    Boolean  roundtrip_tested        = ConvertBamToCram.roundtrip_tested
    Boolean  roundtrip_identical     = ConvertBamToCram.roundtrip_identical
    File?    validation_report       = ValidateCram.report
    File?    validation_log          = ValidateCram.log
    Boolean? validation_passed       = ValidateCram.passed
  }

  ## NOTE on the `meta` block: Terra's Cromwell parser rejects commas between entries
  ## ("Expected rbrace, got ','" when the method is registered) - `meta` is not a WDL
  ## expression map. Keep one entry per line, no trailing commas; miniwdl accepts that too.
  meta {
    description: "BAM -> CRAM for the prospective somatic ES delivery (TDD 2.H.3 / 2.H.4): convert, index, md5, md5-check the source BAM, and prove the reference is the one DRAGEN used by comparing CRAM and BAM records over sampled windows."
    summary: "Convert an aligned BAM to CRAM with index, md5 and round-trip evidence"
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
    File?   input_bam_index
    File?   input_bam_md5
    String  output_basename

    File    ref_fasta
    File    ref_fasta_index
    String  ref_fasta_md5
    File?   reference_m5_table

    Boolean run_roundtrip
    Boolean strict_roundtrip
    Int     roundtrip_window_bp
    Int     roundtrip_windows_per_contig
    Int     roundtrip_contigs
    Int     min_roundtrip_windows

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

    # Cromwell/GCP materializes File inputs at <call_root>/<bucket>/<object prefix>/<name>
    # - NOT in the task's working directory - so `basename "~{input_bam}"` is a dangling
    # RELATIVE path on the cloud even though the identical line works locally, where the
    # runner drops the file into cwd under its own name. That is what killed submission
    # 58698ceb-782a-46a2-9384-89f063edb192 (rc 1, 0-byte stderr, stdout stopping after the
    # samtools version): the first command touching "${BAM}" died under `set -euo pipefail`
    # and its 2>/dev/null ate the message. Keep the localized path and link it into cwd,
    # which is exactly what the reference block below already does.
    BAMSRC='~{input_bam}'
    BAM="$(basename "$BAMSRC")"
    if [ ! -s "$BAMSRC" ]; then
      echo "FATAL: input BAM is not localized at $BAMSRC" >&2
      exit 1
    fi
    # Hardlink first (same filesystem as the reference block above, so a 22 GB BAM costs
    # no second copy), symlink if that is refused, copy as a last resort. The existence
    # guard keeps this correct for runners that already put the file in cwd under its own
    # name (miniwdl does) - `ln -f` there would replace a real file with a link to itself.
    if [ ! -e "./$BAM" ]; then
      ln "$BAMSRC" "./$BAM" 2>/dev/null \
        || ln -s "$BAMSRC" "./$BAM" 2>/dev/null \
        || cp "$BAMSRC" "./$BAM"
    fi
    OUT="~{output_basename}"
    CPU="~{cpu}"

    # Optional inputs: interpolate the File? itself (never a pre-cast String) so Cromwell
    # localizes it and hands back the in-container path. `default=""` is the Cromwell
    # idiom for a File? in a command; miniwdl parses it with a deprecation notice only.
    IDX='~{default="" input_bam_index}'
    MD5IN='~{default="" input_bam_md5}'
    # Loud, not silent, if that ever degrades to a cloud URI again - the exact failure
    # this task already paid 47 minutes of coldline localization for.
    for probe in "$IDX" "$MD5IN"; do
      case "$probe" in
        gs://*|https://*)
          echo "FATAL: optional input arrived as an unlocalized cloud URI: $probe" >&2
          exit 1 ;;
      esac
    done

    echo "== toolchain" | tee integrity.txt
    # `| head -2` is a pipefail landmine (head closes the pipe -> SIGPIPE 141 on samtools);
    # sed reads to EOF, so the pipeline cannot fail on a version banner.
    samtools --version | sed -n '1,2p' | tee -a integrity.txt

    # ---- reference: prove which hg38 this is, before spending an hour on it --------
    # htslib only checks reference md5 when it decodes a slice, so a wrong FASTA can
    # burn the whole conversion before it says anything. Check the file md5 up front.
    #
    # Also re-link it under its canonical name and hand samtools THAT path: htslib
    # writes @SQ UR from the -T argument, so this is what the delivered CRAM will say
    # about where its reference lives. The default would be Cromwell's per-run
    # localization directory, which is meaningless after the job exits (the BI-delivered
    # CRAMs carry M5 and no UR at all; identity there is md5-only, as it is here).
    REFSRC="~{ref_fasta}"
    REFBASE="$(basename "$REFSRC")"
    # Existence-guarded for the same reason as the BAM block: a runner that already put the
    # file in cwd under its own name (miniwdl) makes a bare `ln` fail and the `cp` fallback
    # then refuses to copy a file onto itself.
    if [ ! -e "./$REFBASE" ]; then
      ln "$REFSRC" "./$REFBASE" 2>/dev/null || cp "$REFSRC" "./$REFBASE"
    fi
    if [ ! -e "./${REFBASE}.fai" ]; then
      ln "~{ref_fasta_index}" "./${REFBASE}.fai" 2>/dev/null || cp "~{ref_fasta_index}" "./${REFBASE}.fai"
    fi
    REF="./$REFBASE"
    REF_MD5=$(md5sum "$REF" | awk '{print $1}')
    echo "reference = ${REFBASE} md5=${REF_MD5} size=$(wc -c < "$REF")" >> integrity.txt
    if [ -n "~{ref_fasta_md5}" ] && [ "${REF_MD5}" != "~{ref_fasta_md5}" ]; then
      echo "FATAL: reference md5 ${REF_MD5} != expected ~{ref_fasta_md5}" >&2
      exit 1
    fi
    if [ ! -s "${REF}.fai" ]; then
      echo "FATAL: reference index ${REFBASE}.fai is empty" >&2
      exit 1
    fi
    N_FA=$(cut -f1 "${REF}.fai" | sort -u | wc -l | tr -d ' ')
    # No 2>/dev/null here: a header read that fails must be visible, and under `set -e`
    # a swallowed stderr is how a whole task dies with nothing to read afterwards.
    N_SQ=$(samtools view -H "${BAM}" | awk -F'\t' '$1=="@SQ"{n++}END{print n+0}')
    echo "contigs: fasta=${N_FA} bam_@SQ=${N_SQ}" >> integrity.txt
    if [ "${N_SQ}" -eq 0 ]; then
      echo "FATAL: read 0 @SQ lines from ${BAM} (empty/truncated localization?)" >&2
      exit 1
    fi
    if [ "${N_FA}" != "${N_SQ}" ]; then
      echo "WARNING: reference contig count (${N_FA}) != BAM @SQ count (${N_SQ}) - CRAM encode will fail or skip contigs" >> integrity.txt
    fi

    # ---- index: trust-but-verify the delivered .bai -------------------------------
    # htslib looks for "<bam>.bai" (or an adjacent .csi) and this BAM is a localized
    # copy, so place the supplied index there. BI's index lives in a different
    # datarepo-row prefix, and it was built against BI's copy of the file, so verify
    # with idxstats before believing it; otherwise build our own.
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM}.csi" ]; then
      if [ -n "$IDX" ]; then
        cp "$IDX" "${BAM}.bai"
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
    if [ -n "$MD5IN" ]; then
      SUPPLIED=$(awk 'NR==1{print $1}' "$MD5IN")
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
    # cleanly we want to know. version=3.0 explicitly: it is what the delivered GS CRAMs
    # are and what every consumer here can read; 4.0 would be a unilateral change.
    samtools view -C -@ "${CPU}" -O cram,version=3.0 -T "${REF}" -o "${OUT}.cram" "${BAM}"
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

    # ---- reference identity: what this CRAM now pins ------------------------------
    # The @SQ M5 written at encode time IS the identity claim a consumer will check.
    # If a table of expected M5 is supplied, compare against it and make disagreement
    # fatal: a CRAM pinned to a reference the ARPA-H Commons cannot reproduce is not a
    # deliverable, however faithfully it round-trips locally. Built from the BI-delivered
    # GS CRAM headers by src/scripts/reference_m5_table.py.
    samtools view -H "${OUT}.cram" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $1=="@SQ"{sn="";ln="";m5="";
          for(i=2;i<=NF;i++){t=substr($i,1,2); v=substr($i,4);
            if(t=="SN")sn=v; else if(t=="LN")ln=v; else if(t=="M5")m5=tolower(v)}
          if(sn!="")print sn,ln,m5}' > cram_sq.tsv
    N_SQ_M5=$(awk 'NF>=3 && $3!=""' cram_sq.tsv | wc -l | tr -d ' ')
    echo "cram @SQ=$(wc -l < cram_sq.tsv | tr -d ' ') with M5=${N_SQ_M5}" >> integrity.txt
    # Array name `want`, NOT `exp`: `exp` is awk's builtin exponential function, and using
    # it as an array name is invalid - the first time this block was ever executed (local
    # replay, docs/progress/072) BSD awk died with "awk: illegal statement at source line 2".
    # It had never run on Terra because the canary shipped without reference_m5_table.
    M5_TABLE="~{default="" reference_m5_table}"
    if [ -n "${M5_TABLE}" ]; then
      awk -F'\t' -v t="${M5_TABLE}" '
        NR==FNR { if ($0 ~ /^#/) next; want[$1]=$3; len[$1]=$2; n++; next }
        { e = want[$1]
          if (e == "")      { absent++; if (absent<4) print "  not in table: " $1 > "/dev/stderr" }
          else if ($3 != e) { bad++;  printf "  %s expected=%s got=%s\n", $1, e, $3 > "/dev/stderr" }
          else              { ok++ }
          if (len[$1] != "" && len[$1] != $2) lbad++ }
        END { printf "REFERENCE_M5=%s ok=%d mismatched=%d not_in_table=%d length_mismatch=%d (table has %d contigs)\n",
                (bad+lbad>0 ? "MISMATCH" : "MATCH"), ok+0, bad+0, absent+0, lbad+0, n }' \
        "${M5_TABLE}" cram_sq.tsv > m5_check.txt 2>> m5_check.txt
      cat m5_check.txt >> integrity.txt
      if grep -q 'REFERENCE_M5=MISMATCH' m5_check.txt; then
        echo "FATAL: CRAM @SQ M5 disagrees with the expected reference table" >&2
        cat m5_check.txt >&2
        exit 1
      fi
    else
      echo "REFERENCE_M5=NOT_CHECKED (no reference_m5_table supplied)" >> integrity.txt
    fi

    # ---- round trip: does the CRAM read back as the BAM? ---------------------------
    # Sampled windows keep it cheap (~1-2 min). Restricted to primary contigs because
    # HLA contig NAMES contain ':' and are therefore hostile in a region string; the
    # HLA/ALT/decoy block is covered offline by src/scripts/verify_reference_against_bam.py.
    RT_N=0
    RT_TESTED=0
    RT_EMPTY=0
    RT_FAIL=0
    if [ "~{run_roundtrip}" != "true" ]; then
      echo "ROUNDTRIP=SKIPPED (run_roundtrip=false)" > roundtrip.txt
      # Three states, never two: "not tested" must not arrive as true or false. Written as
      # text and decoded in the output block as TWO booleans (roundtrip_tested /
      # roundtrip_identical), because a portable WDL 1.0 null literal - `None`, or the
      # `read_string?` optional-function form - is precisely what Cromwell and miniwdl
      # parse differently, and a deliverable flag must not depend on that.
      echo "not_tested" > roundtrip_identical.txt
    else
      : > roundtrip.err
      samtools view -H "${OUT}.cram" 2>/dev/null \
        | awk -F'\t' '$1=="@SQ"{n=substr($2,4); l=substr($3,4); if (n ~ /^chr[0-9XY]+$/) print l"\t"n}' \
        | sort -rn | sed -n "1,~{roundtrip_contigs}p" \
        | awk -F'\t' -v W=~{roundtrip_window_bp} -v K=~{roundtrip_windows_per_contig} \
            '{for(k=0;k<K;k++){s=int($1*(k+1)/(K+1)); print $2":"s"-"s+W}}' > windows.txt

      while read -r W; do
        RT_N=$((RT_N + 1))
        # Region goes as a POSITIONAL argument. `samtools view -r` is --read-group, not
        # a region: on these BAMs every record carries an RG tag, so `-r chr1:1-20000`
        # selects records in that read group (none) and the comparison would pass
        # vacuously on two empty files.
        RCA=0; RCB=0
        # Exit status is part of the evidence: htslib reports a reference md5 mismatch
        # per SLICE, after it has already emitted the reads of earlier slices. Swallowing
        # the failure (|| true) and comparing two truncated files would read as a mere
        # record-count difference and hide the diagnosis.
        samtools view -T "${REF}" "${BAM}" "${W}" > a.sam 2>> roundtrip.err                || RCA=$?
        samtools view -T "${REF}" "${OUT}.cram" "${W}" > b.sam 2>> roundtrip.err           || RCB=$?
        NA=$(wc -l < a.sam | tr -d ' '); NB=$(wc -l < b.sam | tr -d ' ')
        DERIVED=0; SOURCE_TAGS=0
        if cmp -s a.sam b.sam; then
          FULL=MATCH; CORE=MATCH; TAGS=MATCH
        else
          FULL=MISMATCH
          # A byte-identical SAM is NOT expected even for a perfect conversion, for two
          # benign reasons measured on the delivered BAM (616,215 records, docs/progress/070):
          #   - htslib re-derives MD:Z on read-back in some slices (the source has none:
          #     DRAGEN ran with generate-md-tags=false);
          #   - CRAM stores tags in its own order.
          # So: (core) the 11 alignment fields, sorted - a wrong reference changes SEQ/QUAL;
          #      (tags) every tag VALUE in the window as a sorted field multiset, minus the
          #              keys htslib derives from the reference itself.
          # NM:i: belongs with MD:Z: here, and that is not a nicety: the delivered BAMs
          # carry neither tag (generate-md-tags=false; Picard's MISSING_TAG_NM complaint in
          # the TDD is the same absence), while CRAM decoding adds BOTH - measured locally
          # as 11 SAM fields in / 13 out, core identical. Comparing them anyway made every
          # window tags=MISMATCH, which strict_roundtrip=true then aborts the task over.
          # Dropping them costs no signal: both are computed from SEQ against the
          # reference, and a wrong reference shows up in the core comparison first.
          #
          # awk, not `cut -f12- | tr '\t' '\n'`: on a record with no optional fields at all
          # (which is what a tagless source BAM gives), `cut -f12-` prints an EMPTY LINE on
          # both BSD and GNU cut - measured locally, 6 blank lines for a 6-read window - so
          # the two multisets differ by one blank per read and the gate reports MISMATCH on
          # a conversion that is actually exact. awk over fields 12..NF has no such case.
          for f in a b; do
            cut -f1-11 $f.sam | LC_ALL=C sort > $f.core
            { cut -f1-11 $f.sam
              awk -F'\t' '{for(i=12;i<=NF;i++) if ($i !~ /^MD:Z:/ && $i !~ /^NM:i:/) print $i}' $f.sam; } \
              | LC_ALL=C sort > $f.tags
          done
          if cmp -s a.core b.core; then CORE=MATCH; else CORE=MISMATCH; fi
          if cmp -s a.tags  b.tags;  then TAGS=MATCH;  else TAGS=MISMATCH;  fi
          # Reported, never gated: the source's own tag count and what CRAM added back.
          # On BI data expect source_tags=0 and cram_derived_tags ~= 2 reads x 2 tags.
          SOURCE_TAGS=$(awk -F'\t' '{if(NF>11) n+=NF-11}END{print n+0}' a.sam)
          DERIVED=$(awk -F'\t' '{for(i=12;i<=NF;i++) if ($i ~ /^MD:Z:/ || $i ~ /^NM:i:/) c++}END{print c+0}' b.sam)
        fi
        echo "window ${W} bam_records=${NA} cram_records=${NB} bam_exit=${RCA} cram_exit=${RCB} full=${FULL} core=${CORE} tags=${TAGS} source_tags=${SOURCE_TAGS} cram_derived_tags=${DERIVED}" >> roundtrip.txt

        # An empty window is neither a pass nor a failure: it is the absence of a test.
        # On capture data a midpoint/quarter-point window lands off-target routinely -
        # canary 0c6ccba5 drew 1 empty window in 24 on ACH0024_LM17 while the other 23
        # came back core=MATCH tags=MATCH, and the strict gate aborted a 2 h conversion
        # over that empty window. Counted separately; min_roundtrip_windows below is what
        # fails a run that could not sample enough reads to compare anything.
        if [ "${NA}" -eq 0 ] && [ "${NB}" -eq 0 ] && [ "${RCA}" -eq 0 ] && [ "${RCB}" -eq 0 ]; then
          RT_EMPTY=$((RT_EMPTY + 1))
        else
          RT_TESTED=$((RT_TESTED + 1))
          if [ "${CORE}" != "MATCH" ] || [ "${TAGS}" != "MATCH" ] || [ "${NA}" != "${NB}" ] \
             || [ "${RCA}" -ne 0 ] || [ "${RCB}" -ne 0 ]; then
            RT_FAIL=$((RT_FAIL + 1))
          fi
        fi
        rm -f a.sam b.sam a.core b.core a.tags b.tags
      done < windows.txt
      echo "roundtrip: windows=${RT_N} tested=${RT_TESTED} empty=${RT_EMPTY} failing=${RT_FAIL}" >> integrity.txt
      if [ -s roundtrip.err ]; then
        echo "roundtrip stderr (first 20 lines):" >> integrity.txt
        head -20 roundtrip.err >> integrity.txt
      fi

      # One verdict, written to both reports, and worded for what was actually measured.
      # NO_WINDOWS means the @SQ parse found no primary contigs; INSUFFICIENT_COVERAGE
      # means the windows held too few reads. Both are inconclusive, therefore failures -
      # but they are reported as themselves, not as "the CRAM disagrees".
      RT_BAD=0
      if [ "${RT_N}" -eq 0 ]; then
        RT_VERDICT="ROUNDTRIP=NO_WINDOWS (no chr[0-9XY] contigs in the CRAM header?)"
        RT_BAD=1
      elif [ "${RT_TESTED}" -lt "~{min_roundtrip_windows}" ]; then
        RT_VERDICT="ROUNDTRIP=INSUFFICIENT_COVERAGE (tested=${RT_TESTED} < min_roundtrip_windows=~{min_roundtrip_windows}, sampled=${RT_N}, empty=${RT_EMPTY})"
        RT_BAD=1
      elif [ "${RT_FAIL}" -gt 0 ]; then
        RT_VERDICT="ROUNDTRIP=FAIL(${RT_FAIL}/${RT_TESTED} tested, ${RT_EMPTY} empty)"
        RT_BAD=1
      else
        RT_VERDICT="ROUNDTRIP=PASS(${RT_TESTED} tested of ${RT_N} sampled, ${RT_EMPTY} empty)"
      fi
      echo "${RT_VERDICT}" >> roundtrip.txt
      echo "${RT_VERDICT}" >> integrity.txt

      # The boolean is written BEFORE any strict exit. On submission 0c6ccba5 the
      # `exit 1` below skipped this file, Cromwell stopped delocalizing at the first
      # missing required output, and the converted CRAM plus integrity.txt - which had
      # existed since the first second of the call - were never uploaded.
      if [ "${RT_BAD}" -eq 0 ]; then echo "true" > roundtrip_identical.txt; else echo "false" > roundtrip_identical.txt; fi

      if [ "${RT_BAD}" -eq 1 ] && [ "~{strict_roundtrip}" = "true" ]; then
        # Deliberately NOT "wrong reference?": a CRAM decodes against whatever reference
        # it is handed and reproduces what it stored even if that reference is the wrong
        # build, so this gate detects a reference that changed between encode and decode,
        # not a wrong build. A wrong build is the @SQ M5 gate above (or the offline join
        # of cram_sq.tsv against docs/RESOURCES.md's M5 table) that detects that.
        echo "FATAL: CRAM does not read back as the source BAM over ${RT_TESTED} tested window(s) - see roundtrip.txt" >&2
        grep -v 'core=MATCH tags=MATCH' roundtrip.txt >&2 || true
        exit 1
      fi
    fi

    # ---- completion manifest -----------------------------------------------------------
    # Every declared output of this task is File? (see the workflow output block), which
    # removes the stranding hazard but also removes Cromwell's own guarantee that a
    # successful call produced them. This manifest is that guarantee: a call that exits 0
    # has produced every artifact the deliverable needs, and a call that is missing one
    # says which, loudly, while uploading whatever does exist.
    MANIFEST="${OUT}.cram ${OUT}.cram.crai ${OUT}.cram.md5 ${OUT}.cram.crai.md5 integrity.txt bam_idxstats.txt input_bam_md5_computed.txt cram_sq.tsv roundtrip_identical.txt"
    MISSING=""
    for f in ${MANIFEST}; do
      [ -s "${f}" ] || MISSING="${MISSING} ${f}"
    done
    if [ -n "${MISSING}" ]; then
      echo "FATAL: task finished without producing:${MISSING}" >&2
      exit 1
    fi
  >>>

  runtime {
    docker: samtools_docker
    memory: machine_mem_mb + " MB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  # All optional: see the completion manifest above and the workflow output comment.
  # `bam_idxstats` stays here as an output too - it is what a window landing off-target
  # gets checked against, so it is part of the evidence a reviewer needs.
  output {
    File? output_cram           = "~{output_basename}.cram"
    File? output_cram_index     = "~{output_basename}.cram.crai"
    File? output_cram_md5       = "~{output_basename}.cram.md5"
    File? output_cram_index_md5 = "~{output_basename}.cram.crai.md5"

    File? integrity_report      = "integrity.txt"
    File? roundtrip_report      = "roundtrip.txt"
    File? reference_m5_check    = "m5_check.txt"
    File? reference_sq_table    = "cram_sq.tsv"
    File? input_bam_md5_computed = "input_bam_md5_computed.txt"
    File? bam_idxstats          = "bam_idxstats.txt"
    # Three states, kept distinct without a null literal (see the shell comment above):
    #   roundtrip_tested=false                -> the comparison never ran
    #   roundtrip_tested=true  identical=false -> tested and disagreed
    #   roundtrip_tested=true  identical=true  -> tested and every window matched
    String  roundtrip_verdict   = read_string("roundtrip_identical.txt")
    Boolean roundtrip_tested    = roundtrip_verdict == "true" || roundtrip_verdict == "false"
    Boolean roundtrip_identical = roundtrip_verdict == "true"

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
    REPORT="~{output_basename}.cram.validation_report"
    LOG="~{output_basename}.cram.validation.log"
    # Picard exits non-zero when it counts ERRORS, which is information, not a reason to
    # throw away the conversion: capture the rc instead of letting `set -e` kill a call
    # that has already produced a valid CRAM (and, as in ConvertBamToCram's manifest
    # comment, a task that dies with a declared-but-absent output strands the other
    # outputs behind it in Cromwell's delocalization list). The rc travels with the
    # report as `passed`, and reading the report stays a human judgement.
    RC=0
    java -Xms~{java_xms_mb}m -Xmx~{java_xmx_mb}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT="~{cram}" \
      OUTPUT="${REPORT}" \
      REFERENCE_SEQUENCE="~{ref_fasta}" \
      MAX_OUTPUT=100000 \
      MODE=VERBOSE \
      SKIP_MATE_VALIDATION=~{skip_mate_validation} \
      IS_BISULFITE_SEQUENCED=false \
      ${IGNORES} ~{picard_extra_args} \
      2> "${LOG}" || RC=$?
    echo "VALIDATE_RC=${RC}" > validate_rc.txt
    echo "ValidateSamFile rc=${RC}" >> "${LOG}"
    # Both declared outputs must exist whatever Picard did: it writes no report at all
    # when it dies on its own (no sequence dictionary, unreadable reference, ...).
    if [ ! -s "${REPORT}" ]; then
      {
        echo "ValidateSamFile produced no report (rc=${RC}). Log tail:"
        tail -40 "${LOG}" 2>/dev/null || true
      } > "${REPORT}"
    fi
  >>>

  runtime {
    docker: picard_docker
    memory: machine_mem_mb + " MB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File    report  = "~{output_basename}.cram.validation_report"
    File    log     = "~{output_basename}.cram.validation.log"
    String exit_note = read_string("validate_rc.txt")
    Boolean passed   = exit_note == "VALIDATE_RC=0"
  }

  meta {
    description: "Picard ValidateSamFile over the converted CRAM (MODE=VERBOSE, MISSING_TAG_NM ignored by default, mate validation on unless skipped)."
  }
}
