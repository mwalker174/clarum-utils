## MOC12 clinical annotation: gnomAD v4.1 join, grouped by contig  (v0.1)
##
## THE CLAUSE
##   "Generate a cloud-based annotation VCF with >40 different metrics from the VCF and external
##   databases (e.g. gnomAD allele frequencies)" -- cloud-based is in the wording, and the metric count
##   is the stated success criterion (docs/progress/094 §7).
##
## TWO MEASURED FACTS THAT DICTATED THIS SHAPE (2026-09-18, gsutil stat / pysam)
##   1. gnomAD v4.1 per-chromosome sites VCFs are huge -- **joint chr1 = 72.06 GB compressed**, joint
##      chr7 = 47.5 GB, genomes chr1 = 44.1 GB, exomes chr1 = 18.8 GB, because they carry age/GQ/DP/AB
##      histograms and per-ancestry-per-sex blocks. So gnomAD must NOT be a localized File input: a
##      contig scatter would try to pull hundreds of GB onto task disks. It is read by **indexed remote
##      tabix**, with a bearer token minted inside the task from the VM identity -- the same mechanism
##      the join uses locally, where 4,409 records took 90 s and pulled a trivial number of bytes.
##   2. Per-shard scatter = ~1,769 VMs for the delivered WES callset. The join is single-threaded and
##      I/O-bound, so work is grouped **by contig** (~24 tasks): groupByContig -> annotateContig ->
##      gatherDeliverable.
##
## WHY SHARD URIS TRAVEL AS STRINGS
##   Passing them as `String` rather than `File` is what stops Cromwell localizing 1,769 shards into one
##   task. Stated cost of that choice rather than hidden: a mistyped URI fails at run time inside the
##   task instead of at localization time, and call caching keys on the string list.
##
## WHAT THE TASK RUNS, AND THE ONE SCRIPT CHANGE v0.1 NEEDS
##   annotate_gnomad_v41.py (src/qc/moc12_variant_annotation/), staged as a File, called through its
##   proven interface only: --vcf (repeatable) --dataset --maf --window --max-records --out-dir.
##   It emits our own `g41_*` columns, so a VEP cache bump (which renames gnomADe_oth -> gnomADe_remaining)
##   cannot silently move a clinical filter, and a three-state PASS / FAIL / NOT MEASURED verdict, so a
##   variant gnomAD does not know is never counted as rare.
##   The script currently imports `paths` from src/qc/lib for its default out-dir; in a task only the
##   script file is present, so it must treat that import as optional when --out-dir is supplied.
##
## v0.1 SCOPE IS DELIBERATELY THE JOIN ONLY
##   ClinVar accession, the 86-gene panel filter and the deliverable README statistics are v0.2. Keeping
##   the filter a separate stage is the economics: Eren's full-callset annotation run measured **$127.90**
##   (docs/progress/094 §8), so re-thresholding the panel or the MAF boundary must not require re-running
##   an annotation. Annotate once, filter many, for cents.
##
## TRAPS RESPECTED HERE (see external/clarum-utils/test/replay_cromwell_layout.sh)
##   - meta{} entries newline-separated, no trailing comma (Cromwell rejects one; miniwdl accepts it).
##   - No basename() on a localized path: output names are written once as task-level Strings.
##   - No `2>/dev/null`, no `head` in a pipeline under `set -euo pipefail`; tool exit codes are captured
##     as data, and stderr survives so a failure is diagnosable from the platform alone.
##   - Anything a gate could skip is an optional output (File?), so a late failure cannot delete what the
##     call already produced (canary #2 lost a 21.9 GB CRAM to a required-output hash cut).
##   - size()/disk units are inside Cromwell's vocabulary; nothing says "Bytes".
##   - Three-state results are a String verdict, never a nullable Boolean (miniwdl rejects None / read_string?).

version 1.1

## ---------------------------------------------------------------------------
## Shared shape: every task needs (a) pysam and (b) a GCS bearer token for gs:// range reads.
## ---------------------------------------------------------------------------

task groupByContig {
  input {
    Array[String] shard_uris
    # A group is sized by wall-clock, not memory: 4,409 records measured 90 s locally, so ~40 shards of
    # the delivered WES callset per task stays comfortably inside a preemptible window.
    Int max_shards_per_group = 40
    File group_script
    String annotate_docker = "python:3.12-slim"
  }

  File shards_in = write_lines(shard_uris)

  command <<<
    set -euo pipefail
    # pysam's manylinux wheel bundles htslib, whose gs:// support links libcurl; the slim image does not
    # necessarily ship libcurl4. Both installs are unpinned-silent on purpose: a missing dependency here
    # later masquerades as a network failure, so let it fail loudly at install time.
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends libcurl4 ca-certificates
    pip install --no-cache-dir pysam==0.24.0

    if [ -z "${GCS_OAUTH_TOKEN:-}" ]; then
      minted="$(curl -fsS -H 'Metadata-Flavor: Google' \
        'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token')"
      GCS_OAUTH_TOKEN="$(printf '%s' "$minted" | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')"
      export GCS_OAUTH_TOKEN
    fi

    python3 '~{group_script}' '~{shards_in}' '~{max_shards_per_group}' > groups.tsv
    wc -l < groups.tsv | tr -d ' ' | sed 's/^/groups=/'

    cat > /dev/null <<'UNUSED'
    UNUSED
  >>>

  runtime {
    docker: annotate_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }

  output {
    # Fixed three columns -- contig, group index, then every shard URI of the group tab-joined. Ragged
    # TSV would make read_tsv's shape depend on shard counts, which is a bad thing for a config to
    # depend on.
    Array[Array[String]] groups = read_tsv("groups.tsv")
  }

  meta {
    calls_caching: true
    purpose: "decide which shards share one task, by contig, without localizing any of them"
  }
  parameter_meta {
    shard_uris: "gs:// URIs of the callset's per-interval shards"
    max_shards_per_group: "group size cap; drives task wall-clock"
  }
}

## ---------------------------------------------------------------------------

task annotateContig {
  input {
    String contig
    # Tab-joined shard URIs for this contig, in coordinate order, as emitted by groupByContig.
    String shard_blob
    File annotate_script
    File gather_contig_script
    String gnomad_dataset = "joint"
    Float maf = 0.01
    Int window = 100000
    Int max_records = 0
    String cohort_prefix = "clarum"
    String annotate_docker = "python:3.12-slim"
    Int cpu_cores = 2
    Int memory_gb = 8
    Int disk_gb = 200
    Int preemptible_tries = 2
  }

  File shards_file = write_lines([shard_blob])
  String out_vcf = "~{cohort_prefix}.~{contig}.g41.vcf.gz"
  String out_summary = "~{cohort_prefix}.~{contig}.gnomad_join_summary.json"
  String out_manifest = "~{cohort_prefix}.~{contig}.metrics_manifest.tsv"

  command <<<
    set -euo pipefail
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends libcurl4 ca-certificates
    pip install --no-cache-dir pysam==0.24.0

    if [ -z "${GCS_OAUTH_TOKEN:-}" ]; then
      minted="$(curl -fsS -H 'Metadata-Flavor: Google' \
        'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token')"
      GCS_OAUTH_TOKEN="$(printf '%s' "$minted" | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')"
      export GCS_OAUTH_TOKEN
    fi

    # One --vcf per shard. The URIs come from a file, not from an interpolated command line, so a
    # large group cannot turn into an E2BIG discovered only at scale.
    tr '\t' '\n' < '~{shards_file}' | sed '/^$/d' > shards.txt
    vcf_args=""
    while read -r uri; do
      vcf_args="$vcf_args --vcf $uri"
    done < shards.txt
    echo "shards_in_group=$(wc -l < shards.txt | tr -d ' ')"

    mkdir -p out
    # `|| rc=$?` keeps a nonzero tool exit as data instead of an unannotated abort, and leaves stderr
    # intact for the platform log.
    rc=0
    python3 '~{annotate_script}' $vcf_args \
      --dataset '~{gnomad_dataset}' \
      --maf '~{maf}' \
      --window '~{window}' \
      --max-records '~{max_records}' \
      --out-dir out || rc=$?
    echo "annotate_rc=$rc"
    if [ "$rc" -ne 0 ]; then
      echo "ANNOTATION FAILED rc=$rc for contig '~{contig}'"
      exit "$rc"
    fi

    python3 '~{gather_contig_script}' \
      --out-dir out \
      --out-vcf '~{out_vcf}' \
      --out-summary '~{out_summary}' \
      --contig '~{contig}' \
      --dataset '~{gnomad_dataset}' \
      --maf '~{maf}'

    ls -l '~{out_vcf}' '~{out_summary}' '~{out_manifest}'
  >>>

  runtime {
    docker: annotate_docker
    cpu: cpu_cores
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    maxRetries: 1
  }

  output {
    File annotated_vcf = out_vcf
    File summary_json = out_summary
    File metrics_manifest = out_manifest
  }

  meta {
    calls_caching: true
    cost_note: "gnomAD is streamed by tabix, never localized; joint chr1 alone is 72.06 GB compressed"
    cache_rationale: "inputs are the callset shards + gnomAD release, neither of which moves when the gene panel or MAF boundary changes"
  }
  parameter_meta {
    contig: "the single contig this group covers, e.g. chr1"
    shard_blob: "tab-joined gs:// shard URIs for that contig, coordinate-ordered"
    annotate_script: "annotate_gnomad_v41.py staged to a readable gs:// prefix"
    gnomad_dataset: "rarity spine: joint | exomes | genomes"
    maf: "rarity boundary on popmax FAF95 (falling back to AF); 0.01 per the 2026-09-17 PM ruling"
  }
}

## ---------------------------------------------------------------------------

task gatherDeliverable {
  input {
    Array[File] contig_vcfs
    Array[File] contig_summaries
    File metrics_manifest
    File gather_deliverable_script
    String deliverable_name = "clarum_annotated"
    String annotate_docker = "python:3.12-slim"
    Int disk_gb = 300
    Int memory_gb = 8
  }

  String out_vcf = "~{deliverable_name}.annotated.vcf.gz"
  String out_acceptance = "~{deliverable_name}.acceptance.json"

  command <<<
    set -euo pipefail
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends libcurl4 ca-certificates
    pip install --no-cache-dir pysam==0.24.0

    python3 '~{gather_deliverable_script}' \
      --parts '~{sep="\t" contig_vcfs}' \
      --out-vcf '~{out_vcf}' \
      --manifest '~{metrics_manifest}' \
      --acceptance '~{out_acceptance}' \
      --summaries '~{sep="\t" contig_summaries}'
  >>>

  runtime {
    docker: annotate_docker
    cpu: 2
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    File annotated_vcf = out_vcf
    File acceptance_json = out_acceptance
  }

  meta {
    calls_caching: false
    cache_rationale: "this is the artifact that states whether the '>40 metrics' clause is met; recompute it whenever the definition of done changes -- it is minutes of work"
  }
}

## ---------------------------------------------------------------------------

workflow moc12AnnotateMetrics {
  input {
    Array[String] callset_shard_uris
    File annotate_script
    File group_script
    File gather_contig_script
    File gather_deliverable_script
    String cohort_prefix = "clarum"
    String gnomad_dataset = "joint"
    Float maf = 0.01
    Int window = 100000
    Int max_shards_per_group = 40
    String deliverable_name = "clarum_annotated"
    String annotate_docker = "python:3.12-slim"
    Int memory_gb = 8
    Int disk_gb = 200
    Boolean gather = true
  }

  call groupByContig {
    input:
      shard_uris             = callset_shard_uris,
      max_shards_per_group   = max_shards_per_group,
      group_script           = group_script,
      annotate_docker        = annotate_docker,
  }

  scatter (grp in groupByContig.groups) {
    call annotateContig {
      input:
        contig          = grp[0],
        shard_blob      = grp[2],
        annotate_script = annotate_script,
        gather_contig_script = gather_contig_script,
        gnomad_dataset  = gnomad_dataset,
        maf             = maf,
        window          = window,
        cohort_prefix   = cohort_prefix,
        annotate_docker = annotate_docker,
        memory_gb       = memory_gb,
        disk_gb         = disk_gb,
    }
  }

  Array[File] per_contig_vcfs = annotateContig.annotated_vcf
  Array[File] per_contig_summaries = annotateContig.summary_json

  if (gather) {
    call gatherDeliverable {
      input:
        contig_vcfs      = per_contig_vcfs,
        contig_summaries = per_contig_summaries,
        metrics_manifest = annotateContig.metrics_manifest[0],
        gather_deliverable_script = gather_deliverable_script,
        deliverable_name = deliverable_name,
        annotate_docker  = annotate_docker,
    }
  }

  output {
    Array[File] contig_annotated_vcfs = per_contig_vcfs
    Array[File] contig_summaries = per_contig_summaries
    Array[File] metrics_manifests = annotateContig.metrics_manifest
    File? annotated_vcf = gatherDeliverable.annotated_vcf
    File? acceptance_json = gatherDeliverable.acceptance_json
  }

  meta {
    annotations: ["gnomAD-v4.1", "MOC12", "route-A", "three-state-verdict", "annotate-once-filter-many"]
  }
}
