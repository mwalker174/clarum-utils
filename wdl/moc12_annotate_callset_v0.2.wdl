## MOC12 2.J.6: gnomAD v4.1 + ClinVar onto a VEP + clinical-layer annotated callset (v0.2)
##
## INPUT: the parts `vepAnnotateHailExtra` wrote (WES: ~99 merged chunks; WGS: 25 JG shards; mosaic: 1).
## Each part runs, in one task:
##   1. drop the clinical layer's own ClinVar INFO tags (CLNSIG/CLNREVSTAT/CLNSIGCONF/GENEINFO, a
##      March 2025 release with no accession) -- owner ruled 2026-09-25: one ClinVar source, ours;
##   2. annotate_gnomad_v41.py  -> g41_* (remote indexed tabix on the public gnomAD bucket);
##   3. annotate_clinvar.py     -> clinvar_* (a localized release, with VCV accession), and VEP's
##      CSQ CLIN_SIG emptied in place (position kept), so exactly one ClinVar source ships;
## then one gather concatenates the parts with bcftools and writes the acceptance JSON.
##
## WHAT CHANGED FROM v0.1, AND WHY -- all six v0.1 submissions on Terra failed (2026-09-22/23):
##   - v0.1 passed shard URIs as String and streamed them from gs:// inside the task; one run died
##     `OSError: truncated file` mid-read. v0.2 takes each part as a localized File: there are
##     at most ~100 parts now, not 1,769, so the reason for streaming is gone.
##   - v0.1 died opening `gnomad.exomes.v4.1.sites.chrM.vcf.bgz`, which gnomAD does not publish.
##     annotate_gnomad_v41.py now proves absence with a 404 from the bucket's metadata API and passes
##     those records through as NOT MEASURED (a transient failure still raises).
##   - one task per part, not per contig group: a part may straddle contigs (docs/progress/096 §17.4),
##     and the script already reads any contig it meets.
##   - gather is `bcftools concat --threads`, not a pysam rewrite: the WGS callset is 47.7M records.
##
## TRAPS RESPECTED (AGENTS.md / external/clarum-utils/test/replay_cromwell_layout.sh):
##   - meta{} entries newline-separated, no trailing comma.
##   - no path built from basename() of a localized input; files the command reads are copied into cwd
##     first, because Cromwell localizes each File to its own mirror directory (so a .tbi is NOT next
##     to its .vcf.gz, and a script's sibling module is NOT next to the script).
##   - no `2>/dev/null`, no `head` in a pipeline; tool rc captured, stderr left intact.
##   - every output a gate could skip is File?, so a late failure keeps what was produced.
##   - size() units inside Cromwell's vocabulary ("GB").
version 1.0

workflow moc12AnnotateCallset {
  input {
    Array[File] part_vcfs
    File annotate_gnomad_script
    File annotate_clinvar_script
    File provenance_script
    File gather_script
    File clinvar_vcf
    File clinvar_vcf_tbi
    Array[String] drop_info = ["CLNSIG", "CLNREVSTAT", "CLNSIGCONF", "GENEINFO"]
    # VEP's own CSQ CLIN_SIG is a third ClinVar source (its cache's release); values emptied,
    # position kept, so one ClinVar source ships (owner, 2026-09-25)
    Array[String] blank_csq = ["CLIN_SIG"]
    String deliverable_name
    String gnomad_dataset = "joint"
    # Owner ruled 2026-09-21: plain AF, max across gnomAD ancestry groups. Passed explicitly so the
    # deliverable records the choice instead of inheriting a script default.
    String spine = "af_grpmax"
    Float maf = 0.01
    Int window = 100000
    String docker = "python:3.12-slim"
    Int part_cpu = 2
    Int part_memory_gb = 8
    Int part_preemptible_tries = 2
    Int gather_cpu = 8
    Boolean gather = true
  }

  scatter (i in range(length(part_vcfs))) {
    call annotatePart {
      input:
        vcf = part_vcfs[i],
        part_name = "~{deliverable_name}.part~{i}",
        annotate_gnomad_script = annotate_gnomad_script,
        annotate_clinvar_script = annotate_clinvar_script,
        provenance_script = provenance_script,
        clinvar_vcf = clinvar_vcf,
        clinvar_vcf_tbi = clinvar_vcf_tbi,
        drop_info = drop_info,
        blank_csq = blank_csq,
        gnomad_dataset = gnomad_dataset,
        spine = spine,
        maf = maf,
        window = window,
        docker = docker,
        cpu = part_cpu,
        memory_gb = part_memory_gb,
        preemptible_tries = part_preemptible_tries
    }
  }

  if (gather) {
    call gatherCallset {
      input:
        part_vcfs = select_all(annotatePart.annotated_vcf),
        part_tbis = select_all(annotatePart.annotated_vcf_tbi),
        gnomad_summaries = select_all(annotatePart.gnomad_summary),
        clinvar_summaries = select_all(annotatePart.clinvar_summary),
        g41_manifest = select_first(annotatePart.g41_manifest),
        clinvar_manifest = select_first(annotatePart.clinvar_manifest),
        gather_script = gather_script,
        drop_info = drop_info,
        blank_csq = blank_csq,
        deliverable_name = deliverable_name,
        docker = docker,
        cpu = gather_cpu
    }
  }

  output {
    Array[File?] part_annotated_vcfs = annotatePart.annotated_vcf
    Array[File?] part_gnomad_summaries = annotatePart.gnomad_summary
    Array[File?] part_clinvar_summaries = annotatePart.clinvar_summary
    File? annotated_vcf = gatherCallset.annotated_vcf
    File? annotated_vcf_tbi = gatherCallset.annotated_vcf_tbi
    File? acceptance_json = gatherCallset.acceptance_json
  }

  meta {
    description: "MOC12 2.J.6 full-callset gnomAD v4.1 + ClinVar join onto VEP-annotated parts"
    annotations: ["gnomAD-v4.1", "ClinVar-VCV", "MOC12", "three-state-verdict"]
  }
}

## ---------------------------------------------------------------------------

task annotatePart {
  input {
    File vcf
    String part_name
    File annotate_gnomad_script
    File annotate_clinvar_script
    File provenance_script
    File clinvar_vcf
    File clinvar_vcf_tbi
    Array[String] drop_info
    Array[String] blank_csq
    String gnomad_dataset
    String spine
    Float maf
    Int window
    String docker
    Int cpu
    Int memory_gb
    Int preemptible_tries
  }

  # input + three rewrites of it (dropped, +g41, +clinvar) must fit at once
  Int disk_gb = ceil(size(vcf, "GB") * 5 + size(clinvar_vcf, "GB") + 20)
  String out_vcf = "~{part_name}.annotated.vcf.gz"
  String out_gnomad = "~{part_name}.gnomad_join_summary.json"
  String out_clinvar = "~{part_name}.clinvar_join_summary.json"

  command <<<
    set -euo pipefail
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends libcurl4 ca-certificates bcftools tabix
    pip install --no-cache-dir -q pysam==0.24.0

    if [ -z "${GCS_OAUTH_TOKEN:-}" ]; then
      # python3, not curl: the python:3.12-slim image ships no curl binary (v0.1 canary, rc 127)
      GCS_OAUTH_TOKEN="$(python3 -c "
import json, urllib.request
r = urllib.request.Request(
    'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token',
    headers={'Metadata-Flavor': 'Google'})
print(json.load(urllib.request.urlopen(r, timeout=30))['access_token'])
")"
      export GCS_OAUTH_TOKEN
    fi

    # Scripts and the ClinVar pair into cwd: each File lands in its own mirror dir, so neither
    # `import provenance` nor the .tbi would be found beside its partner otherwise.
    cp '~{annotate_gnomad_script}' ./annotate_gnomad_v41.py
    cp '~{annotate_clinvar_script}' ./annotate_clinvar.py
    cp '~{provenance_script}' ./provenance.py
    cp '~{clinvar_vcf}' ./clinvar.vcf.gz
    cp '~{clinvar_vcf_tbi}' ./clinvar.vcf.gz.tbi

    n_in=$(bcftools view -H --threads ~{cpu} '~{vcf}' | wc -l | tr -d ' ')
    echo "part=~{part_name} input_records=$n_in"

    # 1. drop the clinical layer's ClinVar tags; only those present are named, so a part without
    #    them (e.g. a layer that skipped ClinVar) is not an error, and what was removed is logged
    # Header to a FILE first. `bcftools view -h | grep -q` under pipefail is false even on a match:
    # grep -q exits at the first hit, bcftools takes SIGPIPE, the pipeline returns non-zero -- the
    # first cloud canary (59a5b298) logged "dropping: <none present>" on parts that carried all four.
    bcftools view -h '~{vcf}' > in.header
    present=""
    for t in ~{sep=" " drop_info}; do
      if grep -q "^##INFO=<ID=$t," in.header; then present="$present,INFO/$t"; fi
    done
    present="${present#,}"
    echo "dropping: ${present:-<none present>}"
    if [ -n "$present" ]; then
      bcftools annotate --threads ~{cpu} -x "$present" -Oz -o part.vcf.gz '~{vcf}'
    else
      cp '~{vcf}' part.vcf.gz
    fi

    # 2. gnomAD v4.1
    rc=0
    python3 ./annotate_gnomad_v41.py --vcf part.vcf.gz \
      --dataset '~{gnomad_dataset}' --spine '~{spine}' --maf '~{maf}' --window '~{window}' \
      --out-dir g41 || rc=$?
    echo "gnomad_rc=$rc"
    if [ "$rc" -ne 0 ]; then exit "$rc"; fi
    rm -f part.vcf.gz

    # 3. ClinVar (ours)
    blank_args=""
    for f in ~{sep=" " blank_csq}; do blank_args="$blank_args --blank-csq $f"; done
    python3 ./annotate_clinvar.py --vcf g41/part.g41.vcf.gz --clinvar clinvar.vcf.gz \
      $blank_args --out-dir cv || rc=$?
    echo "clinvar_rc=$rc"
    if [ "$rc" -ne 0 ]; then exit "$rc"; fi

    mv cv/part.g41.clinvar.vcf.gz '~{out_vcf}'
    tabix -f -p vcf '~{out_vcf}'
    cp g41/gnomad_join_summary.json '~{out_gnomad}'
    cp cv/clinvar_join_summary.json '~{out_clinvar}'

    n_out=$(bcftools index -n '~{out_vcf}')
    echo "output_records=$n_out"
    if [ "$n_out" != "$n_in" ]; then
      echo "RECORD COUNT MISMATCH: in=$n_in out=$n_out"
      exit 7
    fi
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    maxRetries: 1
  }

  output {
    File? annotated_vcf = out_vcf
    File? annotated_vcf_tbi = out_vcf + ".tbi"
    File? gnomad_summary = out_gnomad
    File? clinvar_summary = out_clinvar
    File? g41_manifest = "g41/metrics_manifest.tsv"
    File? clinvar_manifest = "cv/clinvar_metrics_manifest.tsv"
  }
}

task gatherCallset {
  input {
    Array[File] part_vcfs
    Array[File] part_tbis
    Array[File] gnomad_summaries
    Array[File] clinvar_summaries
    File g41_manifest
    File clinvar_manifest
    File gather_script
    Array[String] drop_info
    Array[String] blank_csq
    String deliverable_name
    String docker
    Int cpu
  }

  Int disk_gb = ceil(size(part_vcfs, "GB") * 2.5 + 20)
  String out_vcf = "~{deliverable_name}.annotated.vcf.gz"
  String out_acceptance = "~{deliverable_name}.acceptance.json"

  command <<<
    set -euo pipefail
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends libcurl4 ca-certificates bcftools tabix
    pip install --no-cache-dir -q pysam==0.24.0

    # parts and their indexes into one dir, so each .tbi sits beside its VCF
    mkdir -p parts
    while read -r f; do ln -s "$f" "parts/$(basename "$f")"; done < '~{write_lines(part_vcfs)}'
    while read -r f; do ln -s "$f" "parts/$(basename "$f")"; done < '~{write_lines(part_tbis)}'
    ls -1 parts/*.vcf.gz > parts.list
    echo "parts=$(wc -l < parts.list | tr -d ' ')"

    python3 '~{gather_script}' \
      --parts parts.list \
      --gnomad-summaries '~{write_lines(gnomad_summaries)}' \
      --clinvar-summaries '~{write_lines(clinvar_summaries)}' \
      --g41-manifest '~{g41_manifest}' \
      --clinvar-manifest '~{clinvar_manifest}' \
      --dropped '~{sep="," drop_info}' \
      --blanked-csq '~{sep="," blank_csq}' \
      --out-vcf '~{out_vcf}' \
      --acceptance '~{out_acceptance}' \
      --threads ~{cpu}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "16 GB"
    disks: "local-disk ~{disk_gb} HDD"
    # a gather is not worth a preemption coin flip
    preemptible: 0
    maxRetries: 1
  }

  output {
    File? annotated_vcf = out_vcf
    File? annotated_vcf_tbi = out_vcf + ".tbi"
    File? acceptance_json = out_acceptance
  }
}
