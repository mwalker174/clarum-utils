## MOC12 2.J.6: build the deliverable package for one certified annotated callset, in the cloud (v0.1)
##
## WHY THIS EXISTS: the WGS 620 annotated VCF is 120.6 GB and does not fit on the workstation, so its
## package cannot be built where the WES and mosaic packages were. This runs the SAME script
## (src/qc/moc12_variant_annotation/package_annotated_vcf.py) that built those two, so all three
## deliverables come out of one code path (owner, 2026-10-01: "for reproducibility and to ensure
## consistency").
##
## WHAT IT DOES, all inside the script:
##   - drops any --drop INFO field still present (none, for WGS: g41_vep_consequence never got there);
##   - corrects the stale header descriptions with `bcftools reheader` -- with nothing to drop, the
##     certified file is NOT re-encoded: every compressed block after the header is copied unchanged;
##   - refuses to finish unless the records are byte-identical to the certified source (nothing
##     dropped) or identical apart from the dropped fields, the record and sample counts equal the
##     certifying run's, no header line carries a workstation path, and every INFO / CSQ / FORMAT field
##     has a description in COLUMNS.txt;
##   - writes README.txt, COLUMNS.txt, acceptance.json (sha256 + md5 of the VCF and index).
##
## TRAPS RESPECTED (AGENTS.md / clarum-utils test/replay_cromwell_layout.sh):
##   - every input the command reads is copied/linked into cwd first (Cromwell localizes each File to its
##     own mirror directory, so pkg.py would not sit beside the script otherwise);
##   - no `2>/dev/null`, no `cmd | grep -q` under pipefail; tool rc captured; stderr left intact;
##   - every output is File?, so a late failure keeps what was produced;
##   - size() units are Cromwell's ("GB").
version 1.0

workflow moc12PackageAnnotatedVcf {
  input {
    File vcf
    File source_acceptance
    File package_script
    File pkg_module
    String callset
    Array[String] drop_info = ["g41_vep_consequence"]
    String code_rev
    # the same attributes as vcf / source_acceptance, bound as Strings so the gs:// URIs reach the
    # acceptance record -- a File input only exists in the task as a localized path
    String source_vcf_uri
    String source_acceptance_uri
    String docker = "python:3.12-slim"
    Int cpu = 8
    Int memory_gb = 16
  }

  call packageVcf {
    input:
      vcf = vcf,
      source_acceptance = source_acceptance,
      package_script = package_script,
      pkg_module = pkg_module,
      callset = callset,
      drop_info = drop_info,
      code_rev = code_rev,
      source_vcf_uri = source_vcf_uri,
      source_acceptance_uri = source_acceptance_uri,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb
  }

  output {
    File? package_vcf = packageVcf.package_vcf
    File? package_vcf_tbi = packageVcf.package_vcf_tbi
    File? acceptance_json = packageVcf.acceptance_json
    File? columns_txt = packageVcf.columns_txt
    File? readme_txt = packageVcf.readme_txt
  }

  meta {
    description: "MOC12 2.J.6 deliverable package for one certified annotated VCF (same script as the local packages)"
  }
}

task packageVcf {
  input {
    File vcf
    File source_acceptance
    File package_script
    File pkg_module
    String callset
    Array[String] drop_info
    String code_rev
    String source_vcf_uri
    String source_acceptance_uri
    String docker
    Int cpu
    Int memory_gb
  }

  # the certified input + one output copy + index and headroom
  Int disk_gb = ceil(size(vcf, "GB") * 2.2 + 20)
  String out_vcf = "clarum_~{callset}_annotated.vcf.gz"

  command <<<
    set -euo pipefail
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends bcftools tabix ca-certificates
    bcftools --version | sed -n '1,2p'

    cp '~{package_script}' ./package_annotated_vcf.py
    cp '~{pkg_module}' ./pkg.py
    cp '~{source_acceptance}' ./source_acceptance.json

    drop_args=""
    for d in ~{sep=" " drop_info}; do drop_args="$drop_args --drop $d"; done

    rc=0
    python3 ./package_annotated_vcf.py --callset '~{callset}' --vcf '~{vcf}' \
      --source-acceptance source_acceptance.json $drop_args \
      --out-dir pkg --code-rev '~{code_rev}' --threads ~{cpu} \
      --record-source-vcf '~{source_vcf_uri}' --record-source-acceptance '~{source_acceptance_uri}' || rc=$?
    echo "package_rc=$rc"
    if [ "$rc" -ne 0 ]; then exit "$rc"; fi
    ls -l pkg
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    # SSD: the job is two full sequential reads and one write of a 120 GB file; the WGS gather on HDD
    # spent ~8.5 h on exactly this kind of I/O (docs/progress/097 §4.11)
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 0
    maxRetries: 1
  }

  output {
    File? package_vcf = "pkg/~{out_vcf}"
    File? package_vcf_tbi = "pkg/~{out_vcf}.tbi"
    File? acceptance_json = "pkg/acceptance.json"
    File? columns_txt = "pkg/COLUMNS.txt"
    File? readme_txt = "pkg/README.txt"
  }
}
