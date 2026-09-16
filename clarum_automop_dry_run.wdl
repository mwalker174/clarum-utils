version 1.0

# clarum_automop_dry_run -- a self-contained, DRY-RUN-ONLY re-implementation of
# talkowski-lab/lr-annotation's wdl/tools/Automop.wdl, which is itself a wrapper
# around `fissfc --yes --verbose mop`.
#
# Why this is a wrapper of our own rather than their file:
#   * Theirs imports ../utils/Helpers.wdl (138 KB) and ../utils/Structs.wdl but uses
#     only RuntimeAttr, so the imports are dead weight for registration into AGORA
#     (a single-WDL store). Only RuntimeAttr is inlined here.
#   * Theirs takes a required `automop_docker` with no default, and the image tag is
#     not published anywhere we can read (automop's README: "reach out to Michael
#     Gatzen"). Pinning our own image removes the guess.
#   * The image MUST be Python <= 3.11: firecloud/fccore.py:137 calls
#     configparser.SafeConfigParser(), removed in Python 3.12, so a modern base image
#     gives you a fissfc that dies on import. Measured in this repo's own py3.14 venv.
#   * Theirs writes a `mop_events` row to broad-dsde-methods-automop.automop.mop_events
#     AFTER deleting, so a workspace whose pet SA cannot write that table reports the
#     workflow as Failed once the data is already gone. Removed here.
#
# WHY THIS CANNOT DELETE -- three independent reasons, do not "fix" them:
#   1. python:3.11-slim has no `gsutil`, and `fiss mop` shells out to
#      `gsutil -m rm -I` as its ONLY delete path (firecloud/fiss.py, after the
#      `if args.dry_run ... return 0` guard at line ~1500). No gsutil => no delete.
#   2. `dry_run` defaults to true and the command refuses to run at all unless
#      `confirm_delete` is also set to true, so a mis-typed JSON cannot silently
#      flip it.
#   3. `fissfc mop` has no submission-status gate and no last-copy rule. Measured on
#      this estate (docs/progress/087 §10.12): it would delete 284.72 TiB, including
#      102.54 TiB that is the only copy of its content and 0.43 TiB under a run that
#      was still Submitted at capture. That is why the audited TSV pipeline, not this
#      workflow, is what should ever remove files.
#
# Intended use: run it once per workspace with dry_run=true and read the log. The
# per-file list + "Total Size:" line is an INDEPENDENT recount of what "unreferenced
# inside submissions/" means, computed from the live data tables rather than a
# snapshot, so it cross-checks our inventory instead of confirming it.
#
# Outputs are optional (`File?`): Cromwell delocalizes outputs in hash order and stops
# at the first missing required file, which would throw away the very log you need to
# diagnose a failure (docs/progress/071 trap #8).

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow ClarumAutomopDryRun {
    input {
        String workspace_namespace
        String workspace_name

        # true = print only. Set false ONLY together with confirm_delete, and note
        # reason 1 above: this image has no gsutil, so it cannot delete regardless.
        Boolean dry_run = true
        Boolean confirm_delete = false

        # operator identity, recorded in the log header only (no BigQuery write)
        String user

        String prefix = "clarum_automop_dry_run"

        String fiss_version = "0.16.39"
        String python_image = "python:3.11-slim"

        RuntimeAttr? runtime_attr_run_mop
    }

    call RunMop {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name      = workspace_name,
            user                = user,
            dry_run             = dry_run,
            confirm_delete      = confirm_delete,
            prefix              = prefix,
            fiss_version        = fiss_version,
            docker              = python_image,
            runtime_attr_override = runtime_attr_run_mop
    }

    output {
        File?  fiss_log     = RunMop.fiss_log
        File?  total_size   = RunMop.total_size
        File?  stderr_log   = RunMop.stderr_log
        String mop_status   = RunMop.mop_status
    }
}

task RunMop {
    input {
        String  workspace_namespace
        String  workspace_name
        String  user
        Boolean dry_run
        Boolean confirm_delete
        String  prefix
        String  fiss_version
        String  docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
      set -euo pipefail

      # Gate 2: never run the delete-capable shape without a deliberate second flag.
      if [ "~{dry_run}" = "false" ] && [ "~{confirm_delete}" != "true" ]; then
          echo "REFUSING: dry_run=false requires confirm_delete=true"
          echo "  (and this image ships no gsutil, so fiss mop cannot delete anyway;"
          echo "   see the WDL header -- use the audited TSV pipeline for deletion)"
          exit 3
      fi

      DRY_FLAG="--dry-run"
      if [ "~{dry_run}" = "false" ]; then DRY_FLAG=""; fi

      echo "operator=~{user} mode=auto-mop dry_run=~{dry_run} workspace=~{workspace_namespace}/~{workspace_name}"

      # Pin the interpreter contract, then fail loudly BEFORE analysing anything if
      # fissfc is unusable (Python >= 3.12 + fiss <= 0.16.x => SafeConfigParser error).
      python3 -V
      pip install --quiet "firecloud==~{fiss_version}" google-cloud-storage
      if ! fissfc --help >/dev/null 2>&1; then
          echo "FATAL: fissfc ~{fiss_version} will not run on this interpreter"
          fissfc --help 2>&1 | sed -n '1,12p' || true
          exit 4
      fi
      echo "fissfc import OK on $(python3 -V 2>&1)"

      # `fissfc --yes --verbose mop` is the exact verb the lab's Automop.wdl runs.
      # -p takes the workspace NAMESPACE; --dry-run returns before gsutil rm
      # (firecloud/fiss.py: `if args.dry_run or (not args.yes and not _confirm_prompt(...)): return 0`).
      set +e
      fissfc --yes --verbose mop -w "~{workspace_name}" -p "~{workspace_namespace}" $DRY_FLAG \
          > "~{prefix}.fiss_mop.log" 2>&1
      RC=$?
      set -e
      echo "fissfc_mop_rc=$RC" | tee -a "~{prefix}.fiss_mop.log"

      # Parse the two lines Automop.wdl itself keys on, so the result is machine
      # readable without shipping 100 MB of stdout around.
      if grep -q "^No files to mop in" "~{prefix}.fiss_mop.log"; then
          MOP_STATUS="no_files_to_mop"
          echo "0" > "~{prefix}.total_size.txt"
      elif grep -q "^Operation completed over" "~{prefix}.fiss_mop.log"; then
          MOP_STATUS="operation_completed_DELETE_HAPPENED"
      elif grep -q "^Total Size: " "~{prefix}.fiss_mop.log"; then
          MOP_STATUS="dry_run_total_reported"
      elif [ "$RC" -ne 0 ]; then
          if grep -q "FileNotFoundError.*gsutil\|gsutil: command not found" "~{prefix}.fiss_mop.log"; then
              # Expected if anyone ever sets dry_run=false on this image: fiss's only
              # delete verb is `gsutil -m rm -I`, and this image does not ship gsutil.
              MOP_STATUS="refused_no_gsutil_delete_path"
          else
              MOP_STATUS="fissfc_failed_rc_$RC"
          fi
      else
          MOP_STATUS="unrecognised_fiss_output"
      fi
      grep "^Total Size: " "~{prefix}.fiss_mop.log" | tail -1 > "~{prefix}.total_size.txt" || true
      echo "mop_status=$MOP_STATUS dry_run=~{dry_run}" >> "~{prefix}.total_size.txt"
      echo "$MOP_STATUS" > "~{prefix}.mop_status.txt"
      echo "mop_status=$MOP_STATUS"
    >>>

    output {
        File? fiss_log    = "~{prefix}.fiss_mop.log"
        File? total_size  = "~{prefix}.total_size.txt"
        File? stderr_log  = "~{prefix}.mop_status.txt"
        String mop_status = read_string("~{prefix}.mop_status.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 0,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:       select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:    select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks:     "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker:    select_first([runtime_attr.docker, docker])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
