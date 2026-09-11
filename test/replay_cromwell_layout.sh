#!/usr/bin/env bash
# Replay ConvertBamToCram's command block against a Cromwell-shaped file layout, without
# Cromwell, Docker, or the 21.9 GB BAM.
#
# WHY THIS EXISTS. The task "worked locally" and still failed on Terra with rc=1, an empty
# stderr, and stdout stopping after the samtools version banner. Reason: Cromwell/GCP Batch
# materializes File inputs at <call_root>/<bucket>/<object prefix>/<name> while the task
# runs *in* <call_root>, so a path built with basename() is dangling there but fine in a
# local run, which drops inputs into cwd under its own name. Nothing about `miniwdl check`
# or a local happy-path run can see that (docs/progress/072, 073).
#
# The command was pared down to the bare essentials (convert, index, md5; reference md5 as
# the only gate), so this harness now asserts exactly that surface:
#   1. the slim command passes end-to-end in the cloud layout and the four deliverables
#      are internally consistent (md5s recompute, record counts match the source BAM);
#   2. the reference-md5 gate bites: a wrong expected md5 fails before any conversion;
#   3. the output contract still declares exactly the four deliverables, with units
#      Cromwell can parse (two canaries died in output evaluation on this - cf4cb95f,
#      0c6ccba5 - and miniwdl accepts both offenders).
#
# Usage:  test/replay_cromwell_layout.sh [--wdl path/to/bam_to_cram.wdl]
# Requires: samtools, python3, awk, sed. Reference-free.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WDL="$HERE/../bam_to_cram.wdl"
while [ $# -gt 0 ]; do
  case "$1" in
    --wdl) WDL="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 64 ;;
  esac
done
[ -f "$WDL" ] || { echo "no such WDL: $WDL" >&2; exit 66; }
command -v samtools >/dev/null || { echo "samtools required" >&2; exit 69; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cramreplay.XXXXXX")"
CR="$WORK/cromwell_root"
BIN="$WORK/bin"
mkdir -p "$CR/datarepo-fake/uuid1" \
         "$CR/gcp-public-data--broad-references/hg38/v0" "$BIN"
trap 'rm -rf "$WORK"' EXIT

# macOS has no md5sum (the samtools image ships one): give the replay a faithful stand-in.
if ! command -v md5sum >/dev/null; then
  printf '#!/bin/sh\nmd5 -r "$@" | awk "{print \$1\"  \" \$2}"\n' > "$BIN/md5sum"
  chmod +x "$BIN/md5sum"
fi
MD5=md5sum; command -v md5sum >/dev/null || MD5="$BIN/md5sum"

# A 5 kb single-contig reference and a coordinate-sorted BAM of 47 reads, laid out the way
# Cromwell lays files out: the BAM in its own bucket/prefix dir, the reference in its own.
REF="$CR/gcp-public-data--broad-references/hg38/v0/mini.fasta"
python3 - "$REF" <<'PY'
import random, sys
random.seed(7)
seq = "".join(random.choice("ACGT") for _ in range(5000))
with open(sys.argv[1], "w") as f:
    f.write(">chr1 mini\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
PY
samtools faidx "$REF"
REF_MD5=$($MD5 "$REF" | awk '{print $1}')

Q="$(printf 'I%.0s' $(seq 100))"
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\n'
  for p in $(seq 200 100 4800); do
    printf 'r%s\t99\tchr1\t%s\t60\t100M\t=\t%s\t300\t%s\t%s\tBC:Z:ACTGACTG-1\tBZ:Z:AACC,GGTT,ACGT\tBI:i:1\tQX:i:99\tSM:i:1\tGP:i:100\tNM:i:0\n' \
      "$p" "$p" "$((p+300))" \
      "$(python3 -c "import random;random.seed($p);print(''.join(random.choice('ACGT') for _ in range(100)))")" \
      "$Q"
  done
} > "$WORK/body.sam"
{ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\n'
  grep -v '^@' "$WORK/body.sam" | sort -k3,3 -k4,4n; } > "$WORK/a.sam"
samtools sort -o "$CR/datarepo-fake/uuid1/A.bam" "$WORK/a.sam"

python3 "$HERE/render_command.py" "$WDL" "$WORK/task.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  || { echo "render failed - a placeholder is not covered by the harness" >&2; exit 1; }
# Same command, wrong expected reference md5: the gate must fail before converting.
python3 "$HERE/render_command.py" "$WDL" "$WORK/task_badref.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  --set ref_fasta_md5=00000000000000000000000000000000 || exit 1

fails=0
wipe() {
  cd "$CR" || exit 73
  rm -f A.cram A.cram.crai A.cram.md5 A.cram.crai.md5
}
run() {  # run <rendered-script> -> sets RC, writes $WORK/run.out / $WORK/run.err
  wipe
  bash "$1" > "$WORK/run.out" 2> "$WORK/run.err"
  RC=$?
}
check() {  # check <label> <condition-result>
  if [ "$2" = "0" ]; then echo "  ok   $1"; else echo "  FAIL $1"; fails=$((fails+1)); fi
}

echo "== 1. plain conversion, Cromwell-shaped layout =="
run "$WORK/task.sh"
check "task rc=0 (was rc=1 with the basename() localization bug)" "$([ "$RC" -eq 0 ]; echo $?)"
check "stderr empty" "$([ ! -s "$WORK/run.err" ]; echo $?)"
check "all four deliverables present and non-empty" \
  "$([ -s "$CR/A.cram" ] && [ -s "$CR/A.cram.crai" ] && [ -s "$CR/A.cram.md5" ] \
      && [ -s "$CR/A.cram.crai.md5" ]; echo $?)"
check "delivered cram md5 matches an independent recomputation" \
  "$([ "$(cat "$CR/A.cram.md5")" = "$($MD5 "$CR/A.cram" | awk '{print $1}')" ]; echo $?)"
check "delivered crai md5 matches an independent recomputation" \
  "$([ "$(cat "$CR/A.cram.crai.md5")" = "$($MD5 "$CR/A.cram.crai" | awk '{print $1}')" ]; echo $?)"
BAMN=$(samtools view -c "$CR/datarepo-fake/uuid1/A.bam")
CRAMN=$(samtools view -c "$CR/A.cram")
check "CRAM record count equals the source BAM's ($BAMN vs $CRAMN)" \
  "$([ "$BAMN" = "$CRAMN" ]; echo $?)"
check "CRAM passes quickcheck" "$(samtools quickcheck "$CR/A.cram" >/dev/null 2>&1; echo $?)"
check "reference md5 recorded on stdout" \
  "$(grep -q "reference md5=${REF_MD5}" "$WORK/run.out"; echo $?)"

echo "== 2. wrong expected reference md5 -> FATAL before any conversion =="
run "$WORK/task_badref.sh"
check "task rc!=0" "$([ "$RC" -ne 0 ]; echo $?)"
check "stderr carries the FATAL, not an empty file" \
  "$(grep -q 'FATAL: reference md5' "$WORK/run.err"; echo $?)"
check "no CRAM was started (the gate runs before the convert)" "$([ ! -e "$CR/A.cram" ]; echo $?)"

echo "== 3. the output contract, and the units Cromwell can parse (static) =="
# Both of these are lints, not behaviours, and they exist because two canaries died in
# output evaluation *after* a successful two-hour conversion. miniwdl accepts everything
# here: `size(x, "Bytes")` parsed and checked clean on submission cf4cb95f, and Cromwell
# answered `Bad output 'ConvertBamToCram.cram_bytes': Unit with suffix Bytes was not
# found.` Cromwell's vocabulary is B/KB/MB/GB/TB/PB + KiB/MiB/GiB/TiB/PiB. Same for the
# output list: every extra declared output is another expression that has to evaluate, and
# another required file that can strand the CRAM behind it (0c6ccba5).
TASK_OUT=$(awk '/^task ConvertBamToCram/,0' "$WDL" | awk '/^  output \{/,/^  \}/')
check "task declares exactly the 4 deliverables as outputs" \
  "$([ "$(printf '%s\n' "$TASK_OUT" | grep -c 'File?[[:space:]]*output_cram')" -eq 4 ]; echo $?)"
check "...and nothing else (no reports, no derived Int/String outputs)" \
  "$([ "$(printf '%s\n' "$TASK_OUT" | grep -cE '^[[:space:]]*(Int|String|Boolean|Float|File[^?])')" -eq 0 ]; echo $?)"
check "no read_string() in an output expression (a skipped write would fail output eval)" \
  "$(! grep -vE '^[[:space:]]*#' "$WDL" | grep -qE '=\s*read_string\s*\('; echo $?)"
BAD_UNITS=$(grep -vE '^[[:space:]]*(#|##)' "$WDL" | grep -E 'size\(|memory:' | grep -oE '"[A-Za-z]{1,6}"' | tr -d '"' | sort -u \
  | grep -vE '^(B|KB|MB|GB|TB|PB|EB|KiB|MiB|GiB|TiB|PiB|EiB)$' || true)
if [ -n "$BAD_UNITS" ]; then echo "    Cromwell cannot parse these unit strings: $BAD_UNITS"; fi
check "every size()/memory unit is one Cromwell knows (cf4cb95f died on 'Bytes')" \
  "$([ -z "$BAD_UNITS" ]; echo $?)"
check "no basename() path construction in the command (58698ceb died on it)" \
  "$(! awk '/command <<</,/^  >>>/' "$WDL" | grep -vE '^[[:space:]]*#' | grep -qE '(^|[^[:alnum:]_])basename'; echo $?)"

echo
if [ "$fails" -ne 0 ]; then echo "REPLAY FAILED ($fails check(s))"; exit 1; fi
echo "REPLAY OK - bam_to_cram.wdl verified in a Cromwell-shaped tree"
