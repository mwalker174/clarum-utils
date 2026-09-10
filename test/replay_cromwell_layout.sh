#!/usr/bin/env bash
# Replay ConvertBamToCram's command block against a Cromwell-shaped file layout, without
# Cromwell, Docker, or the 21.9 GB BAM.
#
# WHY THIS EXISTS. The task "worked locally" and still failed on Terra with rc=1, an empty
# stderr, and stdout stopping after the samtools version banner. Reason: Cromwell/GCP Batch
# materializes File inputs at <call_root>/<bucket>/<object prefix>/<name> while the task
# runs *in* <call_root>, so a path built with basename() is dangling there but fine in a
# local run, which drops inputs into cwd under their own name. Nothing about `miniwdl check`
# or a local happy-path run can see that. This script rebuilds the cloud layout with a 5 kB
# synthetic reference and asserts the whole gate chain actually passes through it.
#
# It found three defects in one sitting (docs/progress/072 of the clarum workspace):
#   1. the basename() localization bug above;
#   2. the round-trip tag gate, which compared MD:Z:/NM:i: that CRAM decoding derives and the
#      source BAMs do not carry, and used `cut -f12-`, which emits a blank line per tagless
#      record - together: tags=MISMATCH, and strict_roundtrip=true aborts on that;
#   3. the reference-M5 join, which used `exp` (awk's builtin exponential) as an array name,
#      so the block died with an awk syntax error the first time it executed.
#
# Usage:  test/replay_cromwell_layout.sh [--wdl path/to/bam_to_cram.wdl]
# Requires: samtools (+ python3, awk, sed). Reference-free: the M5 table case is synthesized.
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
mkdir -p "$CR/datarepo-fake/uuid1" "$CR/datarepo-fake/uuid2" "$CR/datarepo-fake/uuid3" \
         "$CR/gcp-public-data--broad-references/hg38/v0" "$BIN"
trap 'rm -rf "$WORK"' EXIT

# macOS has no md5sum and no seq_cache_populate.pl (both ship in the samtools image):
# give the replay a faithful stand-in instead of skipping the steps that need them.
if ! command -v md5sum >/dev/null; then
  printf '#!/bin/sh\nmd5 -r "$@" | awk "{print \$1\"  \" \$2}"\n' > "$BIN/md5sum"
fi
cat > "$BIN/seq_cache_populate.pl" <<'EOS'
#!/bin/bash
set -euo pipefail
root=""; fasta=""
while [ $# -gt 0 ]; do case "$1" in -root) root="$2"; shift 2;; *) fasta="$1"; shift;; esac; done
samtools faidx "$fasta"
cut -f1 "$fasta.fai" | while read -r name; do
  d="$root/$(printf %s "$name" | cut -c1-2)/$(printf %s "$name" | cut -c3-4)"
  mkdir -p "$d"; samtools faidx "$fasta" "$name" > "$d/$name"
done
EOS
chmod +x "$BIN"/* 2>/dev/null
MD5=md5sum; command -v md5sum >/dev/null || MD5="$BIN/md5sum"

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

# One read every 100 bp so that the sampled windows are non-empty: an empty window is a
# failed test, not a passed one, and the task says so.
Q="$(printf 'I%.0s' $(seq 100))"
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\n'
  for p in $(seq 200 100 4800); do
    printf 'r%s\t99\tchr1\t%s\t60\t100M\t=\t%s\t300\t%s\t%s\n' \
      "$p" "$p" "$((p+300))" \
      "$(python3 -c "import random;random.seed($p);print(''.join(random.choice('ACGT') for _ in range(100)))")" \
      "$Q"
  done
} > "$WORK/body.sam"
{ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\n'
  grep -v '^@' "$WORK/body.sam" | sort -k3,3 -k4,4n; } > "$WORK/a.sam"
samtools sort -o "$CR/datarepo-fake/uuid1/A.bam" "$WORK/a.sam"
samtools index "$CR/datarepo-fake/uuid1/A.bam"
# The delivered layout splits these three objects across different datarepo-row prefixes,
# so nothing may assume they sit next to each other.
mv "$CR/datarepo-fake/uuid1/A.bam.bai" "$CR/datarepo-fake/uuid2/A.bam.bai"
$MD5 "$CR/datarepo-fake/uuid1/A.bam" | awk '{print $1"  A.bam"}' > "$CR/datarepo-fake/uuid3/A.bam.md5sum"

# A correct expected-M5 table, and one that disagrees, to prove the gate bites in both ways.
python3 - "$REF" "$WORK/m5_ok.tsv" "$WORK/m5_bad.tsv" <<'PY'
import hashlib, sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
m5 = hashlib.md5(seq.encode()).hexdigest()
open(sys.argv[2], "w").write("#contig\tlength\tm5\nchr1\t5000\t%s\n" % m5)
open(sys.argv[3], "w").write("#contig\tlength\tm5\nchr1\t5000\t%032d\n" % 0)
PY

python3 "$HERE/render_command.py" "$WDL" "$WORK/task.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  || { echo "render failed - a placeholder is not covered by the harness" >&2; exit 1; }

fails=0
wipe() {
  cd "$CR" || exit 73
  rm -rf ref A.bam mini.fasta mini.fasta.fai integrity.txt roundtrip.txt roundtrip.err \
         roundtrip_identical.txt windows.txt m5_check.txt bam_idxstats.txt \
         input_bam_md5_computed.txt A.cram A.cram.crai A.cram.md5 A.cram.crai.md5 cram_sq.tsv
}
run() {  # run <m5-table-or-empty> -> sets RC, writes $WORK/run.out / $WORK/run.err
  local m5="$1"
  wipe
  sed "s|M5_TABLE=\"\"|M5_TABLE=\"$m5\"|" "$WORK/task.sh" > "$WORK/step.sh"
  PATH="$BIN:$PATH" bash "$WORK/step.sh" > "$WORK/run.out" 2> "$WORK/run.err"
  RC=$?
}
check() {  # check <label> <condition-result>
  if [ "$2" = "0" ]; then echo "  ok   $1"; else echo "  FAIL $1"; fails=$((fails+1)); fi
}

echo "== 1. plain conversion (no M5 table), Cromwell-shaped layout =="
run ""
check "task rc=0 (was rc=1 with the basename() localization bug)" "$([ "$RC" -eq 0 ]; echo $?)"
check "stderr empty" "$([ ! -s "$WORK/run.err" ]; echo $?)"
check "supplied .bai found in its own prefix and used" \
  "$([ -f "$CR/A.bam.bai" ]; echo $?)"
check "supplied BAM md5 verified against the localized copy" \
  "$(grep -q 'INPUT_BAM_MD5=MATCH' "$CR/integrity.txt"; echo $?)"
check "BAM/CRAM record counts agree" \
  "$(grep -q 'RECORD_COUNT=MATCH' "$CR/integrity.txt"; echo $?)"
check "CRAM indexed (.crai) and non-empty" "$([ -s "$CR/A.cram.crai" ]; echo $?)"
check "round-trip gate PASSED on sampled windows (tags=MISMATCH bug)" \
  "$(grep -q 'ROUNDTRIP=PASS' "$CR/roundtrip.txt"; echo $?)"
check "roundtrip_identical output = true" \
  "$(grep -qx true "$CR/roundtrip_identical.txt"; echo $?)"
echo "  per-window report:"; sed -n '2,$p' "$CR/roundtrip.txt" 2>/dev/null | sed 's/^/    /'

echo "== 2. reference M5 table agrees -> gate must report MATCH (rc still 0) =="
run "$WORK/m5_ok.tsv"
check "task rc=0 (awk 'exp' array bug made this rc=2)" "$([ "$RC" -eq 0 ]; echo $?)"
check "REFERENCE_M5=MATCH ok=1" \
  "$(grep -q 'REFERENCE_M5=MATCH ok=1' "$CR/m5_check.txt" 2>/dev/null; echo $?)"

echo "== 3. reference M5 table disagrees -> must be FATAL and LOUD =="
run "$WORK/m5_bad.tsv"
check "task rc!=0" "$([ "$RC" -ne 0 ]; echo $?)"
check "stderr carries the FATAL, not an empty file" \
  "$(grep -q 'FATAL: CRAM @SQ M5 disagrees' "$WORK/run.err"; echo $?)"

echo
if [ "$fails" -ne 0 ]; then echo "REPLAY FAILED ($fails check(s))"; exit 1; fi
echo "REPLAY OK - every gate in bam_to_cram.wdl exercised in a Cromwell-shaped tree"
