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
# The first Terra canary that got through it (submission 0c6ccba5-1045-4479-be9e-0b75548801f7,
# 2 h 01 m, docs/progress/073) then exposed three more, which cases 4-7 below guard:
#   4. an empty sampled window counted as a round-trip FAILURE, so one off-target window on
#      capture data aborted a clean conversion (23 of 24 windows matched);
#   5. a strict-gate `exit 1` that skipped writing an output file, which made Cromwell stop
#      delocalizing at the first missing required output and throw away the CRAM behind it -
#      outputs are optional now, and the completion manifest is what enforces completeness;
#   6. roundtrip.txt beginning with a stale "ROUNDTRIP=SKIPPED" line even when the round trip
#      ran, i.e. a deliverable report whose first line contradicted its last.
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

# One read every 100 bp so both sampled windows are non-empty: the tested path needs reads
# in the window, and case 4 needs to be able to make "too few tested windows" fail.
#
# The reads carry SEVEN synthetic tags each. Canary #2 measured the delivered BAM at 7.43 tag
# fields per read (299,010 tags over 40,220 reads, docs/progress/073), while this fixture was
# originally tagless - and a tagless fixture cannot fail the tag-multiset comparison, so the
# one gate that the MD/NM exclusion protects was covered only on Terra. The tag NAMES here are
# invented (plausible DRAGEN-ish flowtag shapes, not a claim about the delivered set); what is
# being reproduced is the shape: a tag-rich source, one of whose tags (NM:i:) is also a key
# htslib re-derives on decode, which the comparison must exclude by name on both sides.
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
# Same command block, different inputs: too few tested windows for the coverage floor, the
# same condition with strict_roundtrip=false, and the round trip switched off altogether.
python3 "$HERE/render_command.py" "$WDL" "$WORK/task_min8.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  --set min_roundtrip_windows=8 || exit 1
python3 "$HERE/render_command.py" "$WDL" "$WORK/task_min8_lenient.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  --set min_roundtrip_windows=8 --set strict_roundtrip=false || exit 1
python3 "$HERE/render_command.py" "$WDL" "$WORK/task_nort.sh" --root "$CR" --ref-md5 "$REF_MD5" \
  --set run_roundtrip=false || exit 1

fails=0
wipe() {
  cd "$CR" || exit 73
  rm -rf ref A.bam mini.fasta mini.fasta.fai integrity.txt roundtrip.txt roundtrip.err \
         roundtrip_identical.txt windows.txt m5_check.txt bam_idxstats.txt \
         input_bam_md5_computed.txt A.cram A.cram.crai A.cram.md5 A.cram.crai.md5 cram_sq.tsv
}
run() {  # run <m5-table-or-empty> [rendered-script] -> sets RC, writes $WORK/run.out / $WORK/run.err
  local m5="$1" src="${2:-$WORK/task.sh}"
  wipe
  sed "s|M5_TABLE=\"\"|M5_TABLE=\"$m5\"|" "$src" > "$WORK/step.sh"
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
  "$(grep -q 'ROUNDTRIP=PASS(2 tested of 2 sampled, 0 empty)' "$CR/roundtrip.txt"; echo $?)"
check "roundtrip.txt does not open with a stale ROUNDTRIP=SKIPPED line" \
  "$(head -1 "$CR/roundtrip.txt" | grep -q '^window '; echo $?)"
check "integrity.txt separates tested / empty / failing windows" \
  "$(grep -q 'roundtrip: windows=2 tested=2 empty=0 failing=0' "$CR/integrity.txt"; echo $?)"
check "a tag-rich source (7 tags/read as delivered) still compares with no tags=MISMATCH" \
  "$(grep -q 'source_tags=[1-9][0-9]* cram_derived_tags=[1-9][0-9]*' "$CR/roundtrip.txt" && ! grep -q 'tags=MISMATCH' "$CR/roundtrip.txt"; echo $?)"
check "CRAM decode re-added MD/NM on top of the source tags" \
  "$(grep -q 'cram_derived_tags=[1-9][0-9]*' "$CR/roundtrip.txt"; echo $?)"
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

echo "== 4. too few tested windows (min_roundtrip_windows=8, fixture has 2) =="
run "" "$WORK/task_min8.sh"
check "task rc!=0 - a round trip that compared too little is not a pass" "$([ "$RC" -ne 0 ]; echo $?)"
check "verdict says INSUFFICIENT_COVERAGE, with the counts" \
  "$(grep -q 'ROUNDTRIP=INSUFFICIENT_COVERAGE (tested=2 < min_roundtrip_windows=8, sampled=2, empty=0)' "$CR/roundtrip.txt"; echo $?)"
check "roundtrip_identical.txt still written (false) before the strict exit" \
  "$(grep -qx false "$CR/roundtrip_identical.txt"; echo $?)"
check "the CRAM and its index survived to be delocalized" \
  "$([ -s "$CR/A.cram" ] && [ -s "$CR/A.cram.crai" ] && [ -s "$CR/integrity.txt" ]; echo $?)"
check "stderr names the verdict, not a bare 'wrong reference?' guess" \
  "$(grep -q 'FATAL: CRAM does not read back as the source BAM over 2 tested window' "$WORK/run.err"; echo $?)"

echo "== 5. same condition, strict_roundtrip=false -> report it, do not abort =="
run "" "$WORK/task_min8_lenient.sh"
check "task rc=0 (non-strict records the disagreement instead of killing the deliverable)" "$([ "$RC" -eq 0 ]; echo $?)"
check "verdict still in roundtrip.txt" \
  "$(grep -q 'ROUNDTRIP=INSUFFICIENT_COVERAGE' "$CR/roundtrip.txt"; echo $?)"
check "roundtrip_identical=false is the machine-readable signal" \
  "$(grep -qx false "$CR/roundtrip_identical.txt"; echo $?)"

echo "== 6. run_roundtrip=false -> skipped is stated, and is not 'identical' =="
run "" "$WORK/task_nort.sh"
check "task rc=0" "$([ "$RC" -eq 0 ]; echo $?)"
check "roundtrip.txt says SKIPPED" \
  "$(grep -q 'ROUNDTRIP=SKIPPED' "$CR/roundtrip.txt"; echo $?)"
check "roundtrip_identical.txt says not_tested, so 'not tested' is neither true nor false" \
  "$(grep -qx not_tested "$CR/roundtrip_identical.txt"; echo $?)"
check "completion manifest still satisfied with the round trip off" \
  "$([ -s "$CR/A.cram" ] && [ -s "$CR/cram_sq.tsv" ]; echo $?)"

echo "== 7. completion manifest alone: teeth, and no false alarm =="
run ""
awk '/# ---- completion manifest/,0' "$WORK/task.sh" > "$WORK/manifest.sh"
check "manifest block extracted from the rendered command" "$([ -s "$WORK/manifest.sh" ]; echo $?)"
(
  cd "$CR" || exit 73
  OUT=A PATH="$BIN:$PATH" bash "$WORK/manifest.sh" > "$WORK/man.out" 2> "$WORK/man.err"
  MRC=$?
  check "complete artifact set -> manifest silent and rc=0" "$([ "$MRC" -eq 0 ]; echo $?)"
  rm -f A.cram.crai.md5
  OUT=A PATH="$BIN:$PATH" bash "$WORK/manifest.sh" > "$WORK/man.out" 2> "$WORK/man.err"
  MRC=$?
  check "one artifact removed -> rc!=0" "$([ "$MRC" -ne 0 ]; echo $?)"
  check "and it names the missing file" \
    "$(grep -q 'task finished without producing: A.cram.crai.md5' "$WORK/man.err"; echo $?)"
)

echo "== 8. the output contract, and the units Cromwell can parse (static) =="
#
# Both of these are lints, not behaviours, and they exist because two canaries died in
# output evaluation *after* a successful two-hour conversion. miniwdl accepts everything
# here: `size(x, "Bytes")` parsed and checked clean on submission cf4cb95f, and Cromwell
# answered `Bad output 'ConvertBamToCram.cram_bytes': Unit with suffix Bytes was not
# found.` Cromwell's vocabulary is B/KB/MB/GB/TB/PB + KiB/MiB/GiB/TiB/PiB - "Bytes" is not
# in it. Same for the output list: every extra declared output is another expression that
# has to evaluate, and another required file that can strand the CRAM behind it (0c6ccba5).
TASK_OUT=$(awk '/^task ConvertBamToCram/,0' "$WDL" | awk '/^  output \{/,/^  \}/')
check "task declares exactly the 4 deliverables as outputs" \
  "$([ "$(printf '%s\n' "$TASK_OUT" | grep -c 'File?[[:space:]]*output_cram')" -eq 4 ]; echo $?)"
check "...and nothing else (no reports, no derived Int/String outputs)" \
  "$([ "$(printf '%s\n' "$TASK_OUT" | grep -cE '^[[:space:]]*(Int|String|Boolean|Float|File[^?])')" -eq 0 ]; echo $?)"
check "no read_string() in an output expression (a skipped write would fail output eval)" \
  "$(! grep -vE '^[[:space:]]*#' "$WDL" | grep -qE '=\s*read_string\s*\('; echo $?)"
check "Picard ValidateCram is gone, not merely unused" \
  "$(! grep -qE 'call ValidateCram|^task ValidateCram' "$WDL"; echo $?)"
BAD_UNITS=$(grep -vE '^[[:space:]]*(#|##)' "$WDL" | grep -E 'size\(|memory:' | grep -oE '"[A-Za-z]{1,6}"' | tr -d '"' | sort -u \
  | grep -vE '^(B|KB|MB|GB|TB|PB|EB|KiB|MiB|GiB|TiB|PiB|EiB)$' || true)
if [ -n "$BAD_UNITS" ]; then echo "    Cromwell cannot parse these unit strings: $BAD_UNITS"; fi
check "every size()/memory unit is one Cromwell knows (cf4cb95f died on 'Bytes')" \
  "$([ -z "$BAD_UNITS" ]; echo $?)"
check "QC evidence is printed to stdout, so trimming outputs loses nothing" \
  "$(grep -q 'cat integrity.txt' "$WDL" && grep -q 'cat roundtrip.txt' "$WDL"; echo $?)"
check "the round-trip gate is still inside the command (trim outputs, not checks)" \
  "$(grep -q 'ROUNDTRIP=PASS' "$WDL" && grep -q 'min_roundtrip_windows' "$WDL"; echo $?)"

echo
if [ "$fails" -ne 0 ]; then echo "REPLAY FAILED ($fails check(s))"; exit 1; fi
echo "REPLAY OK - every gate in bam_to_cram.wdl exercised in a Cromwell-shaped tree"