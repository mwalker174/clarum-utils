# clarum-utils

Utility workflows and scripts for CLARUM retro-WGS processing on Terra.

## Workflows

### `reheader_platform.wdl`
Fixes an invalid read-group `PL` (platform) tag in a BAM or CRAM. Picard
`ValidateSamFile` rejects `PL:NovaSeq X` because `PL` must use the SAM-spec
controlled vocabulary (`ILLUMINA`, `PACBIO`, `ONT`, ...); `NovaSeq X` is an
instrument model and belongs in `PM`. The workflow rewrites every `@RG`:

- `PL:NovaSeq X` → `PL:ILLUMINA`
- adds `PM:NovaSeq X` when no `PM` tag is present

Header-only rewrite via `samtools reheader` — no realignment, no reference
required. Format (BAM vs CRAM) is detected by content (magic bytes), not the
filename, so a mislabeled extension is handled correctly; the output and its
index are written with the matching extension. Run once per affected sample.

**Inputs:** `input_reads` (BAM or CRAM), `sample_name` (SM). Other fields are defaulted.
**Outputs:** `<sample_name>.reheadered.{bam,cram}` + matching `.bai`/`.crai`.
