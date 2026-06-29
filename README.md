# clarum-utils

Utility workflows and scripts for CLARUM retro-WGS processing on Terra.

## Workflows

### `reheader_platform.wdl`
Fixes an invalid read-group `PL` (platform) tag in a CRAM. Picard
`ValidateSamFile` rejects `PL:NovaSeq X` because `PL` must use the SAM-spec
controlled vocabulary (`ILLUMINA`, `PACBIO`, `ONT`, ...); `NovaSeq X` is an
instrument model and belongs in `PM`. The workflow rewrites every `@RG`:

- `PL:NovaSeq X` → `PL:ILLUMINA`
- adds `PM:NovaSeq X` when no `PM` tag is present

Header-only rewrite via `samtools reheader` — no realignment, no reference
required. Run once per affected sample.

**Inputs:** `input_cram`, `sample_name` (SM). Other fields are defaulted.
**Outputs:** `<sample_name>.reheadered.cram` + `.crai`.
