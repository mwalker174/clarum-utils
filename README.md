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

### `strip_sequence_dictionary.wdl`
Removes a stale `@SQ` sequence dictionary from an **unmapped** BAM header. The
wilkinshaug flowcell uBAMs carry a full `@SQ` (3366 contigs) even though every
read is unmapped; a proper uBAM has none. That dictionary makes
`WholeGenomeGermlineSingleSample` fail in `SamToFastqAndBwaMemAndMba` at Picard
`MergeBamAlignment`:

    Do not use this function to merge dictionaries with different sequences...
    Found [] and [chr1, chr2, ...]

Header-only rewrite via `samtools reheader` (drop `@SQ`) — the reads are
already unmapped (`RNAME=*`), so this is safe and lets `MergeBamAlignment` take
its dictionary from the reference. No realignment, no revert, no reference
required; runs in seconds-to-minutes. The task aborts if any read is mapped.

**Inputs:** `input_bam` (an unmapped BAM with a stale `@SQ`).
**Outputs:** `unmapped_bam` (`<basename>.nodict.unmapped.bam`) → wire into
`WholeGenomeGermlineSingleSample.sample_and_unmapped_bams.flowcell_unmapped_bams`.
