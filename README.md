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

### `bam_to_cram.wdl`
Converts an already-aligned BAM to CRAM, with the index, md5s and integrity evidence
the ARPA-H deposit needs. Exists because the prospective somatic ES delivery is BAM-only
(TDD §2.H.3 asks for "generation of sequencing CRAM file", §2.H.4 transfers "≥14 CRAM
files ... integrity confirmed via file size and MD5 checks", but the TDR snapshots behind
`clarum_talkowski_strub_ach_prospective_wes` / `..._smits_bch_prospective_wes` ship
`bam_file`, `bai_file`, `bam_md5_sum_path`, four DRAGEN metrics CSVs and a gVCF — no CRAM
column, and both workspaces have zero method configs).

Deliberately **not** WARP's `tasks/broad/BamToCram.wdl`: that is a sub-workflow of
`ExomeGermlineSingleSample` which requires the pipeline's own uBAM duplication/chimerism
metric maps as inputs. This is the `ConvertToCram` idea plus the checks we need.

**The reference is the whole ballgame.** CRAM pins its reference by MD5, and these BAMs
cannot tell us which hg38 DRAGEN used: `@PG` says `-r /seq/dragen/references/hg38/dragen_komodo/`
and `--ht-reference=…/Hsapiens/hg38/seq/hg38.fa`, the header carries 3,366 `@SQ` and
`M5:` on none of them. 3,366 is also the contig count of Broad's
`hg38/v0/Homo_sapiens_assembly38.fasta`, which makes the two look interchangeable — they
are not: 2,841 name+length entries match, 525 differ, all in the HLA block (`HLA-A*01`,
`HLA-A*01:01:38L`, … vs `HLA-A*01:01:01:01`, `HLA-A*01:02`, …). So `ref_fasta` is a
required, un-defaulted input: pass the FASTA DRAGEN actually used (ask BI for its md5),
and record that md5 with the deposit — contig *names* match Broad's, so joint calling is
unaffected, but only the same file unpacks these CRAMs.

Three checks, two of them gates:

- **md5 of the source BAM** vs the delivered `.md5sum` (§2.H.4's integrity evidence;
  reported as `INPUT_BAM_MD5=MATCH|MISMATCH|NOT_SUPPLIED`).
- **Round trip** over sampled windows on the longest primary contigs: records read back
  out of the CRAM must equal the source BAM's core alignment fields (QNAME…SEQ/QUAL).
  Gated (`strict_roundtrip`, default true). A whole-record `cmp` is *expected* to differ —
  htslib re-derives `MD:Z` on read-back — so the gate compares fields 1–11, sorted.
  Verified negative control: a near-miss reference gives
  `MD5 checksum reference mismatch` → empty output → the task aborts.
- **`ref_agreement`** (reported, not gated): `NM:i:` comes from DRAGEN's alignment against
  *their* reference, `MD:Z:` is re-derived against *ours*; if the two disagree the
  references differ, which is the one thing the round trip cannot see (a CRAM stores
  differences against whatever FASTA it was given, so a wrong one still round-trips).

**Inputs:** `input_bam`, `ref_fasta` + `ref_fasta_index` (required), optional
`bam_index_path` and `bam_md5_path`, `output_basename`. Knobs: `run_roundtrip` /
`strict_roundtrip` / window size+count, `run_validation` / `skip_mate_validation` /
`validation_ignore`, `cpu` (drives `samtools -@`), `mem_gb`, `additional_disk_gb`,
`preemptible`, docker overrides. One exome BAM is ~20 GiB, so disk lands near 85 GiB.
`bam_index_path` matters: BI puts the `.bai` in a *different* datarepo-row prefix than the
`.bam`, and anything assuming `<bam>.bai` beside the BAM will not find it; the task copies
it in, verifies it with `idxstats`, and re-indexes if it is unusable.

**Outputs:** `output_cram`, `output_cram_index`, `output_cram_md5`,
`output_cram_index_md5`, `integrity_report` (md5s, sizes, record counts, ROUNDTRIP +
REF_AGREEMENT lines), `roundtrip_identical` (Boolean), `validation_report`. Wire them to
`sample`-table attributes `output_cram` / `output_cram_index` / `output_cram_md5` as
before.

Validation runs Picard `ValidateSamFile` `picard-cloud:2.26.10` — the same image and
version WARP's `QC.ValidateSamFile` uses and the version the delivered 620-sample gVCFs
came out of — with `MODE=VERBOSE`, `MISSING_TAG_NM` ignored and mate validation on unless
skipped, so these reports stay comparable to the retrospective ones. No `REF_DICT` needed
(WARP declares one and never passes it).

`miniwdl check` clean. The task body was additionally run end-to-end against a synthetic
4-contig / 405-read fixture (two index paths, two md5 paths, correct-reference pass,
mutated-reference abort, `ref_agreement` PLAUSIBLE vs NO_INFORMATIVE_RECORDS) on samtools
1.23; the pinned image is samtools 1.11, the one WARP's `ConvertToCram` runs on, which is
also where `seq_cache_populate.pl` comes from. See `bam_to_cram.inputs.json` for a
canary-shaped input file (Strub `ACH0024_LM17`).
