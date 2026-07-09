version 1.0

## Generate per-sample dipcall confident regions (.dip.bed) + assembly-derived
## calls (.dip.vcf.gz) from HGSVC3 verkko diploid assemblies.
##
## MOC8 Test 3 (long-read concordance) needs a confident-region BED to pass to
## `rtg vcfeval --evaluation-regions`. HGSVC3 ships NO callable BED, so we derive
## one the GIAB-standard way: dipcall aligns each haplotype assembly to the
## reference (minimap2 -x asm5), calls variants, and emits the diploid-covered
## callable BED. See docs/progress/020 and the lr-concordance memory.
##
## Self-contained: the reference and both haplotype FASTAs are fetched by URL
## inside the task (Terra VMs have open internet), so nothing large has to be
## staged into GCS first. Reference = GRCh38 no-ALT analysis set (UCSC chr ids) --
## the same assembly the HGSVC3 truth was called against (validated 0-mismatch),
## and coordinate-compatible on primary contigs with the WARP joint callset
## (Homo_sapiens_assembly38), so the BED drops straight into vcfeval.
##
## The 3 benchmark samples (HGSVC/HPRC trio children) and their verkko assemblies:
##   HG00514 (CHS, male)  -> 20241001_verkko_HG00514_fix (.v2 -- the corrected asm)
##   HG00733 (PUR, male)  -> 20230818_verkko_batch1
##   NA19240 (YRI, female)-> 20240201_verkko_batch3
## Male samples pass -x <PAR.bed> so chrX PAR is treated diploid and non-PAR
## chrX/chrY haploid.
##
## Outputs (one per sample): <sample>.dip.bed, <sample>.dip.vcf.gz -> feed as
##   run_concordance.py --eval-regions-map HG00514=HG00514.dip.bed ...

workflow DipcallConfidentRegions {
  input {
    Array[String] sample_ids
    Array[String] hap1_urls
    Array[String] hap2_urls
    Array[Boolean] is_male

    # GRCh38 no-ALT analysis set (UCSC chr ids), gzipped. Matches the HGSVC3 truth.
    String ref_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

    String dipcall_docker = "quay.io/biocontainers/dipcall:0.3--hdfd78af_0"
    Int cpu = 8
    Int mem_gb = 32
    Int disk_gb = 200
    Int preemptible = 1
  }

  scatter (i in range(length(sample_ids))) {
    call Dipcall {
      input:
        sample_id   = sample_ids[i],
        hap1_url    = hap1_urls[i],
        hap2_url    = hap2_urls[i],
        is_male     = is_male[i],
        ref_url     = ref_url,
        docker      = dipcall_docker,
        cpu         = cpu,
        mem_gb      = mem_gb,
        disk_gb     = disk_gb,
        preemptible = preemptible
    }
  }

  output {
    Array[File] dip_bed = Dipcall.dip_bed
    Array[File] dip_vcf = Dipcall.dip_vcf
  }

  meta { description: "Per-sample dipcall confident-region BED + calls from HGSVC3 verkko assemblies for MOC8 long-read concordance." }
}

task Dipcall {
  input {
    String sample_id
    String hap1_url
    String hap2_url
    Boolean is_male
    String ref_url
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

  command <<<
    set -euo pipefail

    # --- reference (fetch + index; samtools faidx needs an uncompressed FASTA) ---
    curl -fsSL --retry 5 --retry-delay 10 -o ref.fna.gz "~{ref_url}"
    gzip -d ref.fna.gz
    samtools faidx ref.fna

    # --- haplotype assemblies (minimap2 reads gzip directly, no faidx needed) ---
    curl -fsSL --retry 5 --retry-delay 10 -o hap1.fa.gz "~{hap1_url}"
    curl -fsSL --retry 5 --retry-delay 10 -o hap2.fa.gz "~{hap2_url}"

    # --- GRCh38 pseudoautosomal regions on chrX (0-based, chr ids) ---
    printf 'chrX\t10000\t2781479\nchrX\t155701382\t156030895\n' > hs38.PAR.bed

    PAR_OPT=""
    if [ "~{is_male}" = "true" ]; then PAR_OPT="-x hs38.PAR.bed"; fi

    # run-dipcall emits a Makefile; make runs minimap2 -> filter -> sort ->
    # htsbox pileup -> vcfpair, and derives the diploid callable BED via bedtk.
    run-dipcall -t ~{cpu} $PAR_OPT ~{sample_id} ref.fna hap1.fa.gz hap2.fa.gz > ~{sample_id}.mak
    make -j2 -f ~{sample_id}.mak

    # sanity: the two headline outputs must exist and be non-empty
    test -s ~{sample_id}.dip.bed
    test -s ~{sample_id}.dip.vcf.gz
    echo "confident-region span (bp):"
    awk '{s+=$3-$2} END{print s}' ~{sample_id}.dip.bed
  >>>

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }

  output {
    File dip_bed = "~{sample_id}.dip.bed"
    File dip_vcf = "~{sample_id}.dip.vcf.gz"
  }
}
