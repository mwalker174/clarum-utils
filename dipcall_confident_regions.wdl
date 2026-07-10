version 1.0

## Generate per-sample dipcall confident regions (.dip.bed) + assembly-derived
## calls (.dip.vcf.gz) from HGSVC3 verkko diploid assemblies.
##
## MOC8 Test 3 (long-read concordance) needs a confident-region BED to pass to
## `rtg vcfeval --evaluation-regions`. HGSVC3 ships NO callable BED, so we derive
## one the GIAB-standard way: dipcall aligns each haplotype assembly to the
## reference (minimap2 -x asm5), calls variants, and emits the diploid-covered
## callable BED. See docs/progress/020-021 and the lr-concordance memory.
##
## The reference and both haplotype FASTAs are fetched by URL in a separate
## Fetch task (cloud-sdk image -- the dipcall biocontainer has no curl) and
## handed to Dipcall as File inputs, so nothing large has to be staged into GCS
## first. Requires outbound internet on the VM (fails fast in Fetch if the
## workspace is egress-locked; then stage the assemblies to gs:// instead).
##
## Reference = GRCh38 no-ALT analysis set (UCSC chr ids) -- the assembly the
## HGSVC3 truth was called against (validated 0-mismatch), coordinate-compatible
## on primary contigs with the WARP joint callset (Homo_sapiens_assembly38).
##
## The 3 benchmark samples (HGSVC/HPRC trio children, all FEMALE per 1000G ped):
##   HG00514 (CHS) -> 20241001_verkko_HG00514_fix (.v2 -- the corrected asm)
##   HG00733 (PUR) -> 20230818_verkko_batch1
##   NA19240 (YRI) -> 20240201_verkko_batch3
## HG00514 MUST use the fix assembly (the earlier one is data-corrupted).
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

    # HPRC's dipcall image -- ships the dipcall.kit with its own tested
    # minimap2/samtools/k8/htsbox at /opt/dipcall/dipcall.kit. (The biocontainers
    # conda build mangled reference contig names in the @SQ header.)
    String dipcall_docker = "humanpangenomics/hpp_dipcall_v0.3:latest"
    String fetch_docker   = "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"
    Int cpu = 8
    Int mem_gb = 32
    Int disk_gb = 200
    Int preemptible = 1
  }

  call Fetch as FetchRef {
    input: url = ref_url, out_name = "ref.fna.gz",
           docker = fetch_docker, disk_gb = 30, preemptible = preemptible
  }

  scatter (i in range(length(sample_ids))) {
    call Fetch as FetchHap1 {
      input: url = hap1_urls[i], out_name = sample_ids[i] + ".hap1.fa.gz",
             docker = fetch_docker, disk_gb = 30, preemptible = preemptible
    }
    call Fetch as FetchHap2 {
      input: url = hap2_urls[i], out_name = sample_ids[i] + ".hap2.fa.gz",
             docker = fetch_docker, disk_gb = 30, preemptible = preemptible
    }
    call Dipcall {
      input:
        sample_id   = sample_ids[i],
        ref_gz      = FetchRef.out,
        hap1_gz     = FetchHap1.out,
        hap2_gz     = FetchHap2.out,
        is_male     = is_male[i],
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

task Fetch {
  input {
    String url
    String out_name
    String docker
    Int disk_gb
    Int preemptible
  }
  command <<<
    set -euo pipefail
    curl -fSL --retry 5 --retry-delay 10 -o "~{out_name}" "~{url}"
    ls -l "~{out_name}"
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
  output { File out = "~{out_name}" }
}

task Dipcall {
  input {
    String sample_id
    File ref_gz
    File hap1_gz
    File hap2_gz
    Boolean is_male
    String docker
    Int cpu
    Int mem_gb
    Int disk_gb
    Int preemptible
  }

    KIT=/opt/dipcall/dipcall.kit

    # reference: decompress + index. Strip deflines to the bare contig name --
    # the NCBI analysis-set deflines carry trailing "AC:/gi:/LN:/rl:/M5:" fields;
    # bare ">chrN" keeps the @SQ header clean regardless of the minimap2 build.
    # run-dipcall requires ref.fna.fai to exist before it generates the makefile.
    gzip -dc ~{ref_gz} | awk '/^>/{print $1; next} {print}' > ref.fna
    $KIT/samtools faidx ref.fna

    # dipcall's kit reads plain (uncompressed) assembly FASTAs
    gzip -dc ~{hap1_gz} > hap1.fa
    gzip -dc ~{hap2_gz} > hap2.fa

    PAR_OPT=""
    if [ "~{is_male}" = "true" ]; then PAR_OPT="-x $KIT/hs38.PAR.bed"; fi

    # run-dipcall emits a Makefile; make runs minimap2 -> filter -> sort ->
    # htsbox pileup -> vcfpair, and derives the diploid callable BED via bedtk.
    $KIT/run-dipcall -t ~{cpu} $PAR_OPT ~{sample_id} ref.fna hap1.fa hap2.fa > ~{sample_id}.mak
    make -j2 -f ~{sample_id}.mak

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
