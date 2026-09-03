version 1.0

## Single-sample GenotypeGVCFs — MOC8 long-read concordance proxy.
##
## Joint WGS calling is not ready, so genotype each benchmark sample's reblocked
## gVCF on its own (no VQSR) as a stand-in query callset for the LR concordance
## (vs HGSVC3 truth within the dipcall confident regions). Restricted to each
## sample's dipcall .dip.bed so output is confined to the scored regions.
##
## Reference MUST be the same GRCh38 build the gVCF was called against (WARP
## Homo_sapiens_assembly38, full with ALT/decoy) — GATK requires a matching dict.

workflow GenotypeGVCFsProxy {
  input {
    Array[String] sample_ids
    Array[File]   gvcfs
    Array[File]   gvcf_indices
    Array[File]   intervals_bed       # per-sample dipcall .dip.bed
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  scatter (i in range(length(gvcfs))) {
    call GenotypeGVCFs {
      input:
        sample_id   = sample_ids[i],
        gvcf        = gvcfs[i],
        gvcf_index  = gvcf_indices[i],
        intervals   = intervals_bed[i],
        ref_fasta   = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict    = ref_dict,
        gatk_docker = gatk_docker
    }
  }

  output {
    Array[File] genotyped_vcfs        = GenotypeGVCFs.vcf
    Array[File] genotyped_vcf_indices = GenotypeGVCFs.vcf_index
  }
}

task GenotypeGVCFs {
  input {
    String sample_id
    File gvcf
    File gvcf_index
    File intervals
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String gatk_docker
  }
  Int disk_gb = ceil(size(gvcf, "GB") + size(ref_fasta, "GB")) + 40

  command <<<
    set -euo pipefail
    gatk --java-options "-Xmx10g" GenotypeGVCFs \
      -R ~{ref_fasta} \
      -V ~{gvcf} \
      -L ~{intervals} \
      -O ~{sample_id}.geno.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    memory: "13 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 1
  }

  output {
    File vcf       = "~{sample_id}.geno.vcf.gz"
    File vcf_index = "~{sample_id}.geno.vcf.gz.tbi"
  }
}
