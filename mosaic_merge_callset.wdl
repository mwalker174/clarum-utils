version 1.0

## Merge per-sample Mutect2 filtered VCFs into ONE multi-sample joint VCF —
## the MOC10 / TDD 2.J.5 "1 joint VCF" mosaic deliverable.
##
## The mosaic calls themselves come from the canonical Broad somatic pipeline
## (github.com/broadinstitute/gatk scripts/mutect2_wdl):
##   Mutect2_Panel (mutect2_pon.wdl)  -> kit/convention-matched PoN
##   Mutect2       (mutect2.wdl)      -> per-sample tumor-only calling + filtering
## Mutect2 is per-sample and has no joint-genotyping step (somatic calling never
## does), so we gather the per-sample FilterMutectCalls outputs here with
## `bcftools merge` to produce the single multi-sample callset the TDD asks for.
##
## Inputs are the `output_vcf` / `output_vcf_idx` produced by per-sample Mutect2
## runs (e.g. collected from a Terra data-table column across the sample set).

workflow MosaicMergeCallset {
  input {
    Array[File] sample_vcfs
    Array[File] sample_vcf_indices
    String callset_name = "mosaic_callset"
    String bcftools_docker = "staphb/bcftools:1.19"
  }

  call MergeVcfs {
    input:
      sample_vcfs = sample_vcfs,
      sample_vcf_indices = sample_vcf_indices,
      callset_name = callset_name,
      bcftools_docker = bcftools_docker
  }

  output {
    File joint_vcf = MergeVcfs.joint_vcf
    File joint_vcf_index = MergeVcfs.joint_vcf_index
  }
}

task MergeVcfs {
  input {
    Array[File] sample_vcfs
    Array[File] sample_vcf_indices
    String callset_name
    String bcftools_docker
  }
  Int disk_gb = ceil(size(sample_vcfs, "GB")) + 20

  command <<<
    set -euo pipefail
    # FilterMutectCalls emits plain-text VCF + Picard .idx, not bgzip+tabix —
    # bcftools merge requires compressed+tabixed input, so index each first.
    mkdir -p compressed
    i=0
    for f in ~{sep=' ' sample_vcfs}; do
      i=$((i+1))
      # Drop AS_FilterStatus: known GATK bug, Number=A count can mismatch
      # ALT count at multiallelic sites, which htslib's strict parser rejects.
      # Not needed for the merge / IGV accuracy check.
      awk 'BEGIN{FS=OFS="\t"} /^#/{print; next}
           { n=split($8,a,";"); out="";
             for (j=1;j<=n;j++) if (a[j] !~ /^AS_FilterStatus=/)
               out = (out=="") ? a[j] : out";"a[j];
             $8=out; print }' "$f" > compressed/sample_${i}.vcf
      bcftools view -Oz -o compressed/sample_${i}.vcf.gz compressed/sample_${i}.vcf
      bcftools index -t compressed/sample_${i}.vcf.gz
      rm -f compressed/sample_${i}.vcf
    done
    # Per-sample single-sample VCFs -> one multi-sample joint VCF.
    bcftools merge \
      --output-type z \
      --output ~{callset_name}.vcf.gz \
      compressed/*.vcf.gz
    bcftools index --tbi ~{callset_name}.vcf.gz
  >>>

  runtime {
    docker: bcftools_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: 1
  }

  output {
    File joint_vcf = "~{callset_name}.vcf.gz"
    File joint_vcf_index = "~{callset_name}.vcf.gz.tbi"
  }
}
