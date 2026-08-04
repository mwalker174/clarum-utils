version 1.0

## Remove duplicate/redundant samples from a joint-called VCF (SNV/indel or
## SV). Scatters over one or more VCF shards so it works equally for a single
## joint VCF or a sharded callset (per-contig SV VCF, per-interval SNV/indel
## VCF) -- pass a 1-element array for the unsharded case.
##
## Sample-removal logic follows gatk-sv's SubsetVcfBySamplesList task
## (github.com/broadinstitute/gatk-sv/blob/7c3ea4930d90fcbae9dc8b750a5b3a20923af5dc/wdl/Utils.wdl#L852):
## `bcftools view -S ^<list> --force-samples` to drop the listed samples, then
## a second pass to drop sites left with no non-ref genotype among the
## survivors ("private" to the removed samples). For SV VCFs, CNV (depth-only)
## records are exempted from that private-site drop, same as upstream --
## COUNT(GT="alt") is not a meaningful "any carrier" test for those records.
##
## exclude_samples_list: one sample ID per line, no header. This is a fixed
## exclude list, not a duplicate-detection step -- decide which member of
## each duplicate cluster to drop upstream (see docs/progress/049) and pass
## the resulting IDs here.

workflow RemoveDuplicateSamples {
  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File exclude_samples_list
    Boolean is_sv_vcf = false          # true: exempt CNV records from the private-site drop
    Boolean recalculate_af = true      # true (default): let bcftools recompute AC/AN/AF after removal
    Boolean remove_private_sites = true
    String bcftools_docker = "staphb/bcftools:1.19"
  }

  scatter (i in range(length(vcfs))) {
    call RemoveSamplesFromVcf {
      input:
        vcf                   = vcfs[i],
        vcf_idx               = vcf_indices[i],
        exclude_samples_list  = exclude_samples_list,
        is_sv_vcf             = is_sv_vcf,
        recalculate_af        = recalculate_af,
        remove_private_sites  = remove_private_sites,
        bcftools_docker       = bcftools_docker
    }
  }

  output {
    Array[File] dedup_vcfs         = RemoveSamplesFromVcf.dedup_vcf
    Array[File] dedup_vcf_indices  = RemoveSamplesFromVcf.dedup_vcf_index
  }
}

task RemoveSamplesFromVcf {
  input {
    File vcf
    File vcf_idx
    File exclude_samples_list
    Boolean is_sv_vcf
    Boolean recalculate_af
    Boolean remove_private_sites
    String? output_basename
    String bcftools_docker
  }

  String basename_out = select_first([output_basename, basename(basename(vcf, ".gz"), ".vcf") + ".dedup"])
  String out_vcf = basename_out + ".vcf.gz"

  # Two passes over the VCF (remove samples, then drop now-private sites) -> ~2x input size, plus output.
  Int disk_gb = ceil(size(vcf, "GB") * 3) + 20

  String private_site_expr = if is_sv_vcf then "SVTYPE!=\"CNV\" && COUNT(GT=\"alt\")==0" else "COUNT(GT=\"alt\")==0"

  command <<<
    set -euo pipefail

    if ~{if remove_private_sites then "true" else "false"}; then
      bcftools view \
        -S ^~{exclude_samples_list} \
        --force-samples \
        ~{if recalculate_af then "" else "--no-update"} \
        -Ou \
        "~{vcf}" \
      | bcftools view -e '~{private_site_expr}' -Oz -o "~{out_vcf}"
    else
      bcftools view \
        -S ^~{exclude_samples_list} \
        --force-samples \
        ~{if recalculate_af then "" else "--no-update"} \
        -Oz -o "~{out_vcf}" \
        "~{vcf}"
    fi

    bcftools index -t "~{out_vcf}"
  >>>

  runtime {
    docker: bcftools_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }

  output {
    File dedup_vcf = out_vcf
    File dedup_vcf_index = out_vcf + ".tbi"
  }
}
