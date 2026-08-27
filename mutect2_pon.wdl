version 1.0

## Fork of broadinstitute/gatk/mutect2_pon@4.6.1.0 (scripts/mutect2_wdl/mutect2_pon.wdl).
## Upstream source: https://github.com/broadinstitute/gatk/blob/4.6.1.0/scripts/mutect2_wdl/mutect2_pon.wdl
##
## Reason for fork: CreatePanel's GenomicsDBImport task hard-codes machine_mem = 8 (GB)
## with no workflow input to override it. That's fine for small panels (tested at 110/50
## normals) but OOMs on a 1383-sample panel (CLARUM v3 PoN, all 70 CreatePanel shards
## exhausted 3 retries at 8GB with memory_retry_multiplier unset). Added create_panel_mem_gb
## and create_panel_disk_gb workflow inputs (default to upstream's original values, so this
## fork is a no-op unless you set them) so memory/disk can be sized to the panel without
## editing the WDL again. Also switched the `import "mutect2.wdl"` to an absolute GitHub
## URL pinned at the same 4.6.1.0 tag, since this file now lives standalone outside the
## gatk repo's scripts/mutect2_wdl/ directory.
##
## Everything else (Mutect2/SplitIntervals/MergeVCFs calls, CreatePanel's command block)
## is unchanged from upstream.
##
#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  normal_bams: arrays of normal bams
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf
#  m2_extra_args: additional command line parameters for Mutect2.  This should not involve --max-mnp-distance,
#  which the wdl hard-codes to 0 because GenomicsDBImport can't handle MNPs
#  create_panel_mem_gb: memory (GB) for CreatePanel's GenomicsDBImport + CreateSomaticPanelOfNormals; bump for large panels
#  create_panel_disk_gb: local disk (GB) for CreatePanel; bump for large panels

import "https://raw.githubusercontent.com/broadinstitute/gatk/4.6.1.0/scripts/mutect2_wdl/mutect2.wdl" as m2
workflow Mutect2_Panel {
  input {
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    Int scatter_count
    Array[String] normal_bams
    Array[String] normal_bais
    File gnomad
    File gnomad_idx
    String? m2_extra_args
    String? create_pon_extra_args
    Boolean compress = false
    String pon_name

    Int min_contig_size = 1000000
    Int create_panel_scatter_count = 24
    Int create_panel_mem_gb = 8
    Int create_panel_disk_gb = 100

    String? gcs_project_for_requester_pays

    # runtime
    String gatk_docker
    File? gatk_override
    String basic_bash_docker = "ubuntu:16.04"

    Int preemptible = 2
    Int max_retries = 2
    Int small_task_cpu = 2
    Int small_task_mem = 4
    Int small_task_disk = 100
    Int boot_disk_size = 12
  }

  Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries, "preemptible": preemptible, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk, "boot_disk_size": boot_disk_size}

    scatter (normal_bam in zip(normal_bams, normal_bais)) {
        call m2.Mutect2 {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = normal_bam.left,
                tumor_reads_index = normal_bam.right,
                scatter_count = scatter_count,
                m2_extra_args = select_first([m2_extra_args, ""]) + "--max-mnp-distance 0",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    call m2.SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = create_panel_scatter_count,
            split_intervals_extra_args = "--dont-mix-contigs --min-contig-size " + min_contig_size,
            runtime_params = standard_runtime
    }

    scatter (subintervals in SplitIntervals.interval_files ) {
            call CreatePanel {
                input:
                    input_vcfs = Mutect2.filtered_vcf,
                    intervals = subintervals,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    gnomad = gnomad,
                    gnomad_idx = gnomad_idx,
                    output_vcf_name = pon_name,
                    create_pon_extra_args = create_pon_extra_args,
                    machine_mem = create_panel_mem_gb,
                    disk_gb = create_panel_disk_gb,
                    runtime_params = standard_runtime
            }
    }

    call m2.MergeVCFs {
        input:
            input_vcfs = CreatePanel.output_vcf,
            input_vcf_indices = CreatePanel.output_vcf_index,
            compress_vcfs = compress,
            runtime_params = standard_runtime
    }

    output {
        File pon = MergeVCFs.merged_vcf
        File pon_idx = MergeVCFs.merged_vcf_idx
        Array[File] normal_calls = Mutect2.filtered_vcf
        Array[File] normal_calls_idx = Mutect2.filtered_vcf_idx
    }
}

task CreatePanel {
    input {
      File intervals
      Array[String] input_vcfs
      File ref_fasta
      File ref_fai
      File ref_dict
      String output_vcf_name
      File gnomad
      File gnomad_idx
      String? create_pon_extra_args

      Int machine_mem = 8
      Int disk_gb = 100

      # runtime
      Runtime runtime_params
    }

    Int command_mem = machine_mem - 1

        parameter_meta{
          gnomad: {localization_optional: true}
          gnomad_idx: {localization_optional: true}
        }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk GenomicsDBImport --genomicsdb-workspace-path pon_db -R ~{ref_fasta} -V ~{sep=' -V ' input_vcfs} -L ~{intervals}

        gatk --java-options "-Xmx~{command_mem}g"  CreateSomaticPanelOfNormals -R ~{ref_fasta} --germline-resource ~{gnomad} \
            -V gendb://pon_db -O ~{output_vcf_name}.vcf ~{create_pon_extra_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File output_vcf = "~{output_vcf_name}.vcf"
        File output_vcf_index = "~{output_vcf_name}.vcf.idx"
    }
}
