import "local_ancestry_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "local_ancestry_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as BCFTOOLS

task subset_ancestry_samples{
    File sample_file_in
    Array[String] ancestries_to_include
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2

    command <<<
        # Grep anything that matches any of the ancestries
        set -e

        psam=${sample_file_in}

         # Unzip psam if necessary
        if [[ ${sample_file_in} =~ \.gz$ ]]; then
            gunzip -c ${sample_file_in} > ref.psam
            psam=ref.psam
        fi

        grep "${sep="\\|" ancestries_to_include}" $psam | tr -s ' ' '\t' > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File sample_file_out = "${output_filename}"
    }
}

task format_rfmix_genetic_map{
    # Most genetic maps have headers, are space-delimited, and do not have chr at beginning
    # RFMix expects no header, tab-delimited, and chr be the first column
    File genetic_map_in
    String chr
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2

    command <<<

        set -e

        mapfile=${genetic_map_in}

         # Unzip psam if necessary
        if [[ ${genetic_map_in} =~ \.gz$ ]]; then
            gunzip -c ${genetic_map_in} > temp.txt
            mapfile=temp.txt
        fi

        # Remove column, select
        tail -n +2 ${genetic_map_in} | cut -d" " -f1,3 | sed "s/^/${chr} /" | tr -s ' ' '\t' > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File genetic_map_out = "${output_filename}"
    }
}

workflow make_rfmix_refs_wf{
    Array[File] ref_vcfs_in
    Array[File] ref_vcfs_in_tbi
    Array[String] chrs
    String output_basename
    String output_type = "z"

    # Optional regions to subselect within each chromosome
    Array[String]? regions
    # This is really stupid WDL trickery to get around the fact that WDL will not let you index an optoinal array
    # So for instance if I were to do regions[1] it would fail.
    # But if I convert regions into a array using select_first, it can be indexed.
    Array[String] actual_regions = if(defined(regions)) then select_first([regions]) else ["NONE"]


    # Genetic map files by chr
    Array[File] ref_genetic_maps_in

    # Sample info for subsetting to correct reference samples by pop or superpop
    File sample_file_in
    # 1-based index of column containing sample id in sample file
    Int sample_file_sample_id_col
    # 1-based column index of column that will contain ancestry group (usually 2 for pop, 3 for superpop but make sure)
    Int sample_file_pop_id_col
    Array[String] ancestries_to_include = ["AFR", "EUR", "EAS"]

    # Resources used to split each chromosome
    Int bcftools_chr_cpu = 2
    Int bcftools_chr_mem_gb = 8


    # Subset samples to only desired ref populations you want to partition ancestry among
    call subset_ancestry_samples{
        input:
            sample_file_in = sample_file_in,
            ancestries_to_include = ancestries_to_include,
            output_filename = "${output_basename}.samples.txt"
    }

    # Create sample id list for subsetting vcf files
    call UTILS.cut as get_sample_ids{
        input:
             input_file = subset_ancestry_samples.sample_file_out,
             args = "-f ${sample_file_sample_id_col}",
             output_filename = "${output_basename}.samples.idsonly.txt"
    }

    # Create sample id list for subsetting vcf files
    call UTILS.cut as make_rfmix_sample_file{
        input:
             input_file = subset_ancestry_samples.sample_file_out,
             args = "-f ${sample_file_sample_id_col},${sample_file_pop_id_col}",
             output_filename = "${output_basename}.samples.rfmix.txt"
    }

    String? null_region
    scatter(chr_index in range(length(chrs))){

        String chr = chrs[chr_index]

        # This is WDL trickery because you cannot actually index an optioinal array and there is no null literal character
        # Frustrating and weird but whatevs
        String? region = if(defined(regions)) then actual_regions[chr_index] else null_region

        # Format genetic map file
        call format_rfmix_genetic_map{
            input:
                genetic_map_in = ref_genetic_maps_in[chr_index],
                chr = chr,
                output_filename = "${output_basename}.chr${chr}.rfmix_genetic_map.txt"
        }

        # Subset chr vcf to include only samples of interest and possible even just regions of interest
        call BCFTOOLS.view as subset_ref_vcf{
            input:
                vcf_in = ref_vcfs_in[chr_index],
                vcf_tbi_in = ref_vcfs_in_tbi[chr_index],
                samples_file = get_sample_ids.output_file,
                regions = region,
                output_filename = "${output_basename}.chr${chr}.vcf.gz",
                output_type = output_type,
                cpu = bcftools_chr_cpu,
                mem_gb = bcftools_chr_mem_gb
        }

    }

    output{
        Array[File] ref_vcfs_out = subset_ref_vcf.vcf_out
        Array[File] genetic_map_files_out = format_rfmix_genetic_map.genetic_map_out
        File sample_file_out = subset_ancestry_samples.sample_file_out
        File rfmix_sample_file_out = make_rfmix_sample_file.output_file
    }
}

