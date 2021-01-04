import "local_ancestry_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as BCFTOOLS
import "local_ancestry_workflows/biocloud_wdl_tools/rfmix/rfmix.wdl" as RFMIX

workflow rfmix_chr_wf{
    File query_vcf
    File? query_vcf_tbi
    File ref_vcf
    File sample_map
    File genetic_map
    String output_basename
    String chr

    # Optional list of sample files for chunking by sample subroups
    Array[File]? sample_splits

    # Region to analyze <start_pos>-<end_pos>, in Mbp (decimal allowed)
    String? region

    # Random seed
    String? random_seed

    # Optional stuff
    Float? crf_spacing
    Float? rf_window_size
    Float? crf_weight
    Float? generations
    Int? em_iterations
    Boolean reanalyze_reference = false
    Int? node_size
    Int? trees
    Float? max_missing
    Int? bootstrap_mode
    Int? rf_min_snps
    String? debug_flag

    # Resource options
    Int rfmix_cpu = 12
    Int rfmix_mem_gb = 24
    Int bcftools_cpu = 4
    Int bcftools_mem_gb = 8

    # Split file by sample and region
    #if(defined(sample_splits)){
    #    Array[File] actual_sample_splits = select_first([sample_splits])
    #    scatter(split_index in range(len(actual_sample_splits))){
    #        call BCFTOOLS.view as make_sample_chunks{
    #            input:
    #                vcf_in = query_vcf,
    #                vcf_tbi_in = query_vcf_tbi,
    #                samples_file = actual_sample_splits[split_index],
    #                regions = region,
    #                output_filename = "${output_basename}.split${split_index}.vcf.gz",
    #                output_type = "z",
    #                cpu = bcftools_cpu,
    #                mem_gb = bcftools_mem_gb
    #        }
    #    }
    #}

    # Optionally split by region
    #if(defined(region) && !defined(sample_splits)){
    #    call BCFTOOLS.view as subset_region{
    #        input:
    #            vcf_in = query_vcf,
    #            vcf_tbi_in = query_vcf_tbi,
    #            regions = region,
    #            output_filename = "${output_basename}.region_subset.vcf.gz",
    #            output_type = "z",
    #            cpu = bcftools_cpu,
    #            mem_gb = bcftools_mem_gb
    #    }
    #}

    # Loop through splits and do RFMix on each
    #Array[File] split_query_vcfs = select_first([make_sample_chunks.vcf_out, [select_first([subset_region.vcf_out, query_vcf])]])
    #scatter(rfmix_split_index in range(length(split_query_vcfs))){
    call RFMIX.rfmix{
        input:
            query_vcf = query_vcf,
            ref_vcf = ref_vcf,
            sample_map = sample_map,
            genetic_map = genetic_map,
            output_basename = "${output_basename}.rfmix",
            chr = chr,
            random_seed = random_seed,
            crf_spacing = crf_spacing,
            rf_window_size = rf_window_size,
            crf_weight = crf_weight,
            generations = generations,
            em_iterations = em_iterations,
            reanalyze_reference = reanalyze_reference,
            node_size = node_size,
            trees = trees,
            max_missing = max_missing,
            bootstrap_mode = bootstrap_mode,
            rf_min_snps = rf_min_snps,
            debug_flag = debug_flag,
            cpu = rfmix_cpu,
            mem_gb = rfmix_mem_gb
    }

    # Do something to merge these files
    output{
        File msp_out = rfmix.msp_out
        File fb_out = rfmix.fb_out
        File q_out = rfmix.q_out
    }
}