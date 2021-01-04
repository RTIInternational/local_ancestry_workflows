import "local_ancestry_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as BCFTOOLS
import "local_ancestry_workflows/biocloud_wdl_tools/rfmix/rfmix.wdl" as RFMIX
import "local_ancestry_workflows/workflows/rfmix_chr_wf.wdl" as RFMIX_CHR_WF

workflow rfmix_wf{
    Array[File] query_vcfs
    Array[File] ref_vcfs
    Array[File] genetic_maps
    Array[String] chrs
    File sample_map
    String output_basename

    # Splitting parameters
    Boolean split_query_samples = true
    Int max_samples_per_split = 500

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
    Int rfmix_chr_cpu = 12
    Int rfmix_chr_mem_gb = 24
    Int bcftools_chr_cpu = 4
    Int bcftools_chr_mem_gb = 8

    # Loop through chrs and do RFMix on each
    scatter(chr_index in range(length(chrs))){
        String chr = chrs[chr_index]
        call RFMIX_CHR_WF.rfmix_chr_wf{
            input:
                query_vcf = query_vcfs[chr_index],
                ref_vcf = ref_vcfs[chr_index],
                sample_map = sample_map,
                genetic_map = genetic_maps[chr_index],
                output_basename = "${output_basename}.chr${chr}",
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
                rfmix_cpu = rfmix_chr_cpu,
                rfmix_mem_gb = rfmix_chr_mem_gb,
                bcftools_cpu = bcftools_chr_cpu,
                bcftools_mem_gb = bcftools_chr_mem_gb
        }
    }

    output{
        Array[File] msp_out = rfmix_chr_wf.msp_out
        Array[File] fb_out = rfmix_chr_wf.fb_out
        Array[File] q_out = rfmix_chr_wf.q_out
    }
}