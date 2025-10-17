process RUN_GENEID_ORIGINAL {
    tag { "run_geneid_original_${input_file}" }
    
    publishDir "${params.output_dir}/geneid_original_predictions", mode: 'copy'
    cpus 1
    memory '8GB'
    label 'splitfasta'
    
    input:
    path(secis_gff)
    path(input_dir)
    path(param_file)
    
    output:
    path("geneid_original_predictions/*.txt"), optional: true
    
    script:
    """
    # Run geneid on the transcripts
    bash run_geneid_original.sh \
        ${secis_gff} \
        ${input_dir} \
        geneid_original_predictions \
        ${params.geneid_param}

    """
}
