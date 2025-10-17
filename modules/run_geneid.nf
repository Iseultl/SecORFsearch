process RUN_GENEID {
    tag { "run_geneid_${input_file}" }
    
    publishDir "${params.output_dir}/geneid_recoded_predictions", mode: 'copy'
    cpus 1
    memory '8GB'
    label 'splitfasta'
    
    input:
    path(secis_gff)
    path(input_file)
    path(param_file)
    
    output:
    path("geneid_recoded_predictions/*.txt"), optional: true, emit: geneid_txt
    
    script:
    """
    bash run_geneid.sh \
        ${secis_gff} \
        ${input_file} \
        geneid_recoded_predictions \
        ${params.geneid_param}

    """
}
