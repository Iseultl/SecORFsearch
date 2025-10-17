process CREATE_SUMMARY_TABLE {
    tag { "create_summary_table" }
    label 'python'
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 1
    memory '4GB'
    time '2h'
    
    input:
    path "recoded_score.csv"
    path "recoded_longest.csv"
    path "original_score.csv"
    path "original_longest.csv"
    path relocated_gtf
    
    output:
    path("SecORFsearch_result.csv")
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    # Create output folder
    mkdir -p summary_results

    # Run script with output prefix
    create_og_re_table.py \
        --og_score original_score.csv \
        --og_length original_longest.csv \
        --re_score recoded_score.csv \
        --re_length recoded_longest.csv \
        --gff ${relocated_gtf} \
        --output SecORFsearch_result
    """
}
