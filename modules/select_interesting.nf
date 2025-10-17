process SELECT_INTERESTING {
    tag { "select_interesting" }

    publishDir "${params.output_dir}/interesting_predictions", mode: 'copy'
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'

    input:
    path geneid_file
    path relocated_gtf

    output:
    path "recoded_score.csv"
    path "recoded_longest.csv"

    script:
    """
    select_interesting_recodings.py \
        --geneid ${geneid_file} \
        --gff ${relocated_gtf} \
        --score recoded_score.csv \
        --longest recoded_longest.csv
    """
}
