process GET_ORIGINAL_PREDICTIONS {
    tag { "get_original" }
    label 'python'

    publishDir "${params.output_dir}/original_predictions", mode: 'copy'

    cpus 1
    memory '8GB'
    time '1h'

    input:
    path geneid_file

    output:
    path "original_score.csv"
    path "original_longest.csv"

    script:
    """
    get_og_predictions.py \
        --geneid ${geneid_file} \
        --score original_score.csv \
        --longest original_longest.csv
    """
}
