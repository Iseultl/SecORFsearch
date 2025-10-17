process RECODE_TGA {
    tag { "recode_tga" }
    publishDir "${params.output_dir}/recoded_transcripts", mode: 'copy'
    label 'python'
    cpus 1
    memory '8GB'
    time '1h'
    
    input:
    path(transcript_fasta)
    path(filtered_secis_gff)
    
    output:
    path("*.fa")
    
    script:
    def base = transcript_fasta.getBaseName()
    def out_name = "${base}_recoded.fa"
    """
    recode_any_TGA.py \
        --fasta ${transcript_fasta} \
        --gff ${filtered_secis_gff} \
        --recodon TGC \
        --output ${out_name}
    """
}
