/* 
 * NextFlow test secis pipeline
 * @authors
 * Iseult Leahy <iseult.leahy@crg.eu>
 * 
 */

/*
 * Input parameters: genome fasta file, genome gtf file
 * Params are stored in the params.config file
 */

version                 = "1.0"

// Singularity configuration is set in nextflow.config

// Genome and annotation files are expected to be provided by the user,
// for example, in a params.yaml file specified with the `-params-file` flag.
params.debug = params.debug ?: true
params.max_cpus = params.max_cpus ?: 4
params.max_memory = params.max_memory ?: '8GB'
params.dev_mode = params.dev_mode ?: false

// Input files
params.genome_gtf = params.genome_gtf ?: '/Users/iseult/gitlab/secis/test_data/gencode.vM37.annotation.gtf'
params.genome_fasta = params.genome_fasta ?: '/Users/iseult/gitlab/secis/test_data/GRCm39.primary_assembly.genome.fa'
params.output_dir = params.output_dir ?: './output'
params.geneid_param = params.geneid_param ?: '/Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/human3iso.param'

// Set smaller resource limits for development mode
if (params.dev_mode) {
    params.max_cpus = 2
    params.max_memory = '2GB'
}

// this prints the input parameters
log.info """
BIOCORE@CRG - N F TESTPIPE  ~  version ${version}
=============================================
genome_gtf                           : ${params.genome_gtf}
genome_fasta                         : ${params.genome_fasta}
geneid_param                         : ${params.geneid_param}
debug                                : ${params.debug}
max_cpus                             : ${params.max_cpus}
max_memory                           : ${params.max_memory}
"""

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the SECIS search and analysis pipeline'
    log.info 'Features:'
    log.info '  - Extract spliced transcript sequences'
    log.info '  - Search for SECIS elements'
    log.info '  - Recode TGA to TGC in transcripts with SECIS'
    log.info '  - Run geneid predictions on recoded and original transcripts'
    log.info '  - Analyze and compare predictions'
    log.info '\n'
    log.info 'Parameters:'
    log.info '  --dev_mode   : Enable development mode (uses less resources)'
    log.info '  --secis_gff  : Use pre-existing SECIS GFF file instead of running search'
    log.info '\n'
    exit 1
}

Channel
    .fromPath(params.genome_gtf, checkIfExists: true)
    .set { ch_genome_gtf }

Channel
    .fromPath(params.genome_fasta, checkIfExists: true)
    .set { ch_genome_fasta }

include { AGAT_GFF2GTF } from './modules/agat_gff2gtf'
include { AGAT_SPLITGFF } from './modules/agat_splitgff'
include { GFFREAD } from "${projectDir}/modules/gffread" 
include { SPLITFASTA } from "${projectDir}/modules/splitfasta"
include { SPLIT_IF_TOO_LARGE } from "${projectDir}/modules/split_if_too_large"
include { SECISSEARCH } from "${projectDir}/modules/secissearch"
include { FILTER_SECIS } from "${projectDir}/modules/filter_secis"
include { MERGE_GFF } from "${projectDir}/modules/merge_gff"
include { RECODE_TGA } from "${projectDir}/modules/recode_tga"
include { RUN_GENEID } from "${projectDir}/modules/run_geneid"
include { RUN_GENEID_ORIGINAL } from "${projectDir}/modules/run_geneid_original"
include { SELECT_INTERESTING } from "${projectDir}/modules/select_interesting"
include { GET_ORIGINAL_PREDICTIONS } from "${projectDir}/modules/get_original_predictions"
include { CLEAN_GTF } from "${projectDir}/modules/clean_gtf"
include { CONCATENATE_GTFS } from './modules/concatenated_gff'
include { RELOCATE_TRANSCRIPTS } from "${projectDir}/modules/relocate_transcripts"
include { CREATE_SUMMARY_TABLE } from "${projectDir}/modules/create_summary_table"
include { FILTER_FINAL_TABLE } from './modules/filter_final_table'
include { CONCAT_SUMMARY_RESULTS } from './modules/concat_summary_results'
include { CONCATENATE_RE } from './modules/concatenate_re'
include { CONCATENATE_OG } from './modules/concatenate_og'

// Helper function to extract chromosome name (without extension)
def get_chr_name(file) {
    return file.getBaseName().replaceFirst(/\.fa$|\.gtf$/, '')
}

workflow {
    // Step 1: Split GFF by chromosome
    split_results = AGAT_SPLITGFF(params.genome_gtf)

    // Step 2: Collect all .gff files from the output directory
    gff_files_ch = split_results.gff_files
    gff_files_ch = gff_files_ch.flatten()

    // Step 3: Run AGAT_GFF2GTF in parallel to standardise each GFF file
    gtf_files_ch = AGAT_GFF2GTF(gff_files_ch)
    
    // Create a channel of cleaned GTF files for downstream processing
    split_gff_dir_ch = CLEAN_GTF(gtf_files_ch)
    
    // Create paired channels for GTF and FASTA files
    // First, create a channel for GTF files with chromosome names
    base_names_gtf = gtf_files_ch.map { file ->
        def chr = file.name.replaceFirst(/\.gtf$/, '')
        tuple(chr, file)
    }
    
    // Step 4: FASTA Processing Pipeline
    split_fasta_dir_ch = SPLITFASTA(params.genome_fasta)
    base_names_fasta = split_fasta_dir_ch.flatten().map { file ->
        def fname = file.name  // e.g. horse_genome.part_NW_027222397.1.fa
        def chr_match = fname =~ /part_(.+)\.fa/
        def chr = chr_match ? chr_match[0][1] : null
        tuple(chr, file)
    }
    
    // Step 5: Create paired gtf & fasta channel
    paired_ch = base_names_gtf.combine(base_names_fasta, by: 0) 
    // Step 6: Create relocated to transcript gff 
    CONCATENATE_GTFS(split_gff_dir_ch.collect())
    relocated_gtf = RELOCATE_TRANSCRIPTS(CONCATENATE_GTFS.out)
    
    // Step 7: Now pass the paired channel to GFFREAD
    gffread_out = GFFREAD(paired_ch)

    // Step 8: Run SECIS search in parallel
    secis_results = SECISSEARCH(gffread_out)
    secis_results.collect().set { all_gffs }
    
    // Step 9: Merge results
    merged_secis_gff = MERGE_GFF(all_gffs)

    // Step 10: Filter SECIS elements
    filtered_secis_gff = FILTER_SECIS(merged_secis_gff)
        
    // Step 11: Recode TGA to TGC in transcripts with SECIS
    recoded_files = RECODE_TGA(gffread_out, filtered_secis_gff)

    // Step 11. Split recoded transcripts if too large
    split_recoded_files = SPLIT_IF_TOO_LARGE(recoded_files).split_fasta.flatten().filter { file -> file.size() > 0 }
    
    // Step 12: Run geneid predictions on recoded transcripts
    geneid_results_re = RUN_GENEID(filtered_secis_gff, split_recoded_files, params.geneid_param).geneid_txt 
    flat_geneid_results_re = geneid_results_re.flatten().unique()
    all_geneid_results_re = CONCATENATE_RE(flat_geneid_results_re.collect())

    // Step 13: Run geneid predictions on original transcripts, after SECIS GFF is ready
    geneid_results_og = RUN_GENEID_ORIGINAL(filtered_secis_gff, gffread_out, params.geneid_param)
    all_geneid_results_og = CONCATENATE_OG(geneid_results_og.collect())

    // Step 14: Collect and process predictions
    original_predictions = GET_ORIGINAL_PREDICTIONS(all_geneid_results_og)
    
    // Step 15: Select interesting predictions
    interesting_predictions = SELECT_INTERESTING(all_geneid_results_re, relocated_gtf)
    
    // Step 16: Final Output - Pass the combined channel to the summary table process
    SecORFsearch_result = CREATE_SUMMARY_TABLE(interesting_predictions, original_predictions, relocated_gtf)



}

workflow.onComplete { 
    println ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}