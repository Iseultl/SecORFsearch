#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def create_diamond_db(fasta_file, db_prefix):
    """Create a DIAMOND database from a FASTA file."""
    cmd = f'diamond makedb --in {fasta_file} -d {db_prefix}'
    print(f'Creating DIAMOND database: {cmd}')
    subprocess.run(cmd, shell=True, check=True)
    return f'{db_prefix}.dmnd'
"""
def split_sequences_at_u(fasta_file, output_file):
    Split sequences at U position and write to new FASTA file.
    
    Args:
        fasta_file: Path to input FASTA file containing sequences with U
        output_file: Path to output FASTA file with sequences split at U
    
    Returns:
                before_u = seq[:u_pos]
                after_u = seq[u_pos+1:]  # Skip U itself
                # Write to output file
                with open(output_file, 'a') as out:
                    out.write(f'>{current_id}_before_U\n{before_u}\n')
                    out.write(f'>{current_id}_after_U\n{after_u}\n')
    
    print(f"Created split sequences file: {output_file}")
"""
def run_diamond_blastp(query_file, db_file, output_file, threads=16, min_identity=60.0):
    """Run DIAMOND blastp and analyze alignments of split sequences.
    
    Args:
        query_file: Path to FASTA file with sequences split at U positions
        db_file: Path to DIAMOND database file
        output_file: Path to output BLASTP results
        threads: Number of threads to use
        min_identity: Minimum percentage identity threshold
    
    Returns:
        DataFrame with alignment statistics
    """
    # Run DIAMOND blastp
    cmd = f'diamond blastp \
        -d {db_file} \
        -q {query_file} \
        -o {output_file} \
        --outfmt 0 \
        # --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        -p {threads} \
        --sensitive'
    print(f'Running DIAMOND blastp: {cmd}')
    subprocess.run(cmd, shell=True, check=True)

def process_alignments(blastp_file, min_identity=60.0):
    """Process DIAMOND blastp results."""
    alignment_stats = []
    
    with open(blastp_file) as f:
        for line in f:
            if line.startswith('#'):  # Skip header
                continue
                
            fields = line.strip().split('\t')
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            
            if pident < min_identity:
                continue
                
            # Store statistics
            alignment_stats.append({
                'query': qseqid,
                'subject': sseqid,
                'pident': pident,
                'length': length,
                'evalue': float(fields[10]),
                'bitscore': float(fields[11])
            })

    # Create DataFrame and write statistics
    stats_df = pd.DataFrame(alignment_stats)
    stats_output = blastp_file.replace('.blastp', '_stats.csv')
    stats_df.to_csv(stats_output, index=False)
    print(f'Alignment statistics written to: {stats_output}')
    
    return stats_df

def main():
    parser = argparse.ArgumentParser(
        description="Run all-vs-all DIAMOND blastp analysis of human sequences"
    )
    parser.add_argument(
        '--query_fasta', required=True,
        help="Path to query protein sequences FASTA file"
    )
    parser.add_argument(
        '--human_fasta', required=True,
        help="Path to human protein sequences FASTA file"
    )
    parser.add_argument(
        '--output_dir', required=True,
        help="Directory to store results (will be created if it doesn't exist)"
    )
    parser.add_argument(
        '--threads', type=int, default=16,
        help="Number of threads to use for DIAMOND blastp (default: 16)"
    )
    parser.add_argument(
        '--min_identity', type=float, required=True,
        help="Minimum percentage identity threshold for alignments (0-100)"
    )
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Create DIAMOND database
    human_db = create_diamond_db(args.human_fasta, os.path.join(args.output_dir, 'human_db'))
    query_db = create_diamond_db(args.query_fasta, os.path.join(args.output_dir, 'query_db'))

    # Run all-vs-all human vs human comparison
    human_vs_human_out = os.path.join(args.output_dir, 'human_vs_human.blastp')
    run_diamond_blastp(args.human_fasta, human_db, human_vs_human_out, args.threads, args.min_identity)
    #process_alignments(human_vs_human_out, args.min_identity)

    # Run all-vs-all query vs query comparison
    query_vs_query_out = os.path.join(args.output_dir, 'query_vs_query.blastp')
    run_diamond_blastp(args.query_fasta, query_db, query_vs_query_out, args.threads, args.min_identity)
    #process_alignments(query_vs_query_out, args.min_identity)

    # Run all-vs-all human vs query comparison
    human_vs_query_out = os.path.join(args.output_dir, 'human_vs_query.blastp')
    run_diamond_blastp(args.human_fasta, query_db, human_vs_query_out, args.threads, args.min_identity)
    #process_alignments(human_vs_query_out, args.min_identity)

    print("\nDIAMOND blastp all-vs-all analysis completed!")
    print(f"Results written to:")
    print(f"- Human vs Human: {human_vs_human_out}")
    print(f"- Alignment statistics: {human_vs_human_out.replace('.blastp', '_stats.csv')}")

if __name__ == "__main__":
    main()

# Command for running the script
"""
python run_diamond_all_vs_all.py \
    --query_fasta /Users/iseult/Desktop/query_Analysis/geneid_results.fa \
    --human_fasta /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results.fa \
    --output_dir /Users/iseult/Desktop/Geneid_Recoding/geneid_prot_alignments \
    --threads 16 \
    --min_identity 60
    
    # --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        
python run_diamond_all_vs_all.py \
    --query_fasta /Users/iseult/Desktop/Horse_Analysis/inframe_geneid_predictions.fa \
    --human_fasta /Users/iseult/Desktop/Human_Analysis/secis_update_160725/geneid_results.fa \
    --output_dir /Users/iseult/Desktop/Horse_Analysis/ \
    --threads 16 \
    --min_identity 60
"""
