#!/usr/bin/env python3

# The purpose of this script is to take rna 
# transcripts where a secis element was identified, as input 
# and recode any instance of a TGA before the secis element to a TGC, 
# producing a new transcript for each recoding instance

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Load FASTA sequences into a dictionary
def load_fasta(fasta_file):
    return {record.id: record.seq for record in SeqIO.parse(fasta_file, 'fasta')}

# Read secis GFF file 
def read_gff(gff):
    # Read the GFF file and skip the header line
    df = pd.read_csv(gff, sep='\t', comment='#', header=None)
    # Set column names after reading
    df.columns = ['transcript_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    return df

# Need to filter to only transcripts with secis elements 
def filter_sequences(sequences, gff_df):
    target_ids = set(gff_df['transcript_id'])
    filtered = {}

    for seq_id, record in sequences.items():
        # Extract the base transcript id from the FASTA header
        if seq_id in target_ids:
            filtered[seq_id] = record

    return filtered   

# Identify positions of TGA codons in each transcript sequence
def find_tga_occurrences(sequences):
    tga_positions = {}
    
    for transcript_name, seq in sequences.items():
        seq_str = str(seq)  # Convert Seq object to string
        positions = [pos for pos in range(len(seq_str) - 2) if seq_str[pos:pos+3] == "TGA"]
        
        if positions:
            tga_positions[transcript_name] = positions  # Map transcript name to TGA positions
        else:
            print("No TGAs found for ", transcript_name)
    
    return tga_positions

# Filter TGA positions that occur after the secis element
def filter_TGAs(tga_positions, gff):
    for _, row in gff.iterrows():
        transcript = row['transcript_id']
        secis_start = row['start']
        if transcript in tga_positions: 
            pos_list = tga_positions[transcript]
            filtered_pos = [i for i in pos_list if i < secis_start]
            tga_positions[transcript] = filtered_pos
        else:
            continue
    
    return tga_positions 
          

# Perform recoding
def recode(sequences, recodon, tga_positions, output_fasta):
    records = []
    
    for transcript_name, seq in sequences.items():
        seq_str = str(seq)  # Convert Seq object to string
        positions = tga_positions.get(transcript_name, [])  # Get TGA positions or empty list

        # Store the original sequence first
        original_record = SeqRecord(Seq(seq_str), id=f"{transcript_name}_original", description="")
        records.append(original_record)
        
        for pos in positions:
            # Create a new sequence with TGA replaced by the recodon at position `pos`
            new_seq = seq_str[:pos] + recodon + seq_str[pos+3:]  
            # Create a new record for FASTA output
            record_id = f"{transcript_name}_TGA_{pos}_recode"
            record = SeqRecord(Seq(new_seq), id=record_id, description="")
            records.append(record)

    # Write all modified sequences to output FASTA
    with open(output_fasta, "w") as out_fasta:
        SeqIO.write(records, out_fasta, "fasta")

    print(f"Rewritten sequences saved to {output_fasta}")
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recode TGA codons at selenocysteine positions in a genome FASTA file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to genome FASTA file") 
    parser.add_argument('--gff', type=str, required=True, help="Path to secis gff file") 
    parser.add_argument('--recodon', type=str, required=True, choices=['TGT', 'TGC'], help="Codon to replace TGA (TGT or TGC)")
    parser.add_argument('--output', type=str, required=True, help="Output FASTA file of recoded sequences")
     
    args = parser.parse_args()
    
    fasta = load_fasta(args.fasta)  # Load sequences
    gff = read_gff(args.gff)
    
    fasta = filter_sequences(fasta, gff)
    TGAs = find_tga_occurrences(fasta)  # Find TGA positions
    filtered_TGAs = filter_TGAs(TGAs, gff)
    recode(fasta, args.recodon, filtered_TGAs, args.output)  # Recode sequences and save



    
# Command for running script 
'''
python recode_any_TGA.py --fasta gffread_out/transcripts_clean.fa --gff subset_mouse.all_secis.gff --recodon TGC --output recoded_mouse_test.fa
'''  