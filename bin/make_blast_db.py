#!/usr/bin/env python3

# The purpose of this script is to select the in-frame recodings 
# extract their n sequences and make a blast db
# The recoded cysteines need to be converted to U
# The output is a fasta file containing all the in-frame predictions

import pandas as pd
import re
import argparse
import numpy as np



def extract_tga_from_header(header):
    match = re.search(r'_TGA_(\d+)_', header)
    print(match)
    return int(match.group(1)) if match else None

def find_cysteine_position(protein_seq, start_pos, tga_pos):
    """Find the position of the cysteine (C) amino acid within ±1 window of the calculated position.
    Returns the position of the C if found, otherwise None."""
    if tga_pos is None:
        return None
    
    # Calculate base position
    base_pos = (tga_pos - start_pos) // 3
    
    # Search within ±1 window
    for offset in [-1, 0, 1]:
        pos = base_pos + offset
        if 0 <= pos < len(protein_seq) and protein_seq[pos] == 'C':
            return pos
    
    return None

def is_in_frame(tga_pos, start, end):
    if tga_pos is not None:
        return start <= tga_pos <= end and ((tga_pos + 1) - start) % 3 == 0
    else:
        return False

def extract_protein_sequences(geneid_txt):
    """Extract protein sequences from geneid output, including -1 frame sequences.
    Returns a dictionary where keys are sequence IDs with frame information and values are protein sequences.
    Each sequence will be stored separately for each frame."""
    protein_sequences = {}
    
    # Define valid amino acid characters (including U for selenocysteine)
    VALID_AA_CHARS = set('ACDEFGHIKLMNPQRSTVWYUX*')
    
    with open(geneid_txt) as f:
        current_seq_id = None
        current_protein = []
        current_tga_pos = None
        current_start = None
        current_frame = None
        
        for line in f:
            if line.startswith("# Sequence"):
                # Store current sequence if we have one
                if current_seq_id and current_protein:
                    # Join and validate the sequence
                    seq = ''.join(current_protein)
                    if all(c in VALID_AA_CHARS for c in seq):
                        protein_sequences[f"{current_seq_id}+{current_frame}"] = seq
                    else:
                        print(f"Warning: Invalid characters in sequence {current_seq_id}+{current_frame}")
                
                # Start new sequence
                current_seq_id = line.split()[2]
                current_protein = []
                current_tga_pos = extract_tga_from_header(current_seq_id)
                current_frame = 0  # Reset frame counter
            elif line.startswith("  Single"):
                fields = line.strip().split()
                start = int(fields[1])
                end = int(fields[2]) - 3
                
                if current_tga_pos is not None and is_in_frame(current_tga_pos, start, end):
                    aa_seq = fields[-1]
                    if '_TGA_' not in aa_seq and aa_seq:  # Only add valid protein sequences
                        # Validate the sequence
                        if all(c in VALID_AA_CHARS for c in aa_seq):
                            # Store current sequence if we have one
                            if current_protein:
                                protein_sequences[f"{current_seq_id}+{current_frame}"] = ''.join(current_protein)
                                current_frame += 1
                            current_protein = [aa_seq]  # Start new sequence
                            current_start = start
                        else:
                            print(f"Warning: Invalid characters in sequence {current_seq_id}+{current_frame}")
        
        # Add the last sequence if it exists
        if current_seq_id and current_protein:
            # Join and validate the sequence
            seq = ''.join(current_protein)
            if all(c in VALID_AA_CHARS for c in seq):
                protein_sequences[f"{current_seq_id}+{current_frame}"] = seq
            else:
                print(f"Warning: Invalid characters in sequence {current_seq_id}+{current_frame}")
    
    return protein_sequences

def write_fasta_from_df(df, output_file):
    """Write protein sequences from DataFrame to FASTA file.
    
    Args:
        df: DataFrame containing transcript_name and aa columns
        output_file: Path to output FASTA file
    """

    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            header = f'>{row["transcript_name"]}'
            sequence = row["aa"]
            f.write(f'{header}\n{sequence}\n')
            
    print(f'FASTA file written to: {output_file}')
    print(f'Total sequences written: {len(df)}')

def select_target_predictions(geneid_txt):
    data = []
    total_predictions = 0
    out_of_frame = 0
    no_cysteine = 0
    cysteine_count = 0
    invalid_char = 0
    added_empty = set()
    
    # Extract protein sequences
    protein_sequences = extract_protein_sequences(geneid_txt)
    
    with open(geneid_txt) as f:
        current_seq = None
        current_tga_pos = None
        found_prediction = False
        prediction_count = 0
        original_count = 0

        for line in f:
            if line.startswith("# Sequence"):
                # Add empty row for previous sequence if needed
                if current_seq and not found_prediction and "original" not in current_seq and current_seq not in added_empty:
                    data.append([
                        '1', np.nan, np.nan, '+', np.nan, np.nan,
                        current_seq, current_seq, 'protein_coding', np.nan
                    ])
                    added_empty.add(current_seq)

                # Start new sequence
                current_seq = line.split()[2]
                current_tga_pos = extract_tga_from_header(current_seq)
                found_prediction = False
                prediction_count = 0

            elif line.startswith("  Single"):
                total_predictions += 1
                fields = line.strip().split()
                start = int(fields[1])
                end = int(fields[2]) - 3

                if "original" in current_seq:
                    original_count += 1
                    continue  # skip "original" sequences

                if is_in_frame(current_tga_pos, start, end):
                    # Get the sequence ID from the current_seq (without the TGA info)
                    seq_id = current_seq
                    aa = fields[-1]
                    
                    # Validate amino acid sequence
                    VALID_AA_CHARS = set('ACDEFGHIKLMNPQRSTVWYUX*')
                    if not all(c in VALID_AA_CHARS for c in aa):
                        print(f"Warning: Invalid characters in sequence {current_seq}: {set(aa) - VALID_AA_CHARS}")
                        invalid_char += 1
                        continue
                    
                    # Find the cysteine position within ±1 window
                    cysteine_pos = find_cysteine_position(aa, start, current_tga_pos)
                    if cysteine_pos is not None:
                        cysteine_count += 1
                        # Replace the cysteine with a selenocysteine
                        aa = aa[:cysteine_pos] + 'U' + aa[cysteine_pos+1:]
                        print(f"Sequence: {seq_id}, Cysteine position: {cysteine_pos}, AA: {aa}, start: {start}, end: {end}, tga_pos: {current_tga_pos}")
                    else:
                        print(f"Warning: No cysteine found in {current_seq} at position {current_tga_pos}")
                        no_cysteine += 1
                        continue
                    
                    found_prediction = True
                    structure = fields[0]
                    score = float(fields[3])
                    strand = '+'
                    seqname = current_seq + '+' + str(prediction_count)
                    gene_name = current_seq

                    data.append([
                        '1', start, end, strand, structure, score,
                        gene_name, seqname, 'protein_coding', aa
                    ])
                    prediction_count += 1
                else:
                    out_of_frame += 1

        # Handle last sequence if it had no predictions
        if current_seq and not found_prediction and "original" not in current_seq and current_seq not in added_empty:
            data.append([
                '1', np.nan, np.nan, '+', np.nan, np.nan,
                current_seq, current_seq, 'protein_coding', np.nan
            ])
            added_empty.add(current_seq)

    df = pd.DataFrame(data, columns=[
        'seqnames', 'start', 'end', 'strand', 'type', 'score',
        'gene_name', 'transcript_name', 'transcript_biotype', 'aa'   
    ])
    

    print("Total Geneid Predictions:", total_predictions)
    print("Filtered Out of Frame:", out_of_frame)
    print("No cysteine in prediction:", no_cysteine)
    print("Cysteine in prediction:", cysteine_count)
    print("Total Transcripts in Output:", len(df))
    print("Totail invalid seq:", invalid_char)
    x = total_predictions - out_of_frame
    print("In-frame:", x)
    y = x - original_count
    print("Remove originals:", y)
    total_count = y - invalid_char
    print("Recoded predictions in frame:", total_count)
    
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select interesting transcript sequences")
    # Add arguments for input and output directories
    parser.add_argument('--geneid_scores', type=str, required=True, help="Path to geneid scores of all ORFs")
    parser.add_argument('--output_fasta', type=str, required=True, help="Path to output FASTA file")
    args = parser.parse_args()
    
    
    df = select_target_predictions(args.geneid_scores)
    df = df.dropna(subset=["aa"])
    
    # Write FASTA file
    write_fasta_from_df(df, args.output_fasta)
    
    
# Command for running the script 
"""
python make_blast_db.py --geneid_scores /Users/iseult/Desktop/Human_Analysis/secis_update_160725/geneid_recoded_predictions.txt --output_fasta /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results.fa
python make_blast_db.py --geneid_scores /Users/iseult/Desktop/Mouse_Analysis/geneid_all_predictions.txt --output_fasta /Users/iseult/Desktop/Mouse_Analysis/geneid_results.fa
python make_blast_db.py --geneid_scores /Users/iseult/Desktop/Horse_Analysis/geneid_recoded_predictions.txt --output /Users/iseult/Desktop/Horse_Analysis/inframe_geneid_predictions.fa
"""