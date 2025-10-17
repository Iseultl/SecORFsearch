#!/usr/bin/env python3

# The purpose of this script is to select the 
# best predictions on the original transcript sequences
#!/usr/bin/env python3

# The purpose of this script is to select the 
# best predictions on the original transcript sequences

import pandas as pd
import argparse
import numpy as np

def process_predictions(geneid_txt):
    best_by_score = {}
    best_by_length = {}
    current_seq = None

    with open(geneid_txt) as f:
        for line in f:
            if line.startswith("# Sequence"):
                # handle previous seq with no predictions
                if current_seq:
                    transcript_name = current_seq
                    gene_name = transcript_name
                    if transcript_name not in best_by_score:
                        placeholder = {
                            'seq': current_seq,
                            'start': np.nan, 'end': np.nan, 'strand': np.nan,
                            'type': np.nan, 'length': np.nan, 'score': np.nan,
                            'gene_name': gene_name,
                            'transcript_name': transcript_name
                        }
                        best_by_score[transcript_name] = placeholder
                        best_by_length[transcript_name] = placeholder

                current_seq = line.split()[2]  # transcript ID

            elif current_seq and line.startswith("  Single"):
                fields = line.strip().split()
                transcript_name = current_seq
                gene_name = transcript_name

                prediction = {
                    'seq': current_seq,
                    'start': int(fields[1]),
                    'end': int(fields[2]) - 3,
                    'strand': '+',
                    'type': fields[0],
                    'length': abs((int(fields[2]) - 3) - int(fields[1])) + 1,
                    'score': float(fields[3]),
                    'gene_name': gene_name,
                    'transcript_name': transcript_name
                }

                # Update best by score
                if (
                    transcript_name not in best_by_score
                    or pd.isna(best_by_score.get(transcript_name, {}).get('score'))
                    or prediction['score'] > best_by_score[transcript_name]['score']
                ):
                    best_by_score[transcript_name] = prediction

                # Update best by length
                current_best_len = best_by_length.get(transcript_name, {}).get('length')
                if (
                    current_best_len is None
                    or pd.isna(current_best_len)
                    or prediction['length'] > current_best_len
                ):
                    best_by_length[transcript_name] = prediction

        # Handle the very last sequence
        if current_seq:
            transcript_name = current_seq
            gene_name = transcript_name
            if transcript_name not in best_by_score:
                placeholder = {
                    'seq': current_seq,
                    'start': np.nan, 'end': np.nan, 'strand': np.nan,
                    'type': np.nan, 'length': np.nan, 'score': np.nan,
                    'gene_name': gene_name,
                    'transcript_name': transcript_name
                }
                best_by_score[transcript_name] = placeholder
                best_by_length[transcript_name] = placeholder

    score_df = pd.DataFrame(list(best_by_score.values()))
    longest_df = pd.DataFrame(list(best_by_length.values()))

    return score_df, longest_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select best predictions on original transcript sequences"
    )
    parser.add_argument('--score', type=str, required=True, help="Path to output file - csv of predictions (by score)")
    parser.add_argument('--longest', type=str, required=True, help="Path to output file - csv of predictions (by length)")
    parser.add_argument('--geneid', type=str, required=True, help="Path to geneid txt file")

    args = parser.parse_args()

    score_df, longest_df = process_predictions(args.geneid)

    score_df.to_csv(args.score, index=False)
    longest_df.to_csv(args.longest, index=False)

    
# Code for running script
"""
python get_og_predictions.py --geneid /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/geneid_results_originals.txt --output temp_test_og
"""