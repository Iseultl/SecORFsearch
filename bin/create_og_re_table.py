#!/usr/bin/env python3

# As input this script takes the best geneid scores
			#temp_test_longest.csv
			#temp_test_score.csv
			#temp_test_og_longest.csv
            #temp_test_og_score.csv
            
import pandas as pd
import argparse
import re

def safe_read_csv(path, expected_cols):
    """Read CSV, returning empty DataFrame with expected columns if file is empty."""
    try:
        df = pd.read_csv(path)
        if df.empty:
            return pd.DataFrame(columns=expected_cols)
        return df
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=expected_cols)

def read_tables(og_score, og_length, re_score, re_length):
    og_score_df = og_score.rename(columns={
        'start': 'og_score_start',
        'end': 'og_score_end',
        'score': 'og_score_score',
        'length': 'og_score_length'
    })
    og_length_df = og_length.rename(columns={
        'start': 'og_length_start',
        'end': 'og_length_end',
        'score': 'og_length_score',
        'length': 'og_length_length'
    })
    re_score_df = re_score.rename(columns={
        'start': 're_score_start',
        'end': 're_score_end',
        'score': 're_score_score',
        'length': 're_score_length'
    })
    re_length_df = re_length.rename(columns={
        'start': 're_length_start',
        'end': 're_length_end',
        'score': 're_length_score',
        'length': 're_length_length'
    })

    # Set index to transcript_name for all DataFrames
    og_score_df.set_index('transcript_name', inplace=True)
    og_length_df.set_index('transcript_name', inplace=True)
    re_score_df.set_index('transcript_name', inplace=True)
    re_length_df.set_index('transcript_name', inplace=True)
    
    

    # Apply the safe splitting function to create the new columns
    re_score_df['TGA_site_score'] = re_score_df['seqnames'].str.split('_').str[-3]
    re_length_df['TGA_site_longest'] = re_length_df['seqnames'].str.split('_').str[-3]
    # Select only the renamed columns (and drop other duplicate columns if necessary)
    og_score_df = og_score_df[['og_score_start', 'og_score_end', 'og_score_score', 'og_score_length']]
    og_length_df = og_length_df[['og_length_start', 'og_length_end', 'og_length_score', 'og_length_length']]
    re_score_df = re_score_df[['re_score_start', 're_score_end', 're_score_score', 're_score_length', 'TGA_site_score']]
    re_length_df = re_length_df[['re_length_start', 're_length_end', 're_length_score', 're_length_length', 'TGA_site_longest']]

    # Merge all DataFrames on transcript_name index
    merged_df = og_score_df.join([og_length_df, re_score_df, re_length_df], how='outer')
    
    return merged_df

def read_sel(gff):
    df = pd.read_csv(gff, sep='\t', names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    sel_df = df[df['type'] == 'Selenocysteine'].copy()
    if sel_df.empty:
        # Return empty DataFrame with correct columns
        return pd.DataFrame(columns=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']) 
    sel_df['attributes'] = sel_df['attributes'].str.split('_', n=1, expand=True)[1]
    sel_df = sel_df.set_index('attributes')
    sel_df = sel_df['type']

    return sel_df

def read_gff(gff):
    df = pd.read_csv(
        gff,
        sep='\t',
        names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
    
    # Prefer CDS features
    feature_df = df[df['type'] == 'CDS']
    # If no CDS present, fall back to transcript features
    if feature_df.empty:
        feature_df = df[df['type'].isin(['exon', 'transcript'])]
    # If still empty, return an empty frame with the right structure
    if feature_df.empty:
        return pd.DataFrame(columns=['gene_name', 'transcript_name', 'gtf_start', 'gtf_end']).set_index('transcript_name')

    # Group by transcript, keep gene_name as a regular column
    grouped = feature_df.groupby('attributes').agg({
        'start': 'min',
        'end': 'max'
    }).reset_index()
    
    grouped.rename(columns={
        'start': 'gtf_start',
        'end': 'gtf_end'
    }, inplace=True)

    grouped[['gene_name', 'transcript_name']] = grouped['attributes'].str.split('_', n=1, expand=True)
    grouped = grouped.set_index('transcript_name')
    grouped = grouped.drop('attributes', axis=1)
    
    return grouped

def gene_look_up(gff):
    df = pd.read_csv(
        gff,
        sep='\t',
        names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    ) 
    # Split attributes once into gene and transcript
    df[['gene_name', 'transcript_name']] = df['attributes'].str.split('_', n=1, expand=True)
    # Create dictionary: transcript_name -> gene_name
    gene_dict = df.drop_duplicates(subset=['transcript_name']) \
                  .set_index('transcript_name')['gene_name'] \
                  .to_dict()
    return gene_dict
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the coding regions of each transcript and input to coding_potential script")
    # Add arguments for input and output directories
    parser.add_argument('--og_score', type=str, required=True, help="Path to temp_test_og_score.csv")
    parser.add_argument('--og_length', type=str, required=True, help="Path to temp_test_og_longest.csv")  
    parser.add_argument('--re_score', type=str, required=True, help="Path to temp_test_score.csv")
    parser.add_argument('--re_length', type=str, required=True, help="Path to temp_test_longest.csv")  
    parser.add_argument('--gff', type=str, required=True, help="Path to gff relocated_to_transcript.gff")
    parser.add_argument('--output', type=str, required=True, help="Path to output table")
    
    args = parser.parse_args()
    sel_df = read_sel(args.gff)
    gff_df = read_gff(args.gff)

    og_score = safe_read_csv(args.og_score, ['seq','start','end','strand','type','length','score','gene_name','transcript_name'])
    og_length = safe_read_csv(args.og_length, ['seq','start','end','strand','type','length','score','gene_name','transcript_name'])
    re_score = safe_read_csv(args.re_score, ['seq','start','end','strand','type','length','score','gene_name','transcript_name'])
    re_length = safe_read_csv(args.re_length, ['seq','start','end','strand','type','length','score','gene_name','transcript_name'])
    gene_dict = gene_look_up(args.gff)
    merged = read_tables(og_score, og_length, re_score, re_length) 
    sel_df_filtered = sel_df.loc[sel_df.index.intersection(merged.index)]
    gff_df_filtered = gff_df.loc[gff_df.index.intersection(merged.index)]
    merged = merged.join(sel_df_filtered, how='left')
    merged = merged.join(gff_df_filtered, how='left')
    merged['gene_name'] = merged.index.map(gene_dict)
    merged.to_csv(args.output + '.csv', index=True)
    
    
#Command to execute code
"""
python create_og_re_table.py --og_score temp_test_og_score.csv --og_length temp_test_og_longest.csv --re_score temp_test_score.csv --re_length temp_test_longest.csv --gff /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/relocated_to_transcript.gff --secis /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/all_secis_positive.gff
python create_og_re_table.py --og_score temp_test_og_score.csv --og_length temp_test_og_longest.csv --re_score temp_test_score.csv --re_length temp_test_longest.csv --gff /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/relocated_start_stop.gff --secis /Users/iseult/Desktop/Geneid_Recoding/testing_false_positives/all_secis_positive.gff
"""
#!/usr/bin/env python3

# As input this script takes the best geneid scores
			#temp_test_longest.csv
			#temp_test_score.csv
			#temp_test_og_longest.csv
            #temp_test_og_score.csv
            

