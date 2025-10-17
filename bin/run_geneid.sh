#!/bin/bash

# Input parameters from Nextflow
GFF_FILE="$1"        # secis gff file
FASTA_FILE="$2"      # multi-fasta file with transcripts
OUTPUT_DIR="$3"      # geneid_recoded_predictions
PARAM_FILE="$4"      # geneid parameter file

# Create output dir
mkdir -p "$OUTPUT_DIR"

# Build SECIS coordinates list per transcript
declare -A secis_coords
while read -r line; do
    tid=$(echo "$line" | cut -f1)
    start=$(echo "$line" | cut -f4)
    # append coord to list (comma-separated)
    if [[ -n "${secis_coords[$tid]}" ]]; then
        secis_coords["$tid"]="${secis_coords[$tid]},$start"
    else
        secis_coords["$tid"]=$start
    fi
done < <(grep -v '^#' "$GFF_FILE")

echo "------ Loaded SECIS Elements ------"
for key in "${!secis_coords[@]}"; do
    echo "$key => ${secis_coords[$key]}"
done
echo "-----------------------------------"

# Iterate only over transcript IDs with SECIS hits
grep -v '^#' "$GFF_FILE" | cut -f1 | sort -u | while read -r tid; do
    coords="${secis_coords[$tid]}"
    echo "Processing transcript $tid with SECIS coords: $coords"

    tmp_fa=$(mktemp)
    seqkit grep -r -p "^$tid" "$FASTA_FILE" > "$tmp_fa"

    if [[ -s "$tmp_fa" ]]; then
        # Loop over each coord for this transcript
        IFS=',' read -ra coord_list <<< "$coords"
        for coord in "${coord_list[@]}"; do
            output_file="$OUTPUT_DIR/${tid}_SECIS${coord}.geneid.txt"
            echo "  → Running geneid with SECIS at $coord"
            $HOME/bin/geneid-stops -P "$PARAM_FILE" -s -W -k "$coord" "$tmp_fa" > "$output_file"
            echo "    Saved: $output_file"
        done
    else
        echo "Sequence for $tid not found in FASTA — skipping."
    fi
    rm -f "$tmp_fa"
done

