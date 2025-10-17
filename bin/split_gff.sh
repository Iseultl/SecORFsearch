#!/bin/bash

# Input parameters
INPUT_FILE="$1"
OUTPUT_DIR="tmp"

# Create output directory
mkdir -p "$OUTPUT_DIR"

gawk -F '\t' '{ if ($1 !~ /^#/) { print $0 > "$OUTPUT_DIR/"$1".gff" } }' "$INPUT_FILE"

