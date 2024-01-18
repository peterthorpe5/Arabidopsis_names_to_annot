#!/bin/bash

# Input files
transcript_file="genes_feature_lengths.txt"
gene_to_go_file="Go_final.txt"

# Output file
output_file="transcript_to_go.txt"

# Create an associative array to store gene to GO term mappings
declare -A gene_to_go_mapping

# Read the gene to GO term file and populate the array
while IFS=$'\t' read -r gene go_terms; do
    gene_name="${gene%%.*}"  # Remove the .number to get the gene name
    gene_to_go_mapping["$gene_name"]=$go_terms
done < "$gene_to_go_file"

# Iterate through the transcript file and write transcript to GO term mappings
while IFS=$'\t' read -r transcript length; do
    gene_name="${transcript%%.*}"  # Remove the .number to get the gene name
    go_terms="${gene_to_go_mapping[$gene_name]}"
    
    if [ -n "$go_terms" ]; then
        echo -e "${transcript}\t${go_terms}"
    else
        echo -e "${transcript}\t-"
    fi
done < "$transcript_file" > "$output_file"

# Count the lines in the output file and the transcript file
line_count_output=$(wc -l < "$output_file")
line_count_transcript=$(wc -l < "$transcript_file")

# Check if the line counts match
if [ "$line_count_output" -ne "$line_count_transcript" ]; then
    echo "Error: Line counts do not match between output and transcript files."
else
    echo "Transcript to GO term mapping has been written to $output_file"
fi
