#!/bin/bash

# Specify the directory where your result files are located
result_dir="../results_vir/filtered_by_logFC_2.00_FDR_0.001/"

# Iterate over each gene in the "GOI" file
while IFS= read -r gene; do
    # Search for the gene in each result file and print the file names
    #grep -H "$gene" "$result_dir"* | awk '{print $1, $2, $3, $4, $5}' | grep -v "alkbh10c"
    grep -H "$gene" "$result_dir"* | grep -v "alkbh10c"
done < GOI
