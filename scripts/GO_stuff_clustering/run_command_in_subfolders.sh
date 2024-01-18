#!/bin/bash

# Define the command to be executed
command_to_run="python /Users/PThorpe001/Library/CloudStorage/OneDrive-UniversityofDundee/ggs_arabidopsis/expression_clustering_and_GO/vir1_temp/filter_GO.py"







# Define the FDR threshold (optional, default is 0.05)
fdr_threshold=0.05

# Use find to locate all subdirectories and execute the command
find "$input_dir" -type d -exec sh -c 'cd "{}" && '"$command_to_run"' {} --fdr '"$fdr_threshold" \;



# Use find to locate all subdirectories and execute the command
#find /Users/PThorpe001/Library/CloudStorage/OneDrive-UniversityofDundee/ggs_arabidopsis/expression_clustering_and_GO/vir1_temp/ -type d -exec sh -c 'cd "{}" && '"$command_to_run" \;
