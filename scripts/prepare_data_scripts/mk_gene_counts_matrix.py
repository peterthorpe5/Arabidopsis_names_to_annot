import os

# Define the input files
gene_list_file = "gene.model.names"
matrix_file = "./transcript_DE/At.isoform.counts.matrix"

# Define the output file
output_file = "./genes_DE/At.isoform.counts.matrix"

missing_file = "missing.txt"

# Create a set to store gene names from the gene list
gene_names_set = set()

# Read gene names from the gene list file and store them in the set
gene_list_file_path = os.path.join(os.getcwd(), gene_list_file)

with open(gene_list_file_path, 'r') as gene_list:
    for line in gene_list:
        gene_names_set.add(line.strip())  # Include version number

# Initialize counts
expected_count = len(gene_names_set)
found_count = 0

# Initialize a list to store missing gene names
missing_genes = []

# Initialize a set to store novel gene names
novel_genes_set = set()

# Open the matrix file and the output file
matrix_file_path = os.path.join(os.getcwd(), matrix_file)
output_file_path = os.path.join(os.getcwd(), output_file)

with open(matrix_file_path, 'r') as matrix, open(output_file_path, 'w') as output:
    # Read the header (first line) from the matrix and write it to the output
    header = next(matrix)
    output.write(header)

    # Iterate through the matrix and write matching rows to the output file
    for line in matrix:
        parts = line.split('\t')
        gene_name = parts[0]
        if gene_name in gene_names_set:
            parts = line.split()
            gene_id = parts[0].split(".")[0]
            line = line.replace(gene_name, gene_id)
            output.write(line)
            found_count += 1
        elif gene_name.endswith(".novel1"):
            novel_genes_set.add(gene_name)

# Check missing genes against novel genes
for gene in gene_names_set:
    if gene not in novel_genes_set and gene not in missing_genes:
        missing_genes.append(gene)

# Write missing gene names to the missing file
missing_file_path = os.path.join(os.getcwd(), missing_file)

with open(missing_file_path, 'w') as missing_output:
    for gene in missing_genes:
        missing_output.write(gene + '\n')

# Print counts
print(f"Expected gene count: {expected_count}")
print(f"Found gene count: {found_count}")
print(f"Matching rows have been written to {output_file_path}")
print(f"Missing genes have been written to {missing_file_path}")
