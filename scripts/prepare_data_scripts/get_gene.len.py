import os

# Define the input files
gene_list_file = "gene.model.names"
length_file = "genes_feature_lengths.txt"  # Corrected file name

# Define the output file
output_file = "genes_with_lengths.txt"

# Create a set to store gene names from the gene list
gene_names_set = set()

# Read gene names from the gene list file and store them in the set
gene_list_file_path = os.path.join(os.getcwd(), gene_list_file)

with open(gene_list_file_path, 'r') as gene_list:
    for line in gene_list:
        gene_names_set.add(line.strip())  # Include version number

# Initialize a dictionary to store gene names and their corresponding lengths
gene_lengths = {}

# Read gene names and lengths from the length file
length_file_path = os.path.join(os.getcwd(), length_file)

with open(length_file_path, 'r') as length_file:
    for line in length_file:
        parts = line.strip().split()
        gene_name = parts[0]
        length = int(parts[1])
        gene_lengths[gene_name] = length

# Initialize a list to store missing gene names
missing_genes = []

# Initialize a set to store novel gene names
novel_genes_set = set()

# Open the output file
output_file_path = os.path.join(os.getcwd(), output_file)

with open(output_file_path, 'w') as output:
    # Write header to the output
    output.write("Gene\tLength\n")

    # Iterate through the gene names and lengths and write matching genes to the output file
    for gene_name, length in gene_lengths.items():
        if gene_name in gene_names_set:
            # Remove version number from gene name
            gene_id = gene_name.split(".")[0]
            output.write(f"{gene_id}\t{length}\n")
        elif gene_name.endswith(".novel1"):
            novel_genes_set.add(gene_name)

# Check missing genes against novel genes
for gene in gene_names_set:
    if gene not in novel_genes_set and gene not in missing_genes:
        missing_genes.append(gene)

# Print counts
print(f"Found {len(gene_lengths)} genes with lengths.")
print(f"Matching genes and their lengths have been written to {output_file_path}")
print(f"Missing genes: {', '.join(missing_genes)}")
