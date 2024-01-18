# Define a dictionary to store GO terms for each gene
gene_to_go = {}

# Read the input file and populate the dictionary
with open('gen_to_go', 'r') as file:
    for line in file:
        print(line)
        gene, go_term = line.strip().split()
        if gene in gene_to_go:
            gene_to_go[gene].append(go_term)
        else:
            gene_to_go[gene] = [go_term]

# Write the output file with comma-separated GO terms for each gene
with open('Go_final.txt', 'w') as file:
    for gene, go_terms in gene_to_go.items():
        file.write(f"{gene}\t{','.join(go_terms)}\n")
