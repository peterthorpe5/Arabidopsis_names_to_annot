# Define file paths:
vir1_up = 'up_vir1.names'
acd6_up = 'up_acd6.names'
vir1_down = 'down_vir1.names'
acd6_down = 'down_acd6.names'

# Initialize sets for each group of genes
all_up_genes = set()
all_down_genes = set()

vir_up_set = set()
acd6_up_set = set()

# Function to parse genes from a file and add them to a given set
def parse_gene_file(file_path, gene_set):
    with open(file_path, 'r') as file:
        for line in file:
            gene = line.strip()
            gene_set.add(gene)

# Parse genes from the 'up' files
parse_gene_file(vir1_up, all_up_genes)
parse_gene_file(vir1_up, vir_up_set)

parse_gene_file(acd6_up, all_up_genes)
parse_gene_file(acd6_up, all_up_genes)

# Parse genes from the 'down' files
parse_gene_file(vir1_down, all_down_genes)
parse_gene_file(acd6_down, acd6_up_set)

# Find common genes in the 'up' files
common_all_up_genes = all_up_genes.intersection(all_down_genes)

# Find genes different between 'up' and 'down' files
different_all_up_genes = all_up_genes.difference(all_down_genes)
different_all_down_genes = all_down_genes.difference(all_up_genes)

# Display the results
#print("Common genes in 'up' files:", common_all_up_genes)
#print("Genes different in 'up' files:", different_all_up_genes)
#print("Genes different in 'down' files:", different_all_down_genes)

different_up_in_acd6_not_up_vir1 = acd6_up_set.difference(vir_up_set)

#print("different_up_in_acd6_not_up_vir1:") # , different_up_in_acd6_not_up_vir1)
#for i in different_up_in_acd6_not_up_vir1:
#    print(i)



different_up_in_vir1_not_up_acd6 = vir_up_set.difference(acd6_up_set)

print("different_up_in_vir1_not_up_acd6:") # , different_up_in_acd6_not_up_vir1)
for i in different_up_in_vir1_not_up_acd6:
    print(i)

