
input_matrix_file = "At.gene.TMM.EXPR.matrix"  # Replace with your input matrix file
output_filtered_matrix_file = "At.gene.counts.min10.matrix"

min_expression_threshold = 6

with open(input_matrix_file, 'r') as infile, open(output_filtered_matrix_file, 'w') as outfile:
    header = infile.readline()
    outfile.write(header)

    for line in infile:
        parts = line.strip().split('\t')
        gene_name = parts[0]
        expressions = list(map(float, parts[1:]))  # Convert expression values to floats

        if any(expr > min_expression_threshold for expr in expressions):
            outfile.write(line)

print("Filtered matrix written to", output_filtered_matrix_file)


