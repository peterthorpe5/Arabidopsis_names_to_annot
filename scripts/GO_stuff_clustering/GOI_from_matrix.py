import argparse

def parse_gene_counts_matrix(matrix_file, goi_file, output_file):
    # Read the gene of interest (GOI) names from the file
    with open(goi_file, 'r') as goi_file:
        goi_names = set(line.strip() for line in goi_file)

    # Parse the gene counts matrix and write lines for genes of interest
    with open(matrix_file, 'r') as matrix_file, open(output_file, 'w') as output_file:
        # Read the header line
        header = matrix_file.readline().strip().split('\t')
        output_file.write('\t'.join(header) + '\n')

        for line in matrix_file:
            parts = line.strip().split('\t')
            gene_name = parts[0]

            if gene_name in goi_names:
                output_file.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse gene counts matrix for genes of interest.")
    parser.add_argument("-m", "--matrix", required=True, help="Path to the gene counts matrix file")
    parser.add_argument("--goi", required=True, help="Path to the gene of interest file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file")

    args = parser.parse_args()

    parse_gene_counts_matrix(args.matrix, args.goi, args.output)

