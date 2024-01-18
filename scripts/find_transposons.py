import argparse

def check_transposons(tsv_file, transposon_file="transposon.names", output_file=None):
    # Read the list of transposon names
    with open(transposon_file, 'r') as f:
        transposon_names = [line.strip() for line in f]

    # Read the tab-separated file and check for transposon names
    with open(tsv_file, 'r') as f:
        lines = f.readlines()

    found_transposons = []

    for line in lines:
        # Split the line into columns using tab as the delimiter
        columns = line.strip().split('\t')

        # Check both columns for transposon names
        for col in columns[:2]:
            gene_name = col.rstrip('.1')
            if any(transposon_name in gene_name for transposon_name in transposon_names):
                found_transposons.append(gene_name)
                print(f"Transposon found in gene: {gene_name}")
                # Optionally, you can print or process the entire line or other information here

    # Write found transposons to the output file
    if output_file:
        with open(output_file, 'w') as out_file:
            for transposon in found_transposons:
                out_file.write(transposon + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for transposon names in a tab-separated file.")
    parser.add_argument("-t", "--tsv_file", required=True, help="Path to the tab-separated file")
    parser.add_argument("--transposon_file", default="transposon.names", help="Path to the transposon names file (default: transposon.names)")
    parser.add_argument("-o", "--output_file", help="Path to the output file for found transposons")

    args = parser.parse_args()

    check_transposons(args.tsv_file, args.transposon_file, args.output_file)
