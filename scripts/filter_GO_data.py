import os
import argparse
import pandas as pd

def filter_data(input_file, output_dir=os.getcwd(), fdr_threshold=0.05):
    # Load the data from the input file
    data = pd.read_csv(input_file, sep='\t')

    # Determine if the input file is enriched or depleted
    is_enriched = "enriched" in input_file

    # Filter the data based on the FDR threshold
    filtered_data = data[data['over_represented_FDR'] < fdr_threshold]

    # Get the filename without the extension and add "enriched" or "depleted" based on the input
    file_name = os.path.splitext(os.path.basename(input_file))[0]
    file_name += "_enriched" if is_enriched else "_depleted"

    # Create a subdirectory for filtered files
    filtered_dir = os.path.join(output_dir, f'filtered_P{fdr_threshold}')
    if not os.path.exists(filtered_dir):
        os.makedirs(filtered_dir)

    # Save the filtered data to a new file in the filtered directory
    output_file = os.path.join(filtered_dir, f'{file_name}_P{fdr_threshold}.txt')
    filtered_data.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Filter data based on FDR threshold.')
    parser.add_argument('input_dir', help='Input directory containing data files')
    parser.add_argument('--fdr', type=float, default=0.05, help='FDR threshold (default: 0.05)')
    args = parser.parse_args()

    # List all files in the input directory
    input_files = [f for f in os.listdir(args.input_dir) if f.endswith(('GOseq.depleted', 'GOseq.enriched'))]

    for input_file in input_files:
        input_path = os.path.join(args.input_dir, input_file)
        filter_data(input_path, args.input_dir, args.fdr)

if __name__ == '__main__':
    main()
