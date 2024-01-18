import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Parse Excel and DE results files, apply filters, and save common DE results.")
    parser.add_argument("--excel-file", nargs='?', default="41586_2021_3315_MOESM4_ESM.xlsx", help="Path to the Excel file (default: 41586_2021_3315_MOESM4_ESM.xlsx)")
    parser.add_argument("--de-results-file", nargs='?', default="vir1_low_vs_col_vir_high.GLM.edgeR.DE_results", help="Path to DE results file (default: vir1_low_vs_col_vir_high.GLM.edgeR.DE_results)")
    parser.add_argument("--log2fc-threshold", type=float, default=2.0, help="Log2FC threshold (default: 2.0)")
    parser.add_argument("--pval-threshold", type=float, default=0.001, help="P-value threshold (default: 0.001)")
    args = parser.parse_args()

    try:
        # Read the Excel file
        data_excel = pd.read_excel(args.excel_file)
        
        # Read the DE results file
        data_de_results = pd.read_csv(args.de_results_file, sep='\t')
        
        # Apply filters to the Excel file
        positive_log2fc_excel = data_excel[(data_excel['log2FC'] > args.log2fc_threshold) & (data_excel['adj.pval'] < args.pval_threshold)]
        negative_log2fc_excel = data_excel[(data_excel['log2FC'] < -args.log2fc_threshold) & (data_excel['adj.pval'] < args.pval_threshold)]
        
        # Apply filters to the DE results file
        positive_log2fc_de = data_de_results[(data_de_results['logFC'] > args.log2fc_threshold) & (data_de_results['FDR'] < args.pval_threshold)]
        negative_log2fc_de = data_de_results[(data_de_results['logFC'] < -args.log2fc_threshold) & (data_de_results['FDR'] < args.pval_threshold)]
        
        # Collect gene names from both files using intersection
        positive_genes = set(positive_log2fc_excel['Target']).intersection(set(positive_log2fc_de['Row.names']))
        negative_genes = set(negative_log2fc_excel['Target']).intersection(set(negative_log2fc_de['Row.names']))
        
        # Collect all DE genes from both datasets
        all_de_genes = set(positive_genes.union(negative_genes)).intersection(set(positive_log2fc_excel['Target']).union(negative_log2fc_excel['Target']))


        # Collect the names of missing genes
        missing_genes = all_de_genes - (positive_genes.union(negative_genes))
        print(f"Number of missing genes: {len(missing_genes)}")
        # print("Missing genes:", missing_genes)

        # Create a folder for the common results
        common_results_folder = "common_results"
        os.makedirs(common_results_folder, exist_ok=True)


        # Print the number of genes in each intersection set
        print(f"Number of genes in positive intersection: {len(positive_genes)}")
        print(f"Number of genes in negative intersection: {len(negative_genes)}")
        print(f"Number of genes in ALL intersection: {len(positive_genes.union(negative_genes))}")

        print("Total genes in Excel file:", len(data_excel))
        print("Total genes in DE results file:", len(data_de_results))
        print("Number of genes after log2FC and p-value filter in Excel:", len(positive_log2fc_excel))
        print("Number of genes after log2FC and p-value filter in DE results:", len(positive_log2fc_de))
        
        # Write common DE genes to output files
        write_common_results("positive", common_results_folder, positive_genes, data_de_results, args)
        write_common_results("negative", common_results_folder, negative_genes, data_de_results, args)
        write_common_results("ALL", common_results_folder, all_de_genes, data_de_results, args)

    except FileNotFoundError:
        print("Error: File not found")
    except Exception as e:
        print(f"An error occurred: {e}")


def write_common_results(category, folder, common_genes, de_results_data, args):
    # Create a file name based on the category
    file_name = f"{category}_DE_FDR_{args.pval_threshold:.3f}_LFC{args.log2fc_threshold:.1f}.DE_results"
    file_path = os.path.join(folder, file_name)
    
    # Extract column titles from the DE results file
    column_titles = de_results_data.columns
    
    # Write the DE genes and column titles to the output file
    with open(file_path, 'w') as output_file:
        output_file.write('\t'.join(column_titles) + '\n')
        for gene in common_genes:
            try:
                result_row = de_results_data[de_results_data['Row.names'] == gene]
                # Convert all values to strings before writing
                output_file.write('\t'.join(map(str, result_row.values.tolist()[0])) + '\n')
            except IndexError:
                #print(f"Warning: Gene '{gene}' not found in DE results file.")
                # Silently continue if the gene is not found
                pass

if __name__ == "__main__":
    main()
