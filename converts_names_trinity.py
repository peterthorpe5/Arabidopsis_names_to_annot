import os
import gzip
from collections import defaultdict
import pandas as pd

def test_line(line):
    """
    Determine whether a line is valid (not a comment or blank).

    Args:
        line (str): A line of text from the file.

    Returns:
        bool: True if the line is valid, False otherwise.
    """
    if not line.strip():
        return False
    if line.startswith("#"):
        return False
    return True

def parse_flowering_gene():
    """
    Parse the list of flowering time genes.

    Returns:
        set: A set of gene names.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", "flowering_time_genes.gz")
    flowering_genes = set()
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                gene = line.strip()
                flowering_genes.add(gene.upper())
    return flowering_genes

def parse_go():
    """
    Parse the GO terms file and return a dictionary.

    Returns:
        dict: A dictionary with transcripts/genes as keys and GO terms as values.
    """
    transc_to_go = defaultdict(str)
    gene_to_go = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", "Go_final.txt.gz")
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                go = data[1].rstrip()
                transc_to_go[transcript.upper()] = go
                gene_to_go[gene.upper()] = go
    return transc_to_go, gene_to_go

def parse_gene_descriptions():
    """
    Parse the gene descriptions file and return annotation dictionaries.

    Returns:
        tuple: Dictionaries mapping genes and transcripts to functions,
               descriptions, and types.
    """
    gene_to_function = defaultdict(str)
    gene_to_des = defaultdict(str)
    gene_to_type = defaultdict(str)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", "gene_description_20131231.txt.gz")
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                if len(data) < 3:
                    continue
                transcript = data[0].strip()
                gene = transcript.split(".")[0]
                gene_type = data[1].strip()
                des = data[2].strip()
                g_function = "  ".join(data[3:])
                # Populate dictionaries
                gene_to_function[gene.upper()] = g_function
                gene_to_des[gene.upper()] = des
                gene_to_type[gene.upper()] = gene_type
    return gene_to_function, gene_to_des, gene_to_type

def NLR():
    """
    Parse the 'nlrs.txt' file to build a dictionary of NLR genes and their types.

    Returns:
        dict: A dictionary with gene names as keys and NLR types as values.
    """
    gene_to_NLR = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "data", "nlrs.txt")

    with open(input_file, "r") as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                NLR_type = data[1].rstrip()
                gene_to_NLR[gene.upper()] = NLR_type

    return gene_to_NLR

def PTI():
    """
    Parse the PTI_up.txt file to get genes with PTI gene response (upregulated).

    Returns:
        dict: A dictionary with gene IDs as keys and PTI responses as values.
    """
    gene_to_PTI = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "data", "PTI_up.txt")

    with open(input_file, "r") as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                PTI_type = data[1].rstrip()
                gene_to_PTI[gene.upper()] = PTI_type

    return gene_to_PTI

def PTI_down():
    """
    Parse the PTI_down.txt file to get genes with PTI gene response (downregulated).

    Returns:
        dict: A dictionary with gene IDs as keys and PTI responses as values.
    """
    gene_to_PTI = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "data", "PTI_down.txt")

    with open(input_file, "r") as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                PTI_type = data[1].rstrip()
                gene_to_PTI[gene.upper()] = PTI_type

    return gene_to_PTI

def ada6(infile, text_info):
    """
    Parse a file to retrieve gene data for adc6 mutants.

    Args:
        infile (str): Path to the input file.
        text_info (str): Description to associate with each gene.

    Returns:
        dict: A dictionary mapping gene IDs to the associated text info.
    """
    gene_data = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "data", infile)

    with open(input_file, "r") as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                gene_data[gene.upper()] = text_info
                gene_data[transcript.upper()] = text_info

    return gene_data

def parse_names():
    """
    Parse headers.gz file to map transcripts to genes, gene IDs, and functions.

    Returns:
        tuple: Dictionaries for transcript-to-gene, gene-to-transcript, and related mappings.
    """
    trans_to_gene = defaultdict(str)
    gene_to_trans = defaultdict(list)
    tran_to_func = defaultdict(str)
    gene_to_func = defaultdict(str)
    trans_to_gene_id = defaultdict(str)
    gene_to_gene_id = defaultdict(str)
    gene_id_to_func = defaultdict(str)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", "headers.gz")

    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("|")
                transcript = data[0].strip().replace(">", "")
                gene = transcript.split(".")[0]
                gene_id = data[1].split("Symbols: ")[1].strip() if "Symbols: " in data[1] else ""
                annot_func = data[2].strip() if len(data) > 2 else ""

                trans_to_gene[transcript] = gene
                gene_to_trans[gene].append(transcript)
                tran_to_func[transcript] = annot_func
                tran_to_func[gene] = annot_func
                gene_to_func[gene] = annot_func

                trans_to_gene_id[transcript] = gene_id
                trans_to_gene_id[gene] = gene_id
                gene_to_gene_id[gene] = gene_id
                gene_id_to_func[gene_id] = annot_func

    return trans_to_gene, gene_to_trans, tran_to_func, gene_to_func, trans_to_gene_id, gene_to_gene_id, gene_id_to_func

def parse_file(input_file, output_file, annotations):
    """
    Parse the input file, annotate data, and save as an Excel file.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to save the processed Excel file.
        annotations (dict): A dictionary containing parsed annotation data.

    Returns:
        None
    """
    hormone_terms = ["jasmonic", "salicylic", "hormone"]
    pathogen_terms = ["defence", "defense", "pathogenesis", "pathogen", "leucine rich",
                      "immunity"]

    if input_file.endswith('.gz'):
        with gzip.open(input_file, 'rt') as f:
            lines = [line.strip() for line in f if test_line(line)]
    else:
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f if test_line(line)]

    data = [line.split('\t') for line in lines]

    header = data[0]
    if header[0] == "sampleA":
        header.insert(0, "Gene")
        for row in data[1:]:
            row.insert(0, row[0])

    annotation_columns = [
        "Gene", "Gene_ID", "Flowering", "Type", "NLR_Type", "PTI_up", "PTI_down", "acd6", "Function",
        "Description", "GO_Terms", "gene"
    ]
    header = annotation_columns + header[1:]

    annotated_data = []
    for row in data[1:]:
        gene = row[0].upper()
        temp_line = (
            annotations["gene_to_function"].get(gene, "").lower() + " " +
            annotations["gene_to_des"].get(gene, "").lower()
        )
        temp_line = temp_line.lower()

        gene_custom_class = "other"
        for term in hormone_terms:
            if term in temp_line:
                gene_custom_class = "Jasm_sali"

        for term in pathogen_terms:
            if term in temp_line:
                gene_custom_class = "defence"

        gene_id = annotations["trans_to_gene_id"].get(gene, gene)

        annotations_row = [
            gene,
            gene_id,
            "Yes" if gene in annotations["flowering_genes"] else "No",
            gene_custom_class,
            annotations["gene_to_NLR"].get(gene, ""),
            annotations["gene_to_PTI_up"].get(gene, ""),
            annotations["gene_to_PTI_down"].get(gene, ""),
            annotations["adc6_up"].get(gene, "") or annotations["acd6_down"].get(gene, ""),
            annotations["gene_to_function"].get(gene, ""),
            annotations["gene_to_des"].get(gene, ""),
            annotations["gene_to_go"].get(gene, "")
        ]

        annotated_row = annotations_row + row[1:]

        if len(annotated_row) != len(header):
            if len(annotated_row) > len(header):
                annotated_row = annotated_row[:len(header)]
            else:
                annotated_row += [""] * (len(header) - len(annotated_row))

        annotated_data.append(annotated_row)

    df = pd.DataFrame(annotated_data, columns=header)

    output_excel = output_file
    if output_excel.endswith(".txt"):
        output_excel = output_excel.replace(".txt", "_RENAMED.xlsx")
    elif output_excel.endswith(".gz"):
        output_excel = output_excel.replace(".gz", "_RENAMED.xlsx")
    elif output_excel.endswith(".subset"):
        output_excel = output_excel.replace(".subset", ".xlsx")
    else:
        output_excel = f"{output_excel}.xlsx"

    df.to_excel(output_excel, index=False)
    # Save as tab-separated .txt file
    df.to_csv(output_excel.replace(".xlsx", ".txt"), sep="\t", index=False)

    print(f"Processed file saved to: {output_excel}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))

    flowering_genes = parse_flowering_gene()
    transc_to_go, gene_to_go = parse_go()
    gene_to_function, gene_to_des, gene_to_type = parse_gene_descriptions()
    gene_to_NLR = NLR()
    gene_to_PTI_up = PTI()
    gene_to_PTI_down = PTI_down()
    adc6_up = ada6(
        os.path.join(script_dir, "data", "At_adc6_genes.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.acd6-UP.subset"),
        "adc6_up_vs_Col0LFC2_FDR0.001"
    )
    acd6_down = ada6(
        os.path.join(script_dir, "data", "At_adc6_genes.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.col0-UP.subset"),
        "Col0_up_vs_acd6_LFC2_FDR0.001"
    )

    trans_to_gene, gene_to_trans, tran_to_func, gene_to_func, trans_to_gene_id, gene_to_gene_id, gene_id_to_func = parse_names()

    annotations = {
        "flowering_genes": flowering_genes,
        "transc_to_go": transc_to_go,
        "gene_to_go": gene_to_go,
        "gene_to_function": gene_to_function,
        "gene_to_des": gene_to_des,
        "gene_to_type": gene_to_type,
        "gene_to_NLR": gene_to_NLR,
        "gene_to_PTI_up": gene_to_PTI_up,
        "gene_to_PTI_down": gene_to_PTI_down,
        "adc6_up": adc6_up,
        "acd6_down": acd6_down,
        "trans_to_gene": trans_to_gene,
        "trans_to_gene_id": trans_to_gene_id,
    }

    current_directory = os.getcwd()
    for dirpath, _, files in os.walk(current_directory):
        for file in files:
            if file.endswith(".subset"):
                input_file = os.path.join(dirpath, file)
                output_file = os.path.join(dirpath, f"{file}_RENAMED")
                print(f"Processing file: {input_file}")
                parse_file(input_file, output_file, annotations)
