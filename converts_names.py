""" Title: Convert transcript or gene to gene name
"""

# imports

import os
from sys import stdin,argv
import sys
from collections import defaultdict
import gzip



def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def parse_flowering_gene():
    """get the list of flowering time genes
    from: https://www.mpipz.mpg.de/14637/Arabidopsis_flowering_genes 
    data october 2023
    return a set of genes"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", 
                             "flowering_time_genes.gz")
    flowering_genes = set([])
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                gene = line.strip()
                flowering_genes.add(gene)
    return flowering_genes


def NLR():
    """returns a dict [gene] NLR type 	
    e.g.
    At1g57650	NBS-LRR
    At1g17615	TIR-NBS
    """
    gene_to_NLR = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    f_in = open(os.path.join(script_dir, "data", 
                             "nlrs.txt"), "r")
    for line in f_in:
        if test_line(line):
            if line.startswith("sampleA"): continue
            data = line.split("\t")
            transcript = data[0].strip()
            gene = transcript.split(".")[0].rstrip()
            NBR_type = data[1].rstrip()
            gene_to_NLR[gene.upper()] = NBR_type
    return gene_to_NLR


def PTI():
    """returns a dict [gene] PTI_gene_response_LFC 	
    e.g.
    AT1G51820	LRR-RLKs_down_PTI_LFC_-5.935654809
    AT2G19190	LRR-RLKs_down_PTI_LFC_-5.404993937

    """
    gene_to_PTI = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    f_in = open(os.path.join(script_dir, "data", 
                             "PTI.txt"), "r")
    for line in f_in:
        if test_line(line):
            if line.startswith("sampleA"): continue
            data = line.split("\t")
            transcript = data[0].strip()
            gene = transcript.split(".")[0].rstrip()
            PTI_type = data[1].rstrip()
            gene_to_PTI[gene.upper()] = PTI_type
    return gene_to_PTI


def ada6(infile, text_info):
    """first coloumn is the gene name
    """
    DE_acd6 = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    f_in = open(os.path.join(script_dir, "data", 
                             infile), "r")
    for line in f_in:
        if test_line(line):
            data = line.split("\t")
            transcript = data[0].strip()
            gene = transcript.split(".")[0].rstrip()
            DE_acd6[gene.upper()] = text_info
            DE_acd6[transcript.upper()] = text_info
    return DE_acd6


def parse_go():
    """we need to convert the gene id to  go terms	
    """
    transc_to_go = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, "data", 
                             "Go_final.txt.gz")
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("\t")
                transcript = data[0].strip()
                gene = transcript.split(".")[0].rstrip()
                go = data[1].rstrip()
                transc_to_go[transcript] = go
                transc_to_go[gene] = go
    return transc_to_go


def parse_gene_des():
    """we need to convert the gene id to  function:
    EXAMPLE:
    AT1G14820.1	protein_coding	Sec14p-like phosphatidylinositol transfer family protein		Sec14p-like phosphatidylinositol transfer family protein; CONTAINS InterPro DOMAIN/s: Cellular retinaldehyde-binding/triple function, C-terminal (InterPro:IPR001251), Phosphatidylinositol transfer protein-like, N-terminal (InterPro:IPR011074); BEST Arabidopsis thaliana protein match is: Sec14p-like phosphatidylinositol transfer family protein (TAIR:AT1G01630.1); Has 1699 Blast hits to 1699 proteins in 217 species: Archae - 0; Bacteria - 0; Metazoa - 335; Fungi - 476; Plants - 691; Viruses - 0; Other Eukaryotes - 197 (source: NCBI BLink).
    AT1G14690.1	protein_coding	microtubule-associated protein 65-7		microtubule-associated protein 65-7 (MAP65-7); CONTAINS InterPro DOMAIN/s: Microtubule-associated protein, MAP65/ASE1-type (InterPro:IPR007145); BEST Arabidopsis thaliana protein match is: Microtubule associated protein (MAP65/ASE1) family protein (TAIR:AT2G01910.1); Has 3902 Blast hits to 3379 proteins in 390 species: Archae - 69; Bacteria - 268; Metazoa - 1928; Fungi - 367; Plants - 438; Viruses - 6; Other Eukaryotes - 826 (source: NCBI BLink).
    AT1G14700.1	protein_coding	purple acid phosphatase 3	
    """
    gene_to_function = defaultdict(str)
    gene_to_des = defaultdict(str)
    gene_to_type = defaultdict(str)
    transc_to_function = defaultdict(str)
    transc_to_type = defaultdict(str)
    transc_to_des = defaultdict(str)
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    infile = os.path.join(script_dir, 
                             "data", 
                             "gene_description_20131231.txt.gz")
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
                # populate the dictionaries
                gene_to_function[gene] = g_function
                gene_to_des[gene] = des
                gene_to_type[gene] = gene_type
                transc_to_function[transcript] = g_function
                transc_to_function[gene] = g_function
                transc_to_type[transcript] = gene_type
                transc_to_des[transcript] = des
    return gene_to_function, gene_to_des, gene_to_type, \
        transc_to_function, transc_to_type, transc_to_des


def parse_names():
    """we need to convert the trans or gene name to the func name
    >ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|
    OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|

    reruns to default dict e.g.
    T  = transcript
    G = gene
    [ENST] = DDX11L1
    [ENSG] = DDX11L1

    Arabidopsis example:
    >AT1G01020.3 | Symbols: ARV1 |  | chr1:6915-8442 REVERSE LENGTH=236
    >AT1G01060.2 | Symbols: LHY1, LHY | LATE ELONGATED HYPOCOTYL, LATE ELONGATED HYPOCOTYL 1 | 
    chr1:33992-37061 REVERSE LENGTH=64
    5
    """
    trans_to_gene = defaultdict(str)
    gene_to_trans = defaultdict(list)
    tran_to_func = defaultdict(str)
    gene_to_func = defaultdict(str)
    #
    trans_to_gene_id = defaultdict(str)
    gene_to_gene_id = defaultdict(str)
    gene_id_to_func = defaultdict(str)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    infile = os.path.join(script_dir, 
                             "data", "headers.gz")
    with gzip.open(infile, 'rt') as f_in:
        for line in f_in:
            if test_line(line):
                data = line.split("|")
                transcript = data[0].strip()
                transcript = transcript.replace(">", "")
                gene_id = data[1].strip()
                gene_id = gene_id.split("Symbols: ")[1].strip()
                annot_func = data[2].strip()
                position = data[3].strip()
                # populate the dictionaries
                gene = transcript.split(".")[0]
                trans_to_gene[transcript] = gene
                gene_to_trans[gene].append(transcript)
                tran_to_func[transcript] = annot_func
                tran_to_func[gene] = annot_func
                gene_to_func[gene] = annot_func
            #            
                trans_to_gene_id[transcript] = gene_id
                trans_to_gene_id[gene] = gene_id
                gene_to_gene_id[gene] = gene
                gene_id_to_func[gene] = annot_func
    return trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func
            

def parse_file(infile, outfile, trans_to_gene,
               gene_to_trans, tran_to_func,
               gene_to_func, trans_to_gene_id, 
               gene_to_gene_id, gene_id_to_func, 
               transc_to_function, transc_to_go, 
               flowering_genes, pathogen_terms,
               hormone_terms, gene_to_NLR,
               adc6_up, col0_up, adc6_RTD_up,
               col0_RTD_up, gene_to_PTI):
    """ function to parse the input files and return the info
    for the given gene in coloumn 1"""
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    infile = os.path.basename(infile)
    header_out = 0


    for line in f_in:
        if test_line(line):
            gene_custom_class = "other"
            nbl_type = ""
            transcript = ""
            line = line.rstrip()
            acd6 = ""
            acd6_trans = ""
            PTI = ""
            #if line.startswith("sampleA"): 
            #    if "GLM.edgeR.DE" in infile: pass
            #    if "DE_results" in infile: pass
            #    f_out.write("\t" + "gene_id" + "\t" + line)
            #    continue
            data = line.split("\t")
            subject = data[0]
            if "." in subject:
                transcript = subject
                subject = subject.split(".")[0]
            # the next two if are for trinity DE output. 
            # removed
            #if "transcript" in infile:
            #    if "GLM.edgeR.DE" in infile: pass
            #    if "DE_results" in infile: pass
            #    gene_id = trans_to_gene_id[subject]
            #    func = tran_to_func[subject]
            #    out_data = "%s\t%s %s " % (subject, gene_id, func)
            #    line = line.replace(subject, out_data)
            #    f_out.write(line)
                
            #if "gene" in infile:
            #    if "GLM.edgeR.DE" in infile: pass
            #    if "DE_results" in infile: pass
            #    gene_id = gene_to_gene_id[subject]
            #    func = gene_to_func[subject]
            #    out_data = "%s %s %s " % (subject, gene_id, func)
            #    line = line.replace(subject, out_data)
            #    f_out.write(line)

            # write the header for the DE format
            if line.strip().startswith("Row.names"):
                if header_out < 1:
                    if transcript != "":
                        line = line.replace("Row.names", "transcript\tgene\tgeneID\tgene_class\tR_Gene\tacd6_mutants\tacd6_RTD")
                    else:
                        line = line.replace("Row.names", "gene\tgeneID\tgene_class\tR_Gene\tPTI_type\tacd6_mutants")
                    header = line.rstrip() + "\t" + "annot" + "\t" + "full_annot" + "\t" + "GO_terms" +"\n"
                    f_out.write(header)
                    header_out = header_out + 1
                continue
            
            if "GLM.edgeR.DE" in infile or  "DE_results" in infile:
                if not subject.startswith("A"):
                    subject = data[1]
                gene_id = trans_to_gene_id[subject]
                func = tran_to_func[subject]
                func = func.rstrip("\n")
                full_funk = transc_to_function[subject]
                full_funk = full_funk.rstrip()
                go = transc_to_go[subject]

                temp_line = (line.rstrip() + "\t" + gene_custom_class + func + 
                            "\t" + full_funk.strip() + "\t" +
                            go + "\n")
                temp_line = temp_line.lower()
                if subject in flowering_genes:
                    gene_custom_class = "flowering"
                # see if the hormones in the info
                for term in hormone_terms:
                    if term in temp_line:
                        gene_custom_class = "Jasm_sali"
                for term in pathogen_terms:
                    if term in temp_line:
                        gene_custom_class = "defence"
                if gene_to_NLR[subject]:
                    nbl_type = gene_to_NLR[subject]
                
                if gene_to_PTI[subject]:
                    PTI = gene_to_PTI[subject]
                    
                if adc6_up[subject]:
                    acd6 = adc6_up[subject]
                if col0_up[subject]:
                    acd6 = col0_up[subject]
                if acd6.startswith("0.") or acd6.startswith("."):
                    acd6 = ""
                # some mental if to get rid of some weirdness coming through. 
                if transcript != "":
                    if adc6_RTD_up[transcript]:
                        acd6_trans = adc6_RTD_up[transcript]
                        if not acd6_trans.startswith("a") or acd6_trans.startswith("C"):
                            acd6_trans = "" 
                        if acd6_trans.startswith("0.") or acd6_trans.startswith("."):
                            acd6_trans = "" 
                    if col0_RTD_up[transcript]:
                        acd6_trans = col0_RTD_up[transcript]
                        if not acd6_trans.startswith("a") or acd6_trans.startswith("C"):
                            acd6_trans = "" 
                        if acd6_trans.startswith("0.") or acd6_trans.startswith("."):
                            acd6_trans = "" 
                    if acd6_trans.startswith("0.") or acd6_trans.startswith("."):
                        acd6_trans = ""
                    if not acd6_trans.startswith("a") or acd6_trans.startswith("C"):
                        acd6_trans = ""    
                    acd6 = acd6 + "\t" + acd6_trans


                out_data = "%s\t%s\t%s\t%s\t%s\t%s" % (subject, gene_id.rstrip(),  
                                                      gene_custom_class.rstrip(), 
                                                      nbl_type.rstrip(), PTI.rstrip(), 
                                                      acd6)
                if transcript != "":
                    out_data = transcript + "\t" + out_data
                line = line.replace(subject, out_data).rstrip()
                f_out.write(line.rstrip() + "\t" + func.rstrip() + 
                            "\t" + full_funk.strip() + "\t" +
                            go + "\n")
    f_in.close()
    f_out.close()
    

if __name__ == '__main__':
    # to get where the data folder is ...
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # parses headers which is from the Araport11_pep_20220914
    trans_to_gene, gene_to_trans, tran_to_func, \
           gene_to_func, trans_to_gene_id, \
           gene_to_gene_id, gene_id_to_func = parse_names()
    # parses the file: gene_description_20131231.txt from Araport
    gene_to_function, gene_to_des, gene_to_type, \
        transc_to_function, transc_to_type, \
            transc_to_des = parse_gene_des()
    transc_to_go = parse_go()

    # parses the file: contain R genes from A Species-Wide Inventory of NLR Genes and Alleles in Arabidopsis thaliana (2017)
    gene_to_NLR = NLR()
    flowering_genes = parse_flowering_gene()

    gene_to_PTI = PTI()

    # get the gene up in mutant adc6 and col0
    # Initialize dictionaries for adc6_up and col0_up
    adc6_up = {}
    col0_up = {}
    
    adc6_up = ada6(os.path.join(script_dir, "data", 
                                "At_adc6_genes.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.acd6-UP.subset"), 
                   "adc6_up_vs_Col0LFC2_FDR0.001")
    col0_up = ada6(os.path.join(script_dir, "data",
                                "At_adc6_genes.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.col0-UP.subset"), 
                   "Col0_up_vs_acd6_LFC2_FDR0.001")
    
    # Initialize dictionaries for adc6_up and col0_up at_RTD3_data
    adc6_RTD_up = {}
    col0_RTD_up = {}

    adc6_RTD_up = ada6(os.path.join(script_dir, "data",
                                    "At_adc6_atRTD3.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.acd6-UP.subset"), 
                   "adc6_up_vs_Col0_L2_1e-3")
    col0_RTD_up = ada6(os.path.join(script_dir, "data",
                                    "At_adc6_atRTD3.isoform.counts.matrix.acd6_vs_col0.edgeR.DE_results.P1e-3_C2.col0-UP.subset"), 
                   "Col0_up_vs_acd6_L2_1e-3")
    # call the function to get a list of results wanted
    directory = os.getcwd()
    # get all folder names os.walk(directory)
    for dirpath, subdirs, files in os.walk(directory):
        for x in files:
    
            if x.endswith(".pdf"): continue
            if x.endswith("RENAMED"): continue
            if x.endswith("subset") or x.endswith("FDR_0.001") or x.endswith("DE_results") or x.endswith("FDR_0.01") or x.endswith("FDR_0.05"):
                if x.startswith("."): continue
                wanted_file = (os.path.join(dirpath, x))
                print("FILE = ", wanted_file)
                outfile = wanted_file + "_RENAMED"
                #print(wanted_file)
                pathogen_terms = ["disease", "defence", "defense",
                                  "pathogenesis", "leucine-rich", "(lrr)",
                                  "avirulence"]
                hormone_terms = ["jasmonic",  "salicylic"]
                parse_file(wanted_file, outfile, trans_to_gene,
                           gene_to_trans, tran_to_func,
                           gene_to_func, trans_to_gene_id, 
                           gene_to_gene_id, gene_id_to_func,
                           transc_to_function, transc_to_go,
                           flowering_genes, pathogen_terms,
                           hormone_terms, gene_to_NLR,
                           adc6_up, col0_up, adc6_RTD_up,
                           col0_RTD_up, gene_to_PTI)
