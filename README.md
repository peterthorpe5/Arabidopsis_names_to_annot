# Arabidopsis_names_to_annot
DE names to gene_ID plus full annotation

Given a differential expression output, the Arabisopsis names e.g. AT1G09665
Are, well, not enough to form any biological meaning. 

```
├── LICENSE
├── README.md
├── converts_names.py
└── data
    ├── Go_final.txt.gz  (name to GO terms)
    ├── flowering_time_genes.gz  (gene involved with flowring times)
    ├── gene.model.names (gene model names, not used yet)
    ├── gene_description_20131231.txt.gz  (info on what the genes do)
    ├── headers.gz  (harvested from the predicted cds file)
    └── nlrs.txt  (from A Species-Wide Inventory of NLR Genes and Alleles in Arabidopsis thaliana )
```

The script will currently walk through the current directory and all sub directories 
and look for files ending with" x.endswith("subset") or x.endswith("FDR_0.001") or x.endswith("DE_results"):

This is currently how my files are set up. The is easy to change. 

## usage
`python converts_names.py`  

This will print the file name it is abusing and write the same filename with "_RENAMED" at the end. 

Enjoy. 


## input example

```
Row.names	sampleA	sampleB	logFC	logCPM	F	PValue	FDR	mutant_high_1	mutant_high_2	mutant_high_3	mutant_high_4	mutant_low_1	mutant_low_2	mutant_low_3	mutant_low_4	col0_high_1	col0_high_2	col0_high_3	col0_high_4	col0_low_1	col0_low_2	col0_low_3	col0_low_4	vir1_high_1	vir1_high_2	vir1_high_3	vir1_high_4	vir1_low_1	vir1_low_2	vir1_low_3	vir1_low_4
AT4G34190	mutant_high	mutant_low	-13.3534107594459	3.88274825009479	22.5486604709497	0.000264141299520225	0.000925229046147501	0	0	0	0	562.19	1251.97	1505.16	1753.15	1177.27	396.95	914.58	735.25	1776.69	2671.63	787.05	858.91	0	0	0	0	976.43	0	0	0
AT3G09890	mutant_high	mutant_low	-9.90230553221601	-0.0186385434096553	41.3203417869738	1.17957343847607e-05	6.38383172555651e-05	0	0	0	0	64.3	187.08	126.6	90.46	7.15	0	0	0	25.81	48.69	169.25	240.7	0	0	0	0	5.69	2.91	8.87	0
AT1G63590
```

# output example

```
gene	geneID	gene_class	R_Gene   	sampleA	sampleB	logFC	logCPM	F	PValue	FDR	mutant_high_1	mutant_high_2	mutant_high_3	mutant_high_4	mutant_low_1	mutant_low_2	mutant_low_3	mutant_low_4	col0_high_1	col0_high_2	col0_high_3	col0_high_4	col0_low_1	col0_low_2	col0_low_3	col0_low_4	vir1_high_1	vir1_high_2	vir1_high_3	vir1_high_4	vir1_low_1	vir1_low_2	vir1_low_3	vir1_low_4	annot	full_annot	GO_terms
AT4G34190	SEP1	other	 AT4G34190	SEP1	other	 stress enhanced protein 1 	mutant_high	mutant_low	-13.3534107594459	3.88274825009479	22.5486604709497	0.000264141299520225	0.000925229046147501	0	0	0	0	562.19	1251.97	1505.16	1753.15	1177.27	396.95	914.58	735.25	1776.69	2671.63	787.05	858.91	0	0	0	0	976.43	0	0	0	stress enhanced protein 1	Encodes a stress enhanced protein that localizes to the thylakoid membrane and whose mRNA is upregulated in response to high light intensity.  It may be involved in chlorophyll binding.  stress enhanced protein 1 (SEP1); Has 30201 Blast hits to 17322 proteins in 780 species: Archae - 12; Bacteria - 1396; Metazoa - 17338; Fungi - 3422; Plants - 5037; Viruses - 0; Other Eukaryotes - 2996 (source: NCBI BLink).	GO:0009535,GO:0009535,GO:0009611,GO:0009644,GO:0016168,GO:0071486,GO:0071492


```