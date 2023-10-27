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
Row.names	sampleA	sampleB	logFC	logCPM	F	PValue	FDR	mutant_high_1	mutant_high_2	mutant_high_3	mutant_high_4	mutant_low_1	mutant_low_2	mutant_low_3	mutant_low_4	col0_high_1	col0_high_2	col0_high_3	col0_high_4	col0_low_1	col0_low_2	col0_low_3	col0_low_4	mutant2_high_1	mutant2_high_2	mutant2_high_3	mutant2_high_4	mutant2_low_1	mutant2_low_2	mutant2_low_3	mutant2_low_4
AT4G34190	mutant_high	mutant_low	-13.3534107594459	3.88274825009479	22.5486604709497	0.000264141299520225	0.000925229046147501	0	0	0	0	562.19	1251.97	1505.16	1753.15	1177.27	396.95	914.58	735.25	1776.69	2671.63	787.05	858.91	0	0	0	0	976.43	0	0	0
AT3G09890	mutant_high	mutant_low	-9.90230553221601	-0.0186385434096553	41.3203417869738	1.17957343847607e-05	6.38383172555651e-05	0	0	0	0	64.3	187.08	126.6	90.46	7.15	0	0	0	25.81	48.69	169.25	240.7	0	0	0	0	5.69	2.91	8.87	0
AT1G63590
```

# output example

```

gene	geneID	gene_class	R_Gene	sampleA	sampleB	logFC	logCPM	F	PValue	FDR	mutant2_high_1	mutant2_high_2	mutant2_high_3	mutant2_high_4	mutant2_low_1	mutant2_low_2	mutant2_low_3	mutant2_low_4	col0_high_1	col0_high_2	col0_high_3	col0_high_4	col0_low_1	col0_low_2	col0_low_3	col0_low_4	mutant1_high_1	mutant1_high_2	mutant1_high_3	mutant1_high_4	mutant1_low_1	mutant1_low_2	mutant1_low_3	mutant1_low_4	annot	full_annot	GO_terms
AT5G58670	PLC1, ATPLC1, ATPLC	other		mutant1_low	AVR_WO_10C	0.052254597	1.119381761	0.175064409	0.679920781	0.718034798	119.88	106.48	136.77	130.4	146.22	105.9	80.88	109.25	97.46	109.41	78.95	77.9	113.6	127.75	85.99	103.52	34.66	41.98	46.62	43.9	62.27	87.08	80.54	92.2	ARABIDOPSIS THALIANA PHOSPHOLIPASE C, phospholipase C1, phospholipase C 1	phosphatidylinositol-specific phospholipase C is induced to a significant extent under various environmental stresses, such as dehydration, salinity, and low temperature. May play a role in secondary ABA response.  There are two genes called ATPLC1, one corresponding to AT4g38530 and one corresponding ot AT5g58670 (this one).  phospholipase C1 (PLC1); CONTAINS InterPro DOMAIN/s: Phospholipase C, phosphoinositol-specific, EF-hand-like (InterPro:IPR015359), Phospholipase C, phosphatidylinositol-specific , X domain (InterPro:IPR000909), EF-Hand 1, calcium-binding site (InterPro:IPR018247), C2 calcium/lipid-binding domain, CaLB (InterPro:IPR008973), Phospholipase C, phosphoinositol-specific (InterPro:IPR001192), Phospholipase C, phosphatidylinositol-specific, Y domain (InterPro:IPR001711), PLC-like phosphodiesterase, TIM beta/alpha-barrel domain (InterPro:IPR017946), C2 membrane targeting protein (InterPro:IPR018029), C2 calcium-dependent membrane targeting (InterPro:IPR000008); BEST Arabidopsis thaliana protein match is: phospholipase C1 (TAIR:AT4G38530.1); Has 30201 Blast hits to 17322 proteins in 780 species: Archae - 12; Bacteria - 1396; Metazoa - 17338; Fungi - 3422; Plants - 5037; Viruses - 0; Other Eukaryotes - 2996 (source: NCBI BLink).	GO:0004435,GO:0005634,GO:0048015,GO:0051209,GO:0004629,GO:0009737,GO:0009738,GO:0004629,GO:0009409,GO:0009414,GO:0009651
AT2G04550	IBR5, DSPTP1E	other		mutant1_low	AVR_WO_10C	0.052375538	4.279841405	0.79088474	0.383962172	0.431765605	842.87	821.66	847.44	865.86	990.86	827.78	777.73	684.63	982.08	991.2	851.69	900.09	909.37	900.82	654.83	796.02	643.31	801.3	641.62	804.89	967.95	913.62	798.45	882.82	DUAL SPECIFICITY PROTEIN PHOSPHATASE 1E, indole-3-butyric acid response 5	Encodes a protein phosphatase that interacts with MPK12, but not with other MAP kinases. It can dephosphorylate a dually phosphorylated MPK12 in vitro and can inactivate MPK12 in vivo. ibr5 mutants have reduced sensitivity to auxin and abscisic acid.  IBR5 promotes auxin responses, including auxin-inducible transcription, differently than the TIR1 auxin receptor and without destabilizing Aux/IAA repressor proteins.  indole-3-butyric acid response 5 (IBR5); CONTAINS InterPro DOMAIN/s: Protein-tyrosine phosphatase, active site (InterPro:IPR016130), Dual-specific/protein-tyrosine phosphatase, conserved region (InterPro:IPR000387), Dual specificity phosphatase, catalytic domain (InterPro:IPR000340), Dual specificity phosphatase, subgroup, catalytic domain (InterPro:IPR020422); BEST Arabidopsis thaliana protein match is: dual specificity protein phosphatase 1 (TAIR:AT3G23610.2); Has 3359 Blast hits to 3359 proteins in 285 species: Archae - 9; Bacteria - 32; Metazoa - 1975; Fungi - 228; Plants - 239; Viruses - 200; Other Eukaryotes - 676 (source: NCBI BLink).	GO:0004725,GO:0005634,GO:0005634,GO:0005634,GO:0006470,GO:0008138,GO:0033549,GO:0005516,GO:0009733,GO:0009737,GO:0005515,GO:0005634,GO:0009734,GO:0033549,GO:0035556,GO:0046620,GO:0061388,GO:0005515,GO:0005634,GO:0005515,GO:0005515,GO:0005515,GO:0005515,GO:0005515,GO:0005515


```