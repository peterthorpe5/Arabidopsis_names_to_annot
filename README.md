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
