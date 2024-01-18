#!/bin/sh
#SBATCH -J DE_TPM   #jobname
#SBATCH --partition=medium
#SBATCH --mem=10GB



cd /mnt/shared/scratch/pthorpe/private/DE_test/at/transcripts
conda activate R
de_files=*0.01

for de in ${de_files}
do
    cut -f1 ${de}| grep -v "Row.names" | sort -u  > ${de}.names
    /home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl --GO_assignments Go_final.txt --lengths gene_len_GO.txt  --genes_single_factor ${de}.names --background At.isoform.counts.min10.matrix
done 

