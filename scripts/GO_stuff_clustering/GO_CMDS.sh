#!/bin/bash
#SBATCH -J GO   #jobname
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB

cd /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/At.gene.TMM.EXPR.matrix.RData.clusters_fixed_P_20
/home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl --GO_assignments /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/Go_final.txt --lengths /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/gene_len_GO.txt --genes_single_factor subcluster_100_log2_medianCentered_fpkm.matrix.names --background /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/At.gene.counts.min10.matrix
/home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl --GO_assignments /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/Go_final.txt --lengths /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/gene_len_GO.txt --genes_single_factor subcluster_101_log2_medianCentered_fpkm.matrix.names --background /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/At.gene.counts.min10.matrix
cd /home/pthorpe/scratch/DE_test/at/vir1_temp/all_defence/defence.all.TMM.RData.clusters_fixed_P_30
/home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl --GO_assignments /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/Go_final.txt --lengths /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/gene_len_GO.txt --genes_single_factor subcluster_1_log2_medianCentered_fpkm.matrix.names --background /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp/At.gene.counts.min10.matrix
