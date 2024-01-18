#!/bin/bash
#SBATCH -J GO   #jobname
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=14GB
set -e

cd  /mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp

data=/mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp
folder=At.gene.TMM.EXPR.matrix.RData.clusters_fixed_P_*
wd=/mnt/shared/scratch/pthorpe/private/DE_test/at/vir1_temp


for f in ${folder}
do 
    cd "${wd}/${f}" || exit
    echo "cd ${wd}/${f}"
    files=*_log2_medianCentered_fpkm.matrix
    for input_file in ${files}
    do
        # go command -  remove first line and extract first col
        awk 'NR > 1 { print $1 }' "${input_file}" > ${wd}/${f}/${input_file}.names
        wait
        cmd="/home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl 
            --GO_assignments ${data}/Go_final.txt --lengths 
            ${data}/gene_len_GO.txt 
            --genes_single_factor  
            ${input_file}.names
            --background  
            ${data}/At.gene.counts.min10.matrix"
        echo ${cmd}
        eval ${cmd}
        wait
        #rm temp_file.txtc
    done
done    


# defence only

cd /home/pthorpe/scratch/DE_test/at/vir1_temp/all_defence
folder=defence.all.TMM.RData.clusters_fixed_P_*

wd=/home/pthorpe/scratch/DE_test/at/vir1_temp/all_defence

for f in ${folder}
do 
    cd "${wd}/${f}" || exit
    echo "cd ${wd}/${f}"
    files=*_log2_medianCentered_fpkm.matrix
    for input_file in ${files}
    do
        # go command -  remove first line and extract first col
        awk 'NR > 1 { print $1 }' "${input_file}" > ${wd}/${f}/${input_file}.names
        cmd="/home/pthorpe/apps/trinityrnaseq-v2.13.2/Analysis/DifferentialExpression/run_GOseq.pl 
            --GO_assignments ${data}/Go_final.txt --lengths 
        ${data}/gene_len_GO.txt 
        --genes_single_factor  
        ${input_file}.names
        --background  
        ${data}/At.gene.counts.min10.matrix"
        echo ${cmd}
        eval ${cmd}
        #rm temp_file.txt
    done
done 


