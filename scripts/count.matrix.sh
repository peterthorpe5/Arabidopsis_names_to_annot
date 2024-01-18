

conda activate R


trin_path=/homes/pthorpe001/apps/trinityrnaseq-v2.15.1/
v9=/homes/pthorpe001/apps/trinityrnaseq-v2.15.1/


# kallisto
perl ${trin_path}/util/abundance_estimates_to_matrix.pl --est_method kallisto --gene_trans_map  none --name_sample_by_basedir --out_prefix At ./*/abundance.tsv



# saloon
perl ${trin_path}/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map  none --name_sample_by_basedir --out_prefix At ./*/quant.sf
