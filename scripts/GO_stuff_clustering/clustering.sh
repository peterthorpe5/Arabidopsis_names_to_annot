#!/bin/bash


conda activate R

~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/PtR -s samples_described.txt -m At.gene.TMM.EXPR.matrix --heatmap_scale_limits "-2,2" --min_rowSums 20 --heatmap_scale_limits "-2,2" --save --heatmap --log2 --center_rows --imgfmt svg --img_width 8 --img_height 8


~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 20 -R At.gene.TMM.EXPR.matrixRData


mv clusters_fixed_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_P_20
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 30 -R At.gene.TMM.EXPR.matrixRData

mv clusters_fixed_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_P_30
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering --Ptree 50 -R At.gene.TMM.EXPR.matrixRData


mv clusters_fixed_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_P_50
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ptree 70 -R At.gene.TMM.EXPR.matrixRData
mv clusters_fixed_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_P_70
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ktree 4 -R At.gene.TMM.EXPR.matrixRData

mv clusters_fixed_Ktree_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_Ktree_4
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ktree 2 -R At.gene.TMM.EXPR.matrixRData

mv clusters_fixed_Ktree_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_Ktree_2
~/apps/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --lexical_column_ordering  --Ktree 8 -R At.gene.TMM.EXPR.matrixRData

mv clusters_fixed_Ktree_* ./At.gene.TMM.EXPR.matrixRData.clusters_fixed_Ktree_8



