###
#
#
###

#settings -- files
expression_RDS='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/expression/Brain_Hippocampus_Analysis.v6p.normalized.expression.corrected.sharedsamples.RDS'
genotype_txt='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/genotype/Illumina.chr22.pass.maf001rsq08.txt.recode.sharedsamples.txt'
gene_annot_RDS='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/annotation/gene_annot/chr22_gene_pos_GRCh37.p13_mart_export.gene_annot.brain_hippocampus.RDS'
snp_annot_RDS='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/annotation/snp_annot/Illumina.chr22.pass.maf001rsq08.txt.rsids.snp_annot.RDS'   
tf_gene_annot_RDS='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/annotation/tf_gene_annot/gene-tf-pair.big_table.GRC37.positions.brain_hippocampus.genes.RDS'
tf_genotype_dir='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/genotype/'
tf_snp_annot_dir='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/inputs/annotation/tf_snp_annot/exonicSNPs/'

#settings -- study
study='TF-binding'
snpset='Illumina_chr22'
chrom='22'
out_dir='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/example/outputs/'

#settings -- model
n_k_folds=as.numeric('10')
alpha=as.numeric('0.5')
window=as.numeric(1e6)
seed=as.numeric(1345)
