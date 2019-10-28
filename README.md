# TF-TWAS Workflow 


# Requirements

R version 3.5.0

glmnet 2.0-18

# Installation

# Run 

Rscript runBaselineModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./ --chr 22

Rscript runTFModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./ --chr 22 --model tf-binding

# Directory Structure and Naming Convention

.
├── script
├── data                            # Test files
│   ├── expression          
|   |   ├── prefix.expression.txt   # gene expression (gene*sample) 
|   |   ├── prefix.expression.RDS   # gene expression (sample*gene, use to convert RDS)
│   ├── genotype
|   |   ├── prefix.chr1.snp.txt     # genotype (varID*sample)
|   |   ├── prefix.chr2.snp.txt
|   |   ├──....
|   |   ├── prefix.chr22.snp.txt
│   └── annotation
|   |   ├──
└──...


# Reference


TF-TWAS is a set of tools to incorporate transcription factors (TFs) into gene imputation models (e.g., PrediXcan). By using TF-TWAS __run_elasticnetCV__ utility, users can construct a gene imputation model with predictor variables that include both cis-SNPs within 1MB of gene and trans-SNPs within 1MB of TFs. On the contrary, the TF-TWAS __run_elasticnetCV_TF_window_zero__ function is designed to use trans-SNPs only within the TF coding region along with the cis-SNPs within 1MB of gene.

## Installation

1. Copy the TF-TWAS folder to a working directory (e.g., _./pkg_).

2. Open the _./pkg/TF-TWAS/script/setup.sh_ and edit the path accordingly. 

3. Set up the enviroment for running TF-TWAS by typing:
```bash 
        bash ./pkg/TF-TWAS/script/tftwas.sh setup ./pkg/TF-TWAS/script/setup.sh
```
4. Review the Inputs section below for file format required for TF-TWAS. 

__For impatients__: refer to example/configs for three example configuration files. The following commands show how to construct TF-both model, TF-binding model and TF-regulation model using TF-TWAS.

```bash

$ bash tftwas.sh run_elasticnetCV_TF config_tf-both.R chr22
$ bash tftwas.sh run_elasticnetCV_TF_window_zero config_tf-binding.R chr22
$ bash tftwas.sh run_elasticnetCV_TF_window_zero config_tf-regulation.R chr22

```


__[Note]__
TF-TWAS is a collection of R scripts unified by the main BASH shell script called tftwas.sh, but there are utilities require Python such as generate_x_matrix, generate_y_matrix and generate_gene_annot. A better computing environment for running TF-TWAS is: __GNU bash, version 4.2.46(2)-release__ (x86_64-redhat-linux-gnu), __R version 3.4.3__ (2017-11-30) -- "Kite-Eating Tree" and __Python 3.6.0__.

## Dependicies

R -- glmnet

```R
> install.packages("glmnet", dependencies=TRUE)
```
Python -- h5dy, scikit-allel

```bash
$ conda install -c conda-forge scikit-allel
$ conda install -c anaconda h5py 
```

## More utilities

Type __bash tftwas.sh help__ to get the list of utilities available in TF-TWAS.

[USAGE] bash tftwas.sh run_elasticnetCV config.R prefix.str

[USAGE] bash tftwas.sh run_elasticnetCV_parallel config.R prefix.str splits.int

[USAGE] bash tftwas.sh run_elasticnetCV_TF_window_zero config.R prefix.str

[USAGE] bash tftwas.sh run_tf_background config.R prefix.str num_times

[USAGE] bash tftwas.sh run_tf_background_window_zero config.R prefix.str num_times

[USAGE] bash tftwas.sh generate_x_matrix expression.txt genotype.vcf output_path.str

[USAGE] bash tftwas.sh generate_y_matrix expression.txt genotype.vcf output_path.str

[USAGE] bash tftwas.sh generate_gene_annot grc37.pos.txt expression.genelist.txt prefix.str output_path.str


## Inputs

TF-TWAS uses a configuration file to locate all the required files described below. 

* __expression_RDS__: file name string; n by m matrix stored in RDS format, where n is gene id and m is sample id.

* __genotype_txt__: file name string; n by m matrix stored in txt format, where n sample id and m is varID.

* __gene_annot_RDS__: file name string; flat file saved in RDS format, columns include gene_id, start, end, genename

* __snp_annot_RDS__: file name string; flat file saved in RDS format, columns include varID

* __tf_gene_annot_RDS__: file name string; flat file saved in RDS format, columns include 

* __tf_genotype_dir__: directory name string; folder where genotype files (in txt format, file name separated by .) for TFs locate

* __tf_snp_annot_dir__: directory name string; folder where snp annotation files (in RDS format, file name separated by .) for TFs locate

* __study__: prefix string for output file (e.g., Test1).

* __snpset__: prefix string for output file (e.g., Illumina).

* __chrom__: prefix string for output file (e.g., 1).

* __out_dir__: directory name string (e.g., ./).

* __n_k_folds__: integer; number of cross validation (e.g., 10).

* __alpha__: float; number of penality for the Elastic-Net (e.g., 0.5)

* __window__: integrer or scientific notation; distance from gene (e.g., 1e6 for 1MB)

* __seed__: integer; number of seed for robust modeling (e.g., 1345)


## Outputs

There are four output files, each of which has output prefix combined with string of "study" and chr"chrom" or "snpset". We used study="TFTWAS", chr"chrom"="chr22" and snpset="Illumina_chr22" as an example: 

1. TFTWAS_chr22_elasticNet_model_log.txt: the log file for a job run, columns include __chr__, __n_genes__, __seed_for_cv__, __alpha__

2. TFTWAS_chr22_snpset_Illumina_chr22_alpha_0.5_covariances.txt: the pairwise correlation of SNPs, columns include __gene_id__, __rsID_1__ ,__rsID_2__, __correlation__

3. TW_TFTWAS_elasticNet_alpha0.5_Illumina_chr22_weights_chr22.txt: coefficients for non-zero SNPs, columns include __gene__, __rsid__, __ref__, __alt__, __beta__, __alpha__

4. working_TW_TFTWAS_exp_10-foldCV_elasticNet_alpha0.5_Illumina_chr22_chr22.txt: gene-based modeling result, columns includes __gene__, __alpha__, __cvm__, __lambda.iteration__, __lambda.min__, __n.snps__, __R2__, __pval__, __genename__, __w_TF__, __pct_TF__, __total.snps__, __used_TF__



