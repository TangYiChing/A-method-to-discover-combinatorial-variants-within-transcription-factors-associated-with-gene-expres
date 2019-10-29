# TF-TWAS Workflow 

![Workflow](/Figures/Figure2.png)

# Requirements

R       3.5.0

python  3.6.8

glmnet 2.0-18 (R-version)

# Installation

# Run PartA of the workflow: identifying the outliers   

```R
# step1. run the baseline model (change --chr to run other chromosomes)
Rscript runBaselineModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./baseline --chr 22

# step2. run the TF models (change --chr to run other chromosomes)
Rscript runTFModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./TF-regulation/ --chr 22 --model tf-regulation
Rscript runTFModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./TF-binding/ --chr 22 --model tf-binding
Rscript runTFModel.R --scriptPth ./script/ --filePath ./data/ --outFolder ./TF-both/ --chr 22 --model tf-both

# step3. parse results
python ./script/parse_result.py --result_path ./baseline/ --model_name baseline --out_prefix geuvadis
python ./script/parse_result.py --result_path ./TF-regulation/ --model_name tf-regulation --out_prefix geuvadis
python ./script/parse_result.py --result_path ./TF-binding/ --model_name tf-binding --out_prefix geuvadis
python ./script/parse_result.py --result_path ./TF-both/ --model_name tf-both --out_prefix geuvadis
[Note] this step will generate result files from each model containing modeling result of chr1-chr22 (which is the output of step1.)

# step4. identify outliers by comparing R2 to the baseline model 
python ./script/identify_outlier.py --tf_resultFile geuvadis.tf-regulation.result.txt \
                                    --bl_resultFile geuvadis.baseline.result.txt \
                                    --out_prefix geuvadis.tf-regulation
python ./script/identify_outlier.py --tf_resultFile geuvadis.tf-binding.result.txt \
                                    --bl_resultFile geuvadis.baseline.result.txt \
                                    --out_prefix geuvadis.tf-binding
python ./script/identify_outlier.py --tf_resultFile geuvadis.tf-both.result.txt \
                                    --bl_resultFile geuvadis.baseline.result.txt \ 
                                    --out_prefix geuvadis.tf-both
```
# Directory Structure and Naming Convention

```bash
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
|   |   ├──snp_annot
|   |   |  ├──prefix.chr1.snp_annot.RDS
|   |   |  ├──prefix.chr2.snp_annot.RDS
|   |   |  ├──....
|   |   |  ├──prefix.chr22.snp_annot.RDS
|   |   ├──gene_annot
|   |   |  ├──prefix.gene_annot.RDS
|   |   ├──tf_gene_annot
|   |   |  ├──prefix.tfgene_annot.RDS
└──...
```

# Reference


TF-TWAS is a set of tools to incorporate transcription factors (TFs) into gene imputation models (e.g., PrediXcan). By using TF-TWAS __run_elasticnetCV__ utility, users can construct a gene imputation model with predictor variables that include both cis-SNPs within 1MB of gene and trans-SNPs within 1MB of TFs. On the contrary, the TF-TWAS __run_elasticnetCV_TF_window_zero__ function is designed to use trans-SNPs only within the TF coding region along with the cis-SNPs within 1MB of gene.



__[Note]__
TF-TWAS is a collection of R scripts unified by the main BASH shell script called tftwas.sh, but there are utilities require Python such as generate_x_matrix, generate_y_matrix and generate_gene_annot. A better computing environment for running TF-TWAS is: __GNU bash, version 4.2.46(2)-release__ (x86_64-redhat-linux-gnu), __R version 3.4.3__ (2017-11-30) -- "Kite-Eating Tree" and __Python 3.6.0__.







