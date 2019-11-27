# TF-TWAS Workflow 

![Workflow](/Figures/Figure2.png)

# Requirements

R       3.5.0

python  3.6.8

glmnet 2.0-18 (R-version)

# Installation

# Run 

* PartA: identifying the outliers   

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

* PartB: testing randomness for outliers

```R
#step1: run background model by randomly selecting TFs (change --chr to run other chromosomes for the outliers)
Rscript ./script/runBackgroundModel.R --scriptPth ./script/ --filePath ./data/  --filePrefix geuvadis --outFolder ./TF-binding/backgroundmodel/ -c 11 -m tf-binding -r 100

Rscript ./script/runBackgroundModel.R --scriptPth ./script/ --filePath ./data/  --filePrefix geuvadis --outFolder ./TF-regulation/backgroundmodel/ -c 1 -m tf-regulation -r 100

Rscript ./script/runBackgroundModel.R --scriptPth ./script/ --filePath ./data/  --filePrefix geuvadis --outFolder ./TF-both/backgroundmodel/ -c 21 -m tf-both -r 100

#step2: parse background model results
python ./script/parse_backgroundmodel_result.py --outlierFile ./results/geuvadis.tf-binding.43.outliers.txt --modelResultFile ./results/geuvadis.tf-binding.result.txt --bgfolder_path ./TF-binding/backgroundmodel/ --model_name tf-binding --out_prefix geuvadis

 python ./script/parse_backgroundmodel_result.py --outlierFile ./results/geuvadis.tf-regulation.14.outliers.txt --modelResultFile   ./results/geuvadis.tf-regulation.result.txt --bgfolder_path ./TF-regulation/backgroundmodel/ --model_name tf-regulation --out_prefix geuvadis
 
 python ./script/parse_backgroundmodel_result.py --outlierFile ./results/geuvadis.tf-both.3.outliers.txt --modelResultFile ./results/geuvadis.tf-both.result.txt --bgfolder_path ./TF-both/backgroundmodel/ --model_name tf-both --out_prefix geuvadis
```

* PartC: testing roubustness of models by random selecting samples

```R
#step1. run random fold selection for baseline model and for TF models

Rscript ./script/runRandomFoldSelection.R --scriptPth ./script/ --filePath ./data/ --filePrefix geuvadis --outFolder ./Baseline/randomFoldSelection/ -c 21 -m baseline -r 100

Rscript ./script/runRandomFoldSelection.R --scriptPth ./script/ --filePath ./data/ --filePrefix geuvadis --outFolder ./TF-binding/randomFoldSelection/ -c 21 -m tf-bining -r 100

Rscript ./script/runRandomFoldSelection.R --scriptPth ./script/ --filePath ./data/ --filePrefix geuvadis --outFolder ./TF-regulation/randomFoldSelection/ -c 21 -m tf-regulation -r 100

Rscript ./script/runRandomFoldSelection.R --scriptPth ./script/ --filePath ./data/ --filePrefix geuvadis --outFolder ./TF-both/randomFoldSelection/ -c 21 -m tf-both -r 100

#setp2. parse random selection results

python ../script/parse_randomSelection_result.py --outlierFile ./results/tf-binding.outliers.pass.backgroundmodel.txt 
--baseline_resultFolder ./Baseline/randomFoldSelection --tfmodel_resultFolder ./TF-binding/randomFoldSelection --model_name tf-binding
--threshold 0.01 --run 100 --out_prefix tf-binding

python ../script/parse_randomSelection_result.py --outlierFile ./results/tf-regulation.outliers.pass.backgroundmodel.txt 
--baseline_resultFolder ./Baseline/randomFoldSelection --tfmodel_resultFolder ./TF-regulation/randomFoldSelection --model_name tf-regulation --threshold 0.01 --run 100 --out_prefix tf-binding

python ../script/parse_randomSelection_result.py --outlierFile ./results/tf-both.outliers.pass.backgroundmodel.txt 
--baseline_resultFolder ./Baseline/randomFoldSelection --tfmodel_resultFolder ./TF-both/randomFoldSelection --model_name tf-both
--threshold 0.01 --run 100 --out_prefix tf-both
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







