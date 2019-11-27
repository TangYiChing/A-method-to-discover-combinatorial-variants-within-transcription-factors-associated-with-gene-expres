"""

1. set cutoff=Inf and n.max=Inf in the find_outlier.R to retrieve all observations (i.e., without filtering by Bonferroni p-values=0.01)
and obtained files:
    whole_blood.tf-binding.1418.outliers.txt
    whole_blood.tf-regulation.1286.outliers.txt
    whole_blood.tf-both.2628.outliers.txt

2. find geuvadis outliers to see where are they (which ranking) in the whole blood outliers
"""

import sys
import pandas as pd
import subprocess as sp
from gtfparse import read_gtf

RPATH = '~/anaconda3/envs/bmi6319/bin/Rscript'
SCRIPTPATH = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/script/'
PREDICAN = '/data2/GTEx_v6/code/PredictDBPipeline-master/scripts/'

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def find_common_and_ranking(geuvadis, wb, model='tf-binding'):
    """
    :param geuvadis: outlier table, must have headers=[genename]
    :param wb: outlier table, must have headers=[genename]
    :return model, common, ranking: modelName, genename, indexNumber
    """
    # find common 
    common = sorted( list(set(wb['genename']) & set(geuvadis['genename'])) )
    
    # find ranking
    if len(common) == 0:
        print( 'model={:}, common=0'.format(model) )
    else:
       common_in_wb = wb.loc[wb['genename'].isin(common)]
       print(common_in_wb)
       wb.set_index('genename', inplace=True)
       for name in common:
           ranking = wb.index.get_loc(name)
           print( 'model={:}, genename={:}, ranking={:}'.format(model, name, ranking) )

        
if __name__ == "__main__":
    # set args
    debug = True
    outPath = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/'
    # input files
    wbTFbinding = 'whole_blood.tf-binding.1418.outliers.txt'
    wbTFregulation = 'whole_blood.tf-regulation.1286.outliers.txt'
    wbTFboth = 'whole_blood.tf-both.2628.outliers.txt'

    guTFbinding = 'geuvadis.tf-binding.22.outliers.txt'
    guTFregulation = 'geuvadis.tf-regulation.WBeQTL.21.outliers.txt' #'geuvadis.tf-regulation.32.outliers.txt'
    guTFboth = 'geuvadis.tf-both.3.outliers.txt'

    # locate ranking -- TF-binding
    wb = read_as_df(outPath+wbTFbinding)
    gu = read_as_df(outPath+'result/'+guTFbinding)
    find_common_and_ranking(gu, wb, model='tf-binding')
    # locate ranking -- TF-regulation
    wb = read_as_df(outPath+wbTFregulation)
    gu = read_as_df(outPath+'result/'+guTFregulation)
    find_common_and_ranking(gu, wb, model='tf-regulation')
    # locate ranking -- TF-both
    wb = read_as_df(outPath+wbTFboth)
    gu = read_as_df(outPath+'result/'+guTFboth)
    find_common_and_ranking(gu, wb, model='tf-both')
