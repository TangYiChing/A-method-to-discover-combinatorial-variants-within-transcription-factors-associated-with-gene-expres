import numpy as np
import pandas as pd
import glob
import os

genoFolder = './genotype/tf-binding_genotype_index/'
snpAnnotFolder = './annotation/tf_snp_annot/nsSNP/'

genoIdxFiles = glob.glob(genoFolder+'*tf-binding.txt')
n_snps = []
for genoFile in genoIdxFiles:
    chrStr = os.path.basename(genoFile).split('.')[1]
    annotFile = snpAnnotFolder + 'geuvadis.' + chrStr + '.snp_annot.txt'

    geno = pd.read_csv(genoFile, header=0, index_col=0, sep="\t")
    annot = pd.read_csv(annotFile, header=0, index_col=0, sep="\t")
    
    snps = sorted( list(set(geno.index) & set(annot['varID'].values)) )
    sub_geno = geno.loc[geno.index.isin(snps)]
    sub_geno.to_csv(genoFile, header=True, index=True, sep="\t")

    n_snps.append(len(snps))
    print( '{:}, genoFile={:}, annotFile={:}'.format(chrStr, genoFile, annotFile) )
    print( 'found {:} varID, sub_geno={:}'.format(len(snps), sub_geno.shape[0]) )
print( '#nsSNPs per TF={:}+-{:}'.format(np.mean(n_snps), np.std(n_snps)) )
