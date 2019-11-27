"""
given a list of varID, to generate snp_annot.RDS
"""


import sys
import pandas as pd
import subprocess as sp
from gtfparse import read_gtf

RPATH = '/usr/bin/Rscript'
SCRIPTPATH = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'
PREDICAN = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'

snpFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/genotype/'
snpAnnotFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/vcfAnnotateRSID/'
outFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/annotation/tf_snp_annot/both/'
# set up
chrom = sys.argv[1] # e.g. chr1
snp = pd.read_csv(snpFolder+'geuvadis.'+chrom+'.snp.txt', header=0, usecols=['Id'], sep="\t")
snpAnnot = pd.read_csv(snpAnnotFolder+'geuvadis.'+chrom+'.snp_annot.txt', header=0, index_col=0, sep="\t")
# drop duplicates if any
snp_annot = snpAnnot.drop_duplicates(subset='varID', keep='first', inplace=False)
# subsetting to include varID in snp only
varIDList = sorted( list(set(snp_annot['varID']) & set(snp['Id'])) )
df = snp_annot.loc[snp_annot['varID'].isin(varIDList)]
# save to file
df.to_csv(outFolder+'geuvadis.'+chrom+'.snp_annot.txt', header=True, index=True, sep="\t" )
print( 'chr={:}, snpAnnot={:}, #duplicates={:}'.format(chrom, snpAnnotFolder+'geuvadis.'+chrom+'.snp_annot.txt', snpAnnot.shape[0]-snp_annot.shape[0]) )
print( 'found #varID={:}, chr={:}, #snp={:}'.format(len(varIDList), chrom, df.shape[0]) )
# convert snp_annot txt to snp_annot RDS
cmd = RPATH + ' ' + PREDICAN + 'single_snp_annot_to_RDS.R ' + outFolder+'geuvadis.'+chrom+'.snp_annot'
process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
print(cmd)
process.wait()

# check number in snp
if snp.shape[0] > len(varIDList):
    snp = pd.read_csv(snpFolder+'geuvadis.'+chrom+'.snp.txt', header=0, index_col=0, sep="\t")
    df = snp.loc[varIDList, :]
    df.to_csv(snpFolder+'geuvadis.'+chrom+'.snp.txt', header=True, index=True, sep="\t" )
    print( 'found #varID={:}, chr={:}, #snp={:}'.format(len(varIDList), chrom, df.shape[0]) )
