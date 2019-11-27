import sys
import pandas as pd


chrStr = sys.argv[1] # e.g., chr1
nsSNPFile = './vcfSNPeffect/geuvadis.' + chrStr + '.snp_annot.txt'
snpAnnotFile = './annotation/tf_snp_annot/both/geuvadis.' + chrStr + '.snp_annot.txt' 
outFolder = './annotation/tf_snp_annot/nsSNP/'

print( '{:}, snp_annot={:}, nsSNP file={:}'.format(chrStr, snpAnnotFile, nsSNPFile) )

nsSNP = pd.read_csv(nsSNPFile, header=0, index_col=0, sep="\t")
snp_annot = pd.read_csv(snpAnnotFile, header=0, index_col=0, sep="\t")
snps = sorted( list(set(nsSNP['varID']) & set(snp_annot['varID'])) )
qualified_nsSNP = nsSNP.loc[nsSNP['varID'].isin(snps)]
qualified_nsSNP.to_csv(outFolder+'geuvadis.'+chrStr+'.snp_annot.txt', header=True, index=True, sep="\t")
print( 'found={:}, {:}, #nsSNP={:}'.format(len(snps), chrStr, qualified_nsSNP.shape[0]) )
