"""
Geuvadis eQTL data: https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/
Gene eQTL https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/EUR373.gene.cis.FDR5.all.rs137.txt
 

Find how many eQTLs in the snp annotation file, chromosome by chromosome,
which will be used later on for TF-regulation model.
"""
import os
import pandas as pd
import glob 

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def retrieve_files(folder, fextention):
    """
    :param folder: string representing path to the folder where files reside
    :param fextention: string representing file extention. For example txt, RDS
    :return fileList: a list containing files with given extention
    """
    fileList = glob.glob(folder+'*.'+fextention)
    return fileList


if __name__ == "__main__":
    # set path
    outFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/annotation/tf_snp_annot/eQTL/'
    eQTLFile = '/repo4/ytang4/GTEx_v6/geuvadis_data/EUR373.gene.cis.FDR5.all.rs137.txt'
    snpFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/annotation/tf_snp_annot/both/'
    # get snp annotation files and eQTL file
    snpList = retrieve_files(snpFolder, 'txt')
    eQTL = read_as_df(eQTLFile)
    # find overlapping rsids
    for snp in snpList:
        # load snp file
        snpDf = read_as_df(snp)
        # match by chromosome first 
        chrStr = os.path.basename(snp).split('.')[1].split('chr')[1]
        chr_eQTL =  eQTL.loc[eQTL['CHR_SNP']==int(chrStr)].index.tolist() #SNP_ID
        # subsetting snp file to include matched rsids only
        snp_eQTL = snpDf.loc[snpDf['rsid'].isin(chr_eQTL)]
        # save to file
        print( 'chr{:}, #eQTL={:}, #found in snp annotation file={:}({:.2f}%)'.format(chrStr, len(chr_eQTL), snp_eQTL.shape[0], snp_eQTL.shape[0]/len(chr_eQTL)*100 ) )
        snp_eQTL.to_csv(outFolder + os.path.basename(snp), header=True, index=True, sep="\t")
