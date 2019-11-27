"""

Run snpEff and SnpSift to obtain missense_variant
chr1, #nsSNP=1681, #nsSNP(from imputed Geuvadis)=23752
chr2, #nsSNP=1177, #nsSNP(from imputed Geuvadis)=16557
chr3, #nsSNP=894, #nsSNP(from imputed Geuvadis)=11983
chr4, #nsSNP=687, #nsSNP(from imputed Geuvadis)=8192
chr5, #nsSNP=836, #nsSNP(from imputed Geuvadis)=10213
chr6, #nsSNP=941, #nsSNP(from imputed Geuvadis)=12323
chr7, #nsSNP=771, #nsSNP(from imputed Geuvadis)=10992
chr8, #nsSNP=539, #nsSNP(from imputed Geuvadis)=8134
chr9, #nsSNP=662, #nsSNP(from imputed Geuvadis)=9472
chr10, #nsSNP=751, #nsSNP(from imputed Geuvadis)=8886
chr11, #nsSNP=1271, #nsSNP(from imputed Geuvadis)=15635
chr12, #nsSNP=925, #nsSNP(from imputed Geuvadis)=11535
chr13, #nsSNP=302, #nsSNP(from imputed Geuvadis)=3630
chr14, #nsSNP=643, #nsSNP(from imputed Geuvadis)=8150
chr15, #nsSNP=550, #nsSNP(from imputed Geuvadis)=7586
chr16, #nsSNP=710,  #nsSNP(from imputed Geuvadis)=11644
chr17, #nsSNP=41, #nsSNP(from imputed Geuvadis)=13616
chr18, #nsSNP=338, #nsSNP(from imputed Geuvadis)=3606
chr19, #nsSNP=1363, #nsSNP(from imputed Geuvadis)=18495
chr20, #nsSNP=455, #nsSNP(from imputed Geuvadis)=5926
chr21, #nsSNP=207, #nsSNP(from imputed Geuvadis)=2905
chr22, #nsSNP=454, #nsSNP(from imputed Geuvadis)=5672
"""
import os
import glob
import pandas as pd
import subprocess as sp

SNPEFF = '/data2/GTEx_v6/code/snpEff/snpEff.jar'
SNPSIFT =  '/data2/GTEx_v6/code/snpEff/SnpSift.jar'

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
    # set up path
    snpFolder = '/data2/Geuvadis/genotype/imputed_hrc/' #'/data2/Geuvadis/genotype/'
    snpAnnotFolder = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/useImputedGenotype/vcfAnnotateRSID/' #'/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/annotation/tf_snp_annot/both/'
    outFolder = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/useImputedGenotype/vcfSNPeffect/'#'/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/annotation/snp_annot/exonic/'
    # retrieve genotype snp file
    snpList = retrieve_files(snpFolder, 'dose.vcf.gz') #gz
    # annotating snp file
    for snp in snpList:
        chrStr = os.path.basename(snp).split('.')[0].split('chr')[1]
        cmd = 'java -jar ' + SNPEFF + ' -onlyProtein  GRCh37.75 -t ' + snp + ' | java -jar ' + SNPSIFT + ' filter "ANN[0].EFFECT has \'missense_variant\'" - | grep -v "##" | cut -f1,2 > ' + outFolder + 'annot'+ '.chr' + chrStr +'.vcf'
        process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
        process.wait()
        # match
        missenseSNP = read_as_df(outFolder + 'annot'+ '.chr' + chrStr +'.vcf')
        # obtain snp annotation file
        snpAnnot = read_as_df(snpAnnotFolder + 'geuvadis.chr' + chrStr +'.snp_annot.txt')
        # subsetting to include missenseSNP only
        foundSNP = sorted(list(set(missenseSNP['POS'].values)&set(snpAnnot['pos'].values)))
        snp = snpAnnot.set_index('pos')
        sub_annot = snp.loc[foundSNP]
        print( '#missense={:}, #$np={:}, #found={:}'.format(len(missenseSNP['POS'].values), len(snpAnnot['pos'].values), len(foundSNP)) )
        print(sub_annot.shape)
        # save to file
        sub_annot.to_csv(outFolder + 'geuvadis.chr' + chrStr +'.snp_annot.txt', header=True, index=True, sep="\t")
        print( 'chr{:}, #nsSNP={:}'.format(chrStr, sub_annot.shape[0]) )
