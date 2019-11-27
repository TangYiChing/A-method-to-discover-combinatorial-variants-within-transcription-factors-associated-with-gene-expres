"""
Use eQTLs of GTEX whole blood to train TF-regulation model for Geuvadis dataset

chr=chr1, #snp=249341, #found eQTLs=14674, 14674
chr=chr2, #snp=276789, #found eQTLs=12836, 12836
chr=chr3, #snp=203772, #found eQTLs=7983, 7983
chr=chr4, #snp=194171, #found eQTLs=6195, 6195
chr=chr5, #snp=210028, #found eQTLs=7696, 7696
chr=chr6, #snp=213456, #found eQTLs=9711, 9711
chr=chr7, #snp=179211, #found eQTLs=8920, 8920
chr=chr8, #snp=183170, #found eQTLs=5801, 5801
chr=chr9, #snp=152296, #found eQTLs=5633, 5633
chr=chr10, #snp=174532, #found eQTLs=7017, 7017
chr=chr11, #snp=165935, #found eQTLs=7801, 7801
chr=chr12, #snp=158480, #found eQTLs=7798, 7798
chr=chr13, #snp=130020, #found eQTLs=2603, 2603
chr=chr14, #snp=105063, #found eQTLs=4043, 4043
chr=chr15, #snp=90449, #found eQTLs=5149, 5149
chr=chr16, #snp=92592, #found eQTLs=5784, 5784
chr=chr17, #snp=74513, #found eQTLs=7394, 7394
chr=chr18, #snp=98603, #found eQTLs=2101, 2101
chr=chr19, #snp=49570, #found eQTLs=6453, 6453
chr=chr20, #snp=80189, #found eQTLs=4109, 4109
chr=chr21, #snp=42861, #found eQTLs=2186, 2186
chr=chr22, #snp=42456, #found eQTLs=4200, 4200
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
    outFolder = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/annotation/tf_snp_annot/eQTL/gtex_wb/'
    GTEXeQTLFolder = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/gtex_model_TF/tissues/whole_blood/annotation/tf_snp_annot/eQTLs/'
    snpFolder = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/annotation/tf_snp_annot/both/'

    # retrieve files
    chrList = [ 'chr'+str(i) for i in range(1,23) ]
    snpFiles = retrieve_files(snpFolder, 'txt')
    eQTLFiles = retrieve_files(GTEXeQTLFolder, 'txt')
    snpDict = {}
    eQTLDict = {}
    for chrStr in chrList:
        chrNumStr = '.'+chrStr+'.'
        # add to snpDict
        for snpFile in snpFiles:
            if chrNumStr in snpFile:
                snpDict.update({chrStr:snpFile})
        # add to eQTLDict
        for eQTLFile in eQTLFiles:
            if chrNumStr in eQTLFile:
                eQTLDict.update({chrStr:eQTLFile})    
    
    # find intersection
    for chrStr in chrList:
        snps = read_as_df(snpDict[chrStr])
        eQTLs = read_as_df(eQTLDict[chrStr])
        sharedSNPs = sorted(list(set(snps['varID']) & set(eQTLs['varID'])))
        snp_eQTL = snps.loc[snps['varID'].isin(sharedSNPs)]
        print( 'chr={:}, #snp={:}, #found eQTLs={:}/{:}'.format(chrStr, snps.shape[0], len(sharedSNPs), eQTLs.shape[0]) )

        #save to file
        #snp_eQTL.to_csv(outFolder+'geuvadis.'+chrStr+'.snp_annot.txt', header=True, index=True, sep="\t")
    
