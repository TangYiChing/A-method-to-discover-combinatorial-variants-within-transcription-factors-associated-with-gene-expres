"""
snp folder (unprocessed): /data2/Geuvadis/genotype/imputed_hrc/

snp folder (processed): /data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/useImputedGenotype/vcfGenotype
snp annot folder (processed): /data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/useImputedGenotype/vcfAnnotateRSID
"""
import os
import sys
import glob
import pandas as pd
import subprocess as sp

#RPATH = '~/anaconda3/envs/bmi6319/bin/Rscript'
#SCRIPTPATH = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/script/'
#PREDICAN = '/data2/GTEx_v6/code/PredictDBPipeline-master/scripts/'
#SNPEFF = '/data2/GTEx_v6/code/snpEff/snpEff.jar'
#SNPSIFT =  '/data2/GTEx_v6/code/snpEff/SnpSift.jar'
VCFTOOLS = '/usr/local/bin/vcftools'

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t", compression='gzip')
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
    # set args
    debug = True
    outPath = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/genotype/'
    imputedFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/imputed_hrc/'
    genotypeFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/vcfGenotype/'
    snpAnnotFolder = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/vcfAnnotateRSID/'
    prefix = 'geuvadis'
    chrom = sys.argv[1]
    print(chrom)
    # load data -- info.gz
    infoFiles = retrieve_files(imputedFolder, 'info.gz')
    print(infoFiles)
    qualityDf = {}
    for infoFile in infoFiles:
        # get file chr by chr
        fStr = os.path.basename(infoFile) 
        chrStr = fStr.split('.')[0]
        if chrStr == chrom:
            print( 'processing chr={:}, {:}'.format(chrStr, infoFile) )
            info_df = read_as_df(infoFile)
            # subsetting to include quality SNPs only
            qualified = info_df[(info_df['MAF'] >= 0.01) & (info_df['Rsq'] >= 0.8)]
            # add varID column and pos column
            qualified.reset_index(inplace=True)
            varID = qualified['SNP'].str.split(':', expand=True)
            qualified['varID'] = varID[0] + '_' + varID[1] + '_' + varID[2] + '_' + varID[3] + '_b37'
            qualified['pos'] = varID[1]
            qualified['chr'] = varID[0]
            # save pos and varID
            qualityDf.update( {chrStr: qualified[['SNP', 'varID', 'pos']]} )
            qualifiedPosFile = prefix+'.'+chrStr+'.pos.txt'
            qualified[['chr', 'pos']].to_csv(qualifiedPosFile, header=True, index=False, sep="\t")
            print( 'chr={:}, size={:}, {:.2f}% kept (after quality filtering)'.format(chrStr, info_df.shape[0], qualified.shape[0]/info_df.shape[0]*100, ) )
   
            # filter SNP by positions
            qualifiedVcf = prefix+'.'+chrStr
            cmd = VCFTOOLS + ' --gzvcf ' + imputedFolder + chrStr + '.dose.vcf.gz --recode --recode-INFO-all --remove-indels --out '+ qualifiedVcf + ' --positions ' + qualifiedPosFile 
            print(cmd)
            process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
            process.wait()

            # convert to genotype
            inFile = qualifiedVcf + '.recode.vcf'
            cmd = VCFTOOLS + ' --vcf ' + inFile + ' --012 --out ' + qualifiedVcf
            print(cmd)
            process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
            process.wait()
        
            # add columns and index
            indv = pd.read_csv(qualifiedVcf + '.012.indv', header=None, sep="\t")
            pos = pd.read_csv(qualifiedVcf + '.012.pos', header=None, sep="\t")
            gt = pd.read_csv(qualifiedVcf + '.012', header=None, index_col=0, sep="\t")
            gt.index = indv[0].values.tolist()
            gt.columns = pos[1].values.tolist()
            #print('gt={:}'.format(gt.shape))

            # remove redundant column
            gt_unique = gt.loc[:,~gt.columns.duplicated()]
            print('chr={:} #snp={:} (remove duplicates)'.format(chrStr, gt_unique.shape[1]))
        
            # mapping pos to varID
            posList = [str(pos) for pos in gt_unique.columns]
            varIDList = []
            for pos in posList:
                varid = qualified.loc[qualified['pos'] == pos]['varID'].values[0]
                #print(pos, varid)
                varIDList.append(varid)
            print('found varID={:}, pos={:}'.format(len(varIDList), len(posList)))
            if len(posList) == len(varIDList):
                gt_unique.columns = varIDList #pos[1].values.tolist()
            else:
                print('ERROR: varIDList does not match with posList')
            # save to file
            snp = gt_unique.T
            snp.index.name = 'Id' # varID*samples
            snp.to_csv(outPath + qualifiedVcf + '.snp.txt', header=True, index=True, sep="\t") 
            print('chr={:}, #snp={:}'.format(chrStr, snp.shape[0]))
