"""
1. samples
2. expressed genes
3. number of TFs for each TF-models
"""
import os
import glob
import numpy as np
import pandas as pd

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def get_sample_population(sampleFile):
    """
    :param sampleFile: str representing file name. expected columns=['Source Name','Characteristics[population]']
    :return sub_df: df containing sample and population columns
    """
    df = pd.read_csv(sampleFile, header=0, usecols=['Source Name','Characteristics[population]'], sep="\t")
    sub_df = df.drop_duplicates(['Source Name','Characteristics[population]'], keep='first')
    sub_df.columns = ['sample', 'population']
    return sub_df # headers=['sample', 'population']

if __name__ == "__main__":
    # set args
    debug = True
    outPath = './'
    # input files
    dataPATH = '/repo4/ytang4/GTEx_v6/geuvadis_data/useImputedGenotype/'
    expFile = dataPATH + 'expression/geuvadis.expression.txt'
    snpFile = dataPATH + '/genotype/geuvadis.chr22.snp.txt'
    snpAnnotFolder = dataPATH + '/annotation/tf_snp_annot/both/'
    genetfFile = dataPATH + '/annotation/tf_gene_annot/geuvadis.tfgene_annot.txt'
    nsSNPFolder = dataPATH + '/genotype/tf-binding_genotype_index/'
    eQTLFolder = dataPATH + '/genotype/tf-regulation_genotype_index/'
    bothFolder = dataPATH + '/genotype/tf-both_genotype_index/'
    geneAnnotFile = '/repo4/ytang4/GTEx_v6/geuvadis_data/gencode.v12.annotation.gtf'
    sampleInfoFile = '/repo4/ytang4/GTEx_v6/geuvadis_data/E-GEUV-3.sdrf.txt'
    # load data
    exp = read_as_df(expFile) # gene*sample
    snp = read_as_df(snpFile) # varID*sample
    tfgene_annot = read_as_df(genetfFile) # headers=['gene_id', 'TF_id']
    sample_population = get_sample_population(sampleInfoFile) # headers=['sample', 'population']    

    # summary -- samples
    sharedSamples = sorted(list(set(exp.columns) & set(snp.columns)))
    sample_dict = {}
    for population in sample_population['population'].unique():
        sample_dict.update({population:[]})
    for sample in sharedSamples:
        population = sample_population.loc[sample_population['sample']==sample]['population'].values[0]
        sample_dict[population].append(sample)
    for key, value in sample_dict.items():
        print( 'population={:}, #sample={:}'.format(key, len(value)) )
    print( '#samples {:} in exp, {:} in genotype, {:} shared'.format(len(exp.columns), len(snp.columns), len(sharedSamples)) )
    
    # summary -- genes
    foundGenes = sorted(list(set(exp.index)&set(tfgene_annot['gene_id'].unique())))
    TFNames = []
    n_TF = []
    for gene in foundGenes:
        tfs = tfgene_annot.loc[tfgene_annot['gene_id']==gene]['TF_id'].values.tolist()
        tfnames = tfgene_annot.loc[tfgene_annot['gene_id']==gene]['TF_name'].values.tolist()
        TFNames+=tfnames
        n_TF.append(len(tfs))
    print( '# expressed genes={:}, {:} found in tfgene_annot'.format(len(exp.index), len(foundGenes)) )
    print( '# TF per gene: mean={:}, std={:}'.format(np.mean(n_TF), np.std(n_TF)) )

    # summary -- SNPs
    TF_nsSNP = []
    TF_eQTL = []
    TF_both = []
    uniqueTFs = list(set(TFNames))
    for tf in uniqueTFs:
        tf_start = tfgene_annot[tfgene_annot['TF_name']==tf]['TF_start'].values[0]
        tf_end = tfgene_annot[tfgene_annot['TF_name']==tf]['TF_end'].values[0]
        tf_chr = tfgene_annot[tfgene_annot['TF_name']==tf]['TF_chr'].values[0]
        #print( 'tf={:}, chr={:}, start={:}, end={:}'.format(tf, tf_chr, tf_start, tf_end) )
        #snp_annot
        nsSNPAnnotFile = dataPATH + 'annotation/tf_snp_annot/nsSNP/'  + 'geuvadis.chr' + str(tf_chr) + '.snp_annot.txt'
        eQTLAnnotFile = dataPATH + 'annotation/tf_snp_annot/eQTL/'  + 'geuvadis.chr' + str(tf_chr) + '.snp_annot.txt'
        bothAnnotFile = dataPATH + 'annotation/tf_snp_annot/both/'  + 'geuvadis.chr' + str(tf_chr) + '.snp_annot.txt'
        nsSNP_snp_annot = read_as_df(nsSNPAnnotFile)
        eQTL_snp_annot = read_as_df(eQTLAnnotFile)
        both_snp_annot = read_as_df(bothAnnotFile)
        #print( 'nsSNP file={:}, eQTL file={:}, both file={:}'.format(nsSNPAnnotFile, eQTLAnnotFile, bothAnnotFile) ) 
        # find overlaps
        nsSNPs = nsSNP_snp_annot['pos'].values.tolist()    
        eQTLs = eQTL_snp_annot['pos'].values.tolist()
        boths = both_snp_annot['pos'].values.tolist()
        nsSNPs_pos = [] # a list of nsSNP pos    
        eQTLs_pos = []  # a list of eQTL pos
        boths_pos = []  # a list of boths pos
        for pos in nsSNPs:
            if int(pos) >= int(tf_start) and int(pos) <= int(tf_end):                
                nsSNPs_pos.append(pos)  
        for pos in eQTLs:
            if int(pos) >= int(tf_start) and int(pos) <= int(tf_end):
                eQTLs_pos.append(pos)
        for pos in boths:
            if int(pos) >= int(tf_start)-1000000 and int(pos) <= int(tf_end)+1000000:
                boths_pos.append(pos)
        # update number of snps for each TF
        TF_nsSNP.append(len(nsSNPs_pos))  
        TF_eQTL.append(len(eQTLs_pos))
        TF_both.append(len(boths_pos))
        #print( 'tf={:}, #nsSNP={:}, #eQTL={:}, #both={:}'.format(tf, len(nsSNPs), len(eQTLs), len(boths)) ) 
    print( '# nsSNP per TF: mean={:}, std={:}'.format(np.mean(TF_nsSNP), np.std(TF_nsSNP)) )
    print( '# eQTL per TF: mean={:}, std={:}'.format(np.mean(TF_eQTL), np.std(TF_eQTL)) )
    print( '# SNPs within flanking region per TF: mean={:}, std={:}'.format(np.mean(TF_both), np.std(TF_both)) )
