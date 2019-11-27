"""
1. expression: gene*sample, value=PKM, PEER adjusted
"""
import sys
import pandas as pd
import subprocess as sp
from gtfparse import read_gtf

RPATH = '/usr/bin/Rscript'
SCRIPTPATH = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'
PREDICAN = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def process_expression(fpath, outPath):
    """
    :param fpath: string representing path to file
    :param outPath: string representing path for output file
    :return output: df in gene*sample format
    """
    # load data
    expDf = read_as_df(fpath)
   
    # remove columns: GeneSymbol. Chr, Coord
    gene_sample = expDf.iloc[:, 3:]
    
    # rename index
    gene_sample.index.name = 'gene_id'
    
    return gene_sample

def process_snp_annot(fpath, outPath):
    """
    :param fpath: string representing path to file
    :param outPath: string representing path for output file
    :return output: df containing headers=[chr,pos,varID,refAllele,effectAllele,rsid]
    """
    # load data
    snpDf = read_as_df(fpath)

    # select wanted columns
    cols = ['Pos', 'VariantID', 'Ref_b37', 'Alt', 'RSID_dbSNP137']
    sub_snp = snpDf[cols]
  
    # rename columns
    sub_snp.columns = ['pos','varID','refAllele','effectAllele','rsid']
    sub_snp.index.name = 'chr'
  
    # save to files chr by chr
    #for i in range(0, 21):
    #    n_chr = i + 1
    #    chrDf = sub_snp.loc[n_chr, :]
    #    print( 'chr={:}, size={:}'.format(n_chr, chrDf.shape) )
    #    chrDf.to_csv(outPath+'geuvadis.chr'+str(n_chr)+'.snp_annot.txt')

    return sub_snp

def process_gene_annot(fpath, outPath):
    """
    :param fpath: string representing path to file
    :param outPath: string representing path for output file
    :return output: df containing headers=[chr, gene_id, genename, start, end]

    Note:
    -----
    use gtfparse to load gtf file. https://github.com/openvax/gtfparse
    """
    # load data
    geneDf = read_gtf(fpath)

    # retrieve genes only
    df_genes = geneDf[geneDf["feature"] == "gene"]

    # select wanted columns
    cols = ['seqname', 'gene_id',  'transcript_name', 'start', 'end']
    subdf_genes = df_genes[cols]
  
    # retrieve chr str
    chrStr = subdf_genes['seqname'].str.split('chr', n=1, expand=True)
    subdf_genes['chr'] = chrStr[1]

    # drop seqname and keep chr column
    sub_gene = subdf_genes.drop(['seqname'], axis=1)

    # reorder columns
    sub_gene = sub_gene[['chr', 'gene_id',  'transcript_name', 'start', 'end']]
    sub_gene.columns = ['chr', 'gene_id',  'genename', 'start', 'end']

    return sub_gene

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
    MainPath = '/repo4/ytang4/GTEx_v6/geuvadis_data/'
    outPath = MainPath 
    chrom = sys.argv[1] #e.g., chr1
    # input files
    expFile = MainPath + 'GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'
    snpFile = MainPath + 'useImputedGenotype/genotype/geuvadis.' + chrom + '.snp.txt'
    #snpAnnotFile = MainPath + 'useImputedGenotype/annotation/tf_snp_annot/both/geuvadis' + chrom + '.snp_annot.txt'
    geneAnnotFile = MainPath + 'gencode.v12.annotation.gtf'
    sampleInfoFile = MainPath + 'E-GEUV-3.sdrf.txt'
    #genetfFile = '/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/gtex6_pipeline/files/annotation/gene-tf-pair.big_table.GRC37.positions.sorted.chr1-22.txt'
    # process data -- gene annotation, snp annotation, expression
    geneAnnot = process_gene_annot(geneAnnotFile, outPath) # headers=['chr', 'gene_id',  'genename', 'start', 'end']
    #snpAnnot = process_snp_annot(snpAnnotFile, outPath) # headers=['pos','varID','refAllele','effectAllele','rsid']
    exp = process_expression(expFile, outPath)   #gene*sample
    sample_population = get_sample_population(sampleInfoFile) # headers=['sample', 'population']
    # find shared genes in expression and gene annotation, shared samples in expression and snp (i.e. genotype)
    snp = pd.read_csv(snpFile, header=0, index_col=0, sep="\t") # trun off nrows to get whole data
    sharedGenes = sorted(list(set(geneAnnot['gene_id']) & set(exp.index)))
    sharedSamples = sorted(list(set(exp.columns) & set(snp.columns)))
    print( '{:} genes in exp, {:} genes in geneAnnot, shared genes={:}'.format(len(exp.index), len(geneAnnot['gene_id']), len(sharedGenes)) )
    print( '{:} samples in exp, {:} samples in snp, shared samples={:}'.format(len(exp.columns), len(snp.columns), len(sharedSamples)) )
    # count how many population 
    sample_dict = {}
    for population in sample_population['population'].unique():
        sample_dict.update({population:[]})
    for sample in sharedSamples:
        population = sample_population.loc[sample_population['sample']==sample]['population'].values[0]
        sample_dict[population].append(sample)
    for key, value in sample_dict.items():
        print( 'population={:}, #sample={:}'.format(key, len(value)) )
    # excluding samples from YRI population
    samples_woYRI = sorted( list(set(sharedSamples) - set(sample_dict['YRI'])))
    # subsetting to include shared genes, shared samples
    sub_geneAnnot = geneAnnot.loc[geneAnnot['gene_id'].isin(sharedGenes)]
    sub_exp = exp[samples_woYRI]#[sharedSamples]
    sub_snp = snp[samples_woYRI]#[sharedSamples]
    # save file
    sub_snp.to_csv(outPath+'useImputedGenotype/genotype/geuvadis.'+chrom+'.snp.txt', header=True, index=True, sep="\t")
    sub_exp.to_csv(outPath+'useImputedGenotype/expression/geuvadis.expression.txt', header=True, index=True, sep="\t")
    sub_geneAnnot.to_csv(outPath+'useImputedGenotype/annotation/gene_annot/geuvadis.gene_annot.txt', header=True, index=False, sep="\t")
    print( 'SNP: {:}, #snps={:}, #samples={:}'.format(chrom, sub_snp.shape[0], sub_snp.shape[1]) )
    print( 'EXP: {:}, #genes={:}, #samples={:}'.format(chrom, sub_exp.shape[0], sub_exp.shape[1]) )
    # convert snp_annot txt to snp_annot RDS
    #cmd = RPATH + ' ' + PREDICAN + 'single_snp_annot_to_RDS.R ' + outPath+'geuvadis.chr'+str(n_chr)+'.snp_annot'
    #process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    #print(cmd)
    #process.wait()


    # covert exp txt to exp RDS
    cmd = RPATH + ' ' + PREDICAN + 'expr_to_transposed_RDS.R ' + outPath+'useImputedGenotype/expression/geuvadis.expression.txt' + ' ' + outPath+'useImputedGenotype/expression/geuvadis.expression.RDS'
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    #print(cmd)
    process.wait()
    # convert gene_annot txt to gene_annot RDS
    cmd = RPATH + ' ' + PREDICAN + 'geno_annot_to_RDS.R ' + outPath+'useImputedGenotype/annotation/gene_annot/geuvadis.gene_annot.txt' + ' ' + outPath+'useImputedGenotype/annotation/gene_annot/geuvadis.gene_annot.RDS'
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    #print(cmd)
    process.wait()
