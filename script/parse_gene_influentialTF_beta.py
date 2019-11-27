"""
Given a list of genea, generate influentialTF information containing gene, transSNP, cisWeight, beta, influential TF
"""

import os
import sys
import glob
import argparse
import pandas as pd

def read_as_df(input_path):
    """
    Function to read tab-delimited flat file (1st row is header) from imput path.
    :param input_path: file name str
    :return df: pandas dataframe
    """
    if input_path.endswith('.gz'):
        df = pd.read_csv(input_path, compression='gzip', header=0, sep='\t', quotechar='"', error_bad_lines=False)
    else:
        df = pd.read_csv(input_path, sep="\t", header=0)

    return df

def has_headers(df, headerList):
    notIn = []
    for header in headerList:
        if not header in df.columns:
            notIn.append(header)
    if len(notIn) != 0:
        print( 'ERROR: must have headers={:}, missing={:}'.format(headerList, notIn) )
        sys.exit(1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Given a list of genea, generate influentialTF information containing gene, transSNP, beta, influential TF")
    parser.add_argument("-g", "--geneListFile", 
                        required = True,
                        help = "File name representing gene. Must have headers=[gene]")
    parser.add_argument("-w", "--weightFolder",
                        required = True,
                        help = "Folder name representing weight table from modeling training")
    parser.add_argument("-tf", "--tfgeneAnnotFile",
                        help = "Gene-TF annotation file.")
    parser.add_argument("-m", "--model_name",
                        choices = ['tf-binding', 'tf-both', 'tf-regulation'],
                        required = True,
                        help = "Choices of model name: [ tf-binding | tf-both | tf-regulation ]")
    parser.add_argument("-p", "--out_prefix",
                        nargs = "?",
                        default = "result",
                        help = "Prefix of output file.")
    args = parser.parse_args()

    # initialize
    weightPath = args.weightFolder
    tfgeneFile = args.tfgeneAnnotFile
    geneFile = args.geneListFile
    model = args.model_name
    prefix = args.out_prefix

    
    # load data
    gene = read_as_df(geneFile)
    tfgene = read_as_df(tfgeneFile)
    weightFileList = glob.glob(weightPath+'*weight*.txt')

    # check headers
    has_headers(gene, ['gene'])
    has_headers(tfgene, ['gene_id', 'TF_id'])

    # retrieve all weight files
    weightDict={} # chr:filename
    for weightFile in weightFileList:
        chrStr = weightFile.split('_chr')[1].split('.txt')[0]
        weightDict[chrStr]=weightFile

    # retrieve gene list
    geneList = gene['gene'].values.tolist()

    # define searching window
    if model == 'tf-both':
         window = 1000000
    else:
        window = 0
      
    # influential-SNPs
    gene_transWeight = {}
    dfs = []
    n_found = 0
    n_gene_wo_influentialSNP = 0
    n_gene_wo_transSNP = 0
    for gene in geneList:
        geneChr = tfgene.loc[tfgene['gene_id']==gene]['gene_chr'].values[0].astype(str)
        # load weight file
        weight = read_as_df(weightDict[geneChr])
        weight_woNA = weight.dropna(axis = 0)
        has_headers(weight_woNA, ['gene',  'rsid', 'beta'])
        weight_woNA.loc[:, 'chr'] = weight_woNA['varID'].str.split('_').str[0]
        # sanity check if gene in the weight table
        if gene in weight_woNA['gene'].values.tolist():
            n_found += 1
        # subsetting to include the current gene only
        geneWeight = weight_woNA.loc[weight_woNA['gene']==gene]
        geneCisWeight = geneWeight.loc[geneWeight['chr']==geneChr]
        geneTransWeight = geneWeight.loc[geneWeight['chr']!=geneChr]
        # skip gene if it does not have transSNP
        if geneTransWeight.shape[0] != 0:
            maximalCis = geneCisWeight['beta'].abs().max()
            influentialTrans = geneTransWeight[ (geneTransWeight['beta'].abs() > maximalCis) ]
            if influentialTrans.shape[0] != 0:
                influentialTrans.loc[:, 'maximalCis'] = maximalCis
                dfs.append(influentialTrans)
                gene_transWeight.update( {gene: influentialTrans} )
            else:
                n_gene_wo_influentialSNP += 1
                #print( 'gene={:} does not have influential trans-SNP'.format(gene) )   
        else:
            n_gene_wo_transSNP += 1 
            #print( 'gene={:} does not have trans-SNP'.format(gene) )
    print( 'total genes={:}, found={:}, #genes has NO influential SNP={:}, #genes has NO trans-SNP={:}'.format(len(geneList), n_found, n_gene_wo_influentialSNP, n_gene_wo_transSNP) )

    # loop through dict to add TF for each  transSNP
    TFdfs = []
    for key, value in gene_transWeight.items():
        gene = key #gene id
        transSNP = value
        genename = tfgene.loc[tfgene['gene_id']==gene]['genename'].values[0]
        for i in range(len(transSNP)):
            df = transSNP.iloc[i]
            snp_chr = int(df['varID'].split('_')[0])
            snp_pos = int(df['varID'].split('_')[1])
            tfs = tfgene.loc[tfgene['gene_id']==gene]['TF_name'].values.tolist()
            # grasp all tfs
            if len(tfs) > 0:
                for tf in tfs:
                    tf_chr = int(tfgene.loc[tfgene['TF_name']==tf]['TF_chr'].values[0])
                    if tf_chr == snp_chr:
                             tf_start = int(tfgene.loc[tfgene['TF_name']==tf]['TF_start'].values[0])
                             tf_end = int(tfgene.loc[tfgene['TF_name']==tf]['TF_end'].values[0])
                             # searching for the right TF
                             if (snp_pos>= tf_start-window and snp_pos<=tf_end+window):
                                 #print( 'gene={:}, snpPos={:}, tf={:}, tfStart={:}, tfEnd={:}'.format(gene, snp_pos, tf, tf_start, tf_end) )
                                 df['TF'] = tf
                                 df['genename'] = genename
                                 TFdfs.append(df)



    # concatenate results 
    if len(TFdfs) != 0:
         influentialSNP = pd.concat(TFdfs, axis=1)
         influentialTFs = []
         print( 'total influential snps={:}\n{:}'.format(influentialSNP.T.shape[0], influentialSNP.T) ) 

         # add influential TF
         #genes = influentialSNP['gene'].values.tolist()
         #for gene in genes:
         #    snp_chr = int(influentialSNP.loc[influentialSNP['gene']==gene]['chr'].values[0])
         #    snp_pos = int(influentialSNP.loc[influentialSNP['gene']==gene]['varID'].str.split('_').str[1].values[0])
         #    tfs = tfgene.loc[tfgene['gene_id']==gene]['TF_name'].values.tolist()
         #    if len(tfs) != 0:
         #        for tf in tfs:
         #            tf_chr = int(tfgene.loc[tfgene['TF_name']==tf]['TF_chr'].values[0])
         #            if tf_chr == snp_chr:
         #                tf_start = int(tfgene.loc[tfgene['TF_name']==tf]['TF_start'].values[0])
         #                tf_end = int(tfgene.loc[tfgene['TF_name']==tf]['TF_end'].values[0])
         #                if (snp_pos>= tf_start-window and snp_pos<=tf_end+window):
                             #print( 'gene={:}, snpPos={:}, tf={:}, tfStart={:}, tfEnd={:}'.format(gene, snp_pos, tf, tf_start, tf_end) )
         #                    influentialTFs.append(tf)
    
         # concate TF
         #influentialSNP.loc[:, 'influential TF'] = influentialTFs
         #print( 'influentialSNPs={:}\n{:}'.format(influentialSNP.shape[0], influentialSNP) )
         # save to file
         influentialSNP.T.to_csv(prefix+'.influentialSNP.txt', header=True, index=False, sep="\t")
