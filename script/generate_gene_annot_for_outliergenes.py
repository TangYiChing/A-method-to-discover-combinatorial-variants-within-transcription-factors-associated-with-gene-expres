"""
generate gene_annot.txt for outliergenes


output: [Note: same chromsome in the same file, each file needs to have at least 2 genes]
gene_pos_chr1.txt
gene_pos_chr2.txt
...
"""

import argparse
import pandas as pd
import subprocess as sp

RPATH = '/usr/bin/Rscript'
SCRIPTPATH = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'
PREDICAN = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Generate gene_annot_file for outliers.")
    parser.add_argument("--outlier_file",
                        required = True,
                        help = "File name of outlierTest. Must have Headers=[gene]")
    parser.add_argument("--genePos_file",
                        required = True,
                        help = "Gene annotation file with genename. Must have Headers=[chr, start, end, gene_id, genename]")
    parser.add_argument("--out_prefix",
                        nargs = "?",
                        default = "result",
                        help = "Prefix of output file.")
    args = parser.parse_args()

    # initialize
    outlierFile = args.outlier_file
    genePosFile = args.genePos_file
    prefix = args.out_prefix

    # load data
    outliers = read_as_df(outlierFile)
    genePos = read_as_df(genePosFile)
    outliers_genePos = genePos.loc[genePos['gene_id'].isin(outliers['gene'])]

    # save to file
    outliers_genePos.to_csv(prefix+'.outliers.gene_annot.txt', header=True, index=False, sep="\t")
        
    # convert gene_annot txt to gene_annot RDS
    cmd = RPATH + ' ' + PREDICAN + 'geno_annot_to_RDS.R ' + prefix+'.outliers.gene_annot.txt ' + prefix+'.outliers.gene_annot.RDS' 
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    print(cmd)
    process.wait()

