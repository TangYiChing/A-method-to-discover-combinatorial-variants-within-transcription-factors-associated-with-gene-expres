#!/data/anaconda3/bin/python3.6

import os
import sys
import glob
import argparse
import pandas as pd
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(filename)s %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S')
logger = logging.getLogger(__name__)

# output: tfName_modelName_genotype.txt



# algorithm
"""
get_filelist to collect all the genotype files in the given folder,
get_snp_file to collect genotype file by chromosome,
get_snp_pos to collect snp_pos by chromosome,
get_tf_geno to extract genotype by snp_pos for each tf

"""
def get_filelist(geno_folder): 
    """
    :param: geno_folder: directory path to genotype files from chr1 to chr22
    :return: a list of file with file extention ".txt"
    """
    file_list =[]
    for file in glob.glob(os.path.join(geno_folder,'*.txt')): #*.genotype*.txt or '*.sharedsamples.txt'
        file_list.append(file)

    return file_list

def get_snp_file(file_list, snp_file_dict):
    """
    :param: file_list. a list of file
    :param: snp_file_dict. dictionary of snp file having the string of '.chrXX.' in the file name
    :return: snp_file_dict. {'chr1': 'Illumina.chr1.genotype.txt', 'chr2': 'Illumina.chr2.genotype.txt'}
    """
    for chr in snp_file_dict.keys():
        key = '.' + chr + '.'  #e.g.,  .chr1. ('.chrXX.' should be in filename)
        filepath = [f for f in file_list if key in f][0] # find corresponding file with the string of '.chrXX.'
        snp_file_dict[chr] = filepath #snp_file_dict['chr1']:'../../illumina/Illumina.chr1.txt' 
    return snp_file_dict

def get_snp_pos(snp_file_dict, snp_pos_dict):
    """
    :param: snp_file_dict. a dictionary of snp files. e.g., {'chr1': 'Illumina.chr1.genotype.txt', 'chr2': 'Illumina.chr2.genotype.txt'}
    :param: snp_pos_dict. a dictionary of snp position (varID). e.g., {'chr1': ['1_111_G_A_bc37', ...., '1_222_A_T_bc37']}
    """
    for chr in snp_pos_dict.keys():
        snp_file = snp_file_dict[chr]
        df = pd.read_table(snp_file, sep="\t", header=0, usecols=['Id']) #from genotype file to extract snp position(varID) at the Id column.
        snp_pos_dict[chr] = df['Id'].values.tolist() #snp_pos_dict['chr1'] = ['1_111_G_A_bc37', ...., '1_222_A_T_bc37']         
    return snp_pos_dict

def get_tf_geno(tf_table, snp_pos_dict, snp_file_dict, window_size, outPath, suffix):
    """
    :param: tf_table: txt file of gene-tf-annot.
    :param: snp_pos_dict: dict. a dict of snp position(varID) from chr1 to chr22.
    :param: snp_file_dict: dict. a dict of genotype file from chr1 to chr22.
    :param: window_size: int. searching region. 0 for tf-binding/tf-regulation, while 1 for tf-both model.
    :param: outPath: str. output folder.
    :param: suffix: str. suffix of output file
    :return: txt file of tf-genotype-index
    """
    df = pd.read_table(tf_table, sep="\t", header=0, usecols=['TF_chr', 'TF_id', 'TF_name', 'TF_start', 'TF_end'])
    unique_tfs = list(df['TF_name'].unique())
    tf_snp_dict = {tf:None for tf in unique_tfs}
    logger.debug("Number of tfs={:}".format(len(unique_tfs)))
    for tf in unique_tfs:
        tf_snp = []
        tf_df = df.loc[df['TF_name'] == tf].head(1)
        start = int(tf_df['TF_start'].values[0]) - window_size
        end = int(tf_df['TF_end'].values[0]) + window_size
        chrom = 'chr' + str(tf_df['TF_chr'].values[0])
    
        if chrom not in ['chrX', 'chrY']:  #collect snps
            for snp_pos in snp_pos_dict[chrom]: #snp_pos = '1_1000_C_G_b37'
                pos = int(snp_pos.split('_')[1])
                if pos >= start and pos <= end:
                    tf_snp.append(snp_pos)
        tf_snp_dict[tf] = tf_snp
        logger.debug("TF={:}, number of snps={:}".format(tf, len(tf_snp)))
        if chrom not in ['chrX', 'chrY']:  #extract genotype
            geno_file = snp_file_dict[chrom]
            geno_df = pd.read_table(geno_file, sep="\t", header=0, index_col=0)
            tf_geno_df = geno_df.loc[tf_snp_dict[tf], :]
            fname = outPath + tf + '.' + chrom + '.' + suffix + '.txt'
            tf_geno_df.to_csv(fname, sep="\t", header=True, index=True)
            logger.debug("TF={:}, chr={:}, geno_file={:}".format(tf, chrom, geno_file))
            logger.debug("tf_geno_file={:}, {:}".format(fname,tf_geno_df.shape))
def parse_parameter():
    """
    """
    parser = argparse.ArgumentParser(description = "Extracting TF-genotype.")
    parser.add_argument("--gene_tf_annot",
                        required = True,
                        default = "/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/gtex6_pipeline/files/annotation/gene-tf-pair.big_table.GRC37.positions.sorted.chr1-22.txt",
                        help = "input: gene-tf-big-table.txt. Default: gene-tf-pair.big_table.GRC37.positions.sorted.chr1-22.txt")
    parser.add_argument("--genotype_dir",
                        required = True,
                        help = "input: a folder of genotype file")
    parser.add_argument("--window_size",
                        type = int,
                        required = True,
                        help = "input: 1000000 for TF-both model, 0 for TF-regulation or TF-binding model")
    parser.add_argument("--out_folder",
                        required = True,
                        help = "Output folder.")

    parser.add_argument("--suffix",
                        required = True,
                        help = "suffix of output file.")
    return parser.parse_args()
# algorithm
if __name__ == "__main__":
    # initialize
    args = parse_parameter()
    tf_table = args.gene_tf_annot #sys.argv[1]  #input: gene-tf-big-table.txt
    geno_folder = args.genotype_dir #sys.argv[2] #input: a folder of genotype file
    window_size = args.window_size #int(sys.argv[3]) #input: 1000000 for TF-both model, 0 for TF-regulation or TF-binding model
    outPath = args.out_folder #sys.argv[4]
    suffix = args.suffix

    chr_list = list(range(1, 23)) # chr_list = [1,2,...,22]
    snp_file_dict = {'chr'+str(chr):None for chr in chr_list} #snp_file_dict = {'chr1': None}
    snp_pos_dict = {'chr'+str(chr):None for chr in chr_list}  #snp_pos_dict = {'chr1': None}
    logger.debug("tf table={:}".format(tf_table))
    logger.debug("genotype folder={:}".format(geno_folder))
    logger.debug("window size={:}".format(window_size))


    #indexing genotype files using chrXX as key
    geno_file_list = get_filelist(geno_folder)
    snp_file_dict = get_snp_file(geno_file_list, snp_file_dict)
    logger.debug("Total files in folder={:}".format(len(geno_file_list)))
    for chrom in snp_file_dict.keys():
        if snp_file_dict[chrom] == None:
            logger.debug("ERROR, EMPTY!!!")
            print(chrom, snp_file_dict[chrom])
            sys.exit(1)

    #retrieve snp position from genotype file
    snp_pos_dict = get_snp_pos(snp_file_dict, snp_pos_dict)


    #extract genotype for each tf
    get_tf_geno(tf_table, snp_pos_dict, snp_file_dict, window_size, outPath, suffix)

