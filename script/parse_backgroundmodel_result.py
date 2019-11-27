"""
"""


import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests # for bonferroni p-value

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

def has_headers(df, header_list):
    """
    Function to check whether the header_list is in file.
    :param df: pd.dataframe
    :param header_list: list of str
    :return True/False
    """
    header_in_df = df.columns.tolist()
    if header_in_df == header_list:
        pass
    else:
        print("header should have:{:} (!!!case and position sensitive!!!)".format(header_list))
        sys.exit(1)
    return

def summary_dict(result_dict):
    """
    Function to get statistics of a dict
    :param: result_dict: dictionary
    :return summary
    """
    keys = len(result_dict.keys())
    values = len(result_dict.values())
    vkeys = [key for key, values in result_dict.items() if len(values)>0]
    #print("number of keys={:}".format(keys))
    #print("number of values={:}".format(values))
    #print("number of non-zero keys={:}/{:}".format(len(vkeys), keys))
    #print("Non-zero keys={:}".format(vkeys))
    return vkeys

def collect_working_files(bgfolder_str, outlier_genes_list):
    """
    Function to collect working files and concate into a big merged file for chr1 to chr22
    """
    # collecting results from all chromsomes
    files = glob.glob(bgfolder_str+'working*.txt')
    bgresult_dict = {'chr'+str(chr):[] for chr in range(1,23)}
    for f in files:
        # retrive chr number
        nchr = os.path.basename(f).split('.txt')[0].split('chr')[-1] 
        chr_key_str = 'chr' + nchr
        bgresult_dict[chr_key_str].append(f)
    # summarize bgresult_dict
    chrList = summary_dict(bgresult_dict)
    print( 'outliers were from chromosome(s)={:}'.format(chrList) )
    bgresult_list = []
    for key, values in bgresult_dict.items():
        print( "{:} has {:} results".format(key, len(values)) )
        for item in values:
            working_df = read_as_df(item)
            has_headers(working_df, ["gene","alpha","cvm","lambda.iteration","lambda.min",
                                     "n.snps","R2","pval","genename","w_TF","pct_TF","total.snps","used_TF"])
            keep_list = list(set(working_df['gene'].values.tolist()) & set(outlier_genes_list))
            selected_df = working_df.loc[working_df['gene'].isin(keep_list)]
            n_chr = key.split('chr')[1]
            selected_df['chr'] = [n_chr] * selected_df.shape[0]
            #print(selected_df.head())
            bgresult_list.append(selected_df)
            #print("{:}, working_file={:}, gene={:}".format(key, item, keep_list))
    bgresult_df = pd.concat(bgresult_list)
    return bgresult_df

def combine_model_results(model_df, outlier_df, bg_df, gene_str_list, model_name_str, tissue_name_str):
    """
    Function to merge results from modeling and background modeling for an outlier gene
    :param: model_df: pd.dataframe, header = ["gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval","genename","w_TF","pct_TF","total.snps","used_TF"]
    :param: outlier_df: pd.dataframe, header = ["gene","PrediXcan_R2","TF_R2","genename","rstudent","unadjusted.p_value","Bonferonni.p_value"]
    :param: bg_df: pd.dataframe, header = ["gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval","genename","w_TF","pct_TF","total.snps","used_TF"]
    :param: model_name_str: str, e.g., ['tf-binding', 'tf-regulation', 'tf-both']
    :return: mg_df: pd.dataframe, header = ["TF-model", "gene", "genename", "rstudent", "unadjusted.p_value", "Bonferonni.p-value",
                                            "PrediXcan_R2", "TF_R2", "bg_avg_R2", "tf_used_TFs", "tf_total_snps", "num_of_right", "bg_p-value"]
    """

    df_list = []
    for gene_str in gene_str_list:
        gene_model_df = model_df.loc[model_df['gene']==gene_str]
        gene_outlier_df = outlier_df.loc[outlier_df['gene']==gene_str]
        gene_bg_df = bg_df.loc[bg_df['gene']==gene_str]

        if gene_bg_df.shape[0] == 0: #and gene_model_df.shape[0] == 0:
            #print(gene_str, gene_bg_df.head())
            pass
            #print("gene={:} has no background results".format(gene_str))
        else:
            #print("gene_model_df=\n{:}".format(gene_model_df.head()))
            #print("gene_outlier_df=\n{:}".format(gene_outlier_df.head()))
            #print("gene_bg_df=\n{:}".format(gene_bg_df.head()))
            # calculate num of right and evaluate its significant by adjusted p-value
            tf_total_snps = pd.to_numeric(gene_model_df['total.snps']).values[0]
            R2 = gene_outlier_df['TF_R2'].values[0]
            #print(gene_str, gene_model_df['R2'].values)
            tf_total_snps = float(tf_total_snps)
            num_of_right = 0
            for r2 in gene_bg_df['R2'].values.tolist():
               if float(r2) > R2:
                    num_of_right += 1
            bg_pvalue = float(num_of_right/100) #100 iterations of background modeling
            #print(gene_str,gene_bg_df)
            gene_records = {'TF-model': [model_name_str],
                        'gene': [gene_str],
                        'genename': [gene_model_df['genename'].values[0]],
                        'rstudent': [gene_outlier_df['rstudent'].values[0]],
                        'unadjusted.p_value': [gene_outlier_df['unadjusted.p_value'].values[0]],
                        'Bonferonni.p_value': [gene_outlier_df['Bonferonni.p_value'].values[0]],
                        'baseline_R2': [gene_outlier_df['baseline_R2'].values[0]],
                        'TF_R2': [gene_outlier_df['TF_R2'].values[0]],
                        'bg_avg_R2': np.mean([float(str(v).lstrip(' ')) for v in gene_bg_df['R2'].values]),   #[gene_bg_df['R2'].mean()]
                        'tf_total-snps': [tf_total_snps],
                        'num_of_right': [num_of_right],
                        'bg_p-value': [bg_pvalue],
                        'chr': [gene_bg_df['chr'].values[0]],
                        'tissue': [tissue_name_str]
                       }
            df = pd.DataFrame.from_dict(gene_records)
            #print(df.head())
            df_list.append(df)
    mg_df = pd.concat(df_list)
    cols = df.columns.tolist()
    #print("gene_merged_df=\n{:}".format(mg_df.head()))
    return mg_df

def p_adjust(p, method='fdr_tsbh'):
    """Bonferroni correction of p-value for multiple hypothesis testing.
    [ref]https://stackoverflow.com/questions/41517159/bonferroni-correction-of-p-values-from-hypergeometric-analysis
    """
    #method = fdr_tsbh: B&H
    #method = b: Bonferroni
    p_adjusted = multipletests(p, method=method)[1:2]

    return p_adjusted[0]

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    [ref]https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Parse background model result. Given the path of result folders, return a concatenated result file of N runs")
    parser.add_argument("--outlierFile",
                        required = True,
                        help = "result file from find_outlier.R. Must have headers = []")
    parser.add_argument("--modelResultFile",
                        required = True,
                        help = "model result file obtained from parse_model_result.py.")
    parser.add_argument("--bgfolder_path",
                        required = True,
                        help = "Location to results from background model.")
    parser.add_argument("--model_name",
                        choices = ['tf-binding', 'tf-both', 'tf-regulation'],
                        help = "Choice of TF-model: [ tf-binding | tf-both | tf-regulation ]")
    parser.add_argument("--out_prefix",
                        nargs = "?",
                        default = "result",
                        help = "Prefix of output file.")
    args = parser.parse_args()


    # initialize
    outliers_fname = args.outlierFile
    bgfolder_str = args.bgfolder_path
    model_result_str = args.modelResultFile
    tfmodel_str = args.model_name
    prefix = args.out_prefix

    # read files and checking headers
    outliers_df = read_as_df(outliers_fname)
    has_headers(outliers_df, ["gene","baseline_R2","TF_R2","genename","rstudent","unadjusted.p_value","Bonferonni.p_value"])
    outlier_genes_list = outliers_df['gene'].values.tolist()
    model_result_df = read_as_df(model_result_str)
    has_headers(model_result_df, ["gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval","genename","w_TF","pct_TF","total.snps","used_TF"])
    print( 'model={:}, #outliers={:}'.format(tfmodel_str, len(outlier_genes_list)) )

    # collecting results from all chromsomes
    bgresult_df = collect_working_files(bgfolder_str, outlier_genes_list)
    print(bgresult_df.head())

    # combing model result and background model result for each outlier
    merged_df = combine_model_results(model_result_df, outliers_df, bgresult_df, outlier_genes_list, tfmodel_str, prefix)
    merged_df.chr = pd.to_numeric(merged_df.chr, errors='coerce') # correct for sorting
    merged_df.sort_values(['chr'], inplace=True)
    # correct pvalue
    p = merged_df['bg_p-value'].values
    adjp = p_adjust(p, method='b') # method='tfdr_tsbh', will use B&H; method='b', will use Bonferroni
    merged_df['adj_p-value'] = adjp
    print( '{:} genes/outliers, {:} outlier(s) <= 0.01'.format(merged_df.shape[0], len(merged_df.loc[merged_df['adj_p-value']<0.01])) )
    # write to file
    outfname = prefix + '.' + tfmodel_str + '.bgresult.txt'
    headers = ['TF-model','gene','genename','rstudent','unadjusted.p_value','Bonferonni.p_value',
               'baseline_R2','TF_R2','bg_avg_R2','tf_total-snps','num_of_right', 'bg_p-value','chr', 'tissue', 'adj_p-value']
    merged_df = merged_df[headers]
    merged_df.to_csv(outfname, sep="\t", header=True, index=False)
    print("find result:{:}".format(outfname))
    

