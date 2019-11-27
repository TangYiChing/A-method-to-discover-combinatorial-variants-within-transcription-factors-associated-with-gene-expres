"""
Run Wilconsin Ranked Sum for each outliers that passed background models
Report outliers passing Wilconsin paired test (p<0.01)
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import scipy.stats as scistats

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
    missing_header = []
    for header in header_list:
        if not header in header_in_df:
            missing_header.append(header)
        
    if len(missing_header) > 0:
        print("header should have:{:}, missing={:}".format(header_list, missing_header))
        sys.exit(1)
    return

def collect_workingFiles(resultFolder, chrList):
    # collecting results from all chromsomes
    files = glob.glob(resultFolder+'working*.txt')
    resultDict = {'chr'+str(chr):[] for chr in range(1,23)}
    for f in files:
        # retrive chr number
        nchr = os.path.basename(f).split('.txt')[0].split('chr')[-1]  
        # append to dict  
        if int(nchr) in chrList:
            chr_key_str = 'chr' + str(nchr)
            resultDict[chr_key_str].append(f)

    return resultDict

def chrDict2runDict(baselineDict, tfDict, chrStr, run=100):
    runDict = {'run'+str(run+1):[] for run in range(int(run))}
    for f in baselineDict[chrStr]:
        run_key_str = os.path.basename(f).split('_')[-2].split('.')[-1]
        runDict[run_key_str].append(f)
    for f in tfDict[chrStr]:
        run_key_str = os.path.basename(f).split('_')[-2].split('.')[-1]
        runDict[run_key_str].append(f)

    return runDict
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Parse random fold selection result.")
    parser.add_argument("--outlierFile",
                        required = True,
                        help = "result file from parse_backgroundmodel.py Must have headers = [adj_p-value]")
    parser.add_argument("--baseline_resultFolder",
                        required = True,
                        help = "Location to results from random fold selection model. (e.g., ./baseline/randomSelection/)")
    parser.add_argument("--tfmodel_resultFolder",
                        required = True,
                        help = "Location to results from random fold selection model.")
    parser.add_argument("--model_name",
                        choices = ['tf-binding', 'tf-both', 'tf-regulation'],
                        help = "Choice of TF-model: [ tf-binding | tf-both | tf-regulation ]")
    parser.add_argument("--threshold",
                        required = False,
                        type = float, 
                        default = 0.01,
                        help = "threshold to filter outliers. default: 0.01")
    parser.add_argument("--run",
                        type = int,
                        help = "the number of files from N runs.")
    parser.add_argument("--out_prefix",
                        nargs = "?",
                        default = "result",
                        help = "Prefix of output file.")
    args = parser.parse_args()

    # initialize
    outlierFile = args.outlierFile
    baselineFolder = args.baseline_resultFolder
    tfFolder = args.tfmodel_resultFolder
    model = args.model_name
    threshold = args.threshold
    runs = args.run
    prefix = args.out_prefix

    # get significant outliers (those passed background model)
    outliers = read_as_df(outlierFile)
    has_headers(outliers, ['chr', 'gene', 'genename', 'adj_p-value'])
    sig_outliers = outliers.loc[outliers['adj_p-value'] < threshold]
    print( 'model={:}, outliers={:}, threshold={:}, significant outliers={:}'.format(model, outliers.shape[0], threshold, sig_outliers.shape[0]) )

    # retrieve result from baseline and tfmodel
    chrList = list(sig_outliers['chr'].unique())
    baseline_fileDict = collect_workingFiles(baselineFolder, chrList)
    tfmodel_fileDict = collect_workingFiles(tfFolder, chrList)
    
    # retrieve paired-result from all runs, for each outliers
    qualifiedGeneList = []
    wilcoxinPvalue = []
    for chrom in chrList:
        # get chromsome and corresponding run files
        chromDf = sig_outliers.loc[sig_outliers['chr']==chrom]
        runDict = chrDict2runDict(baseline_fileDict, tfmodel_fileDict, 'chr'+str(chrom), run=100)
        # for each gene/outlier, combine run files to calculate Wilconsin ranked sum
        if chromDf.shape[0]>0:
            print( 'chr={:}, #outliers={:}'.format(chrom, chromDf.shape[0]) )
            genes = list( chromDf['gene'].unique() )
            baselineR2 = []
            tfR2 = []
            for i in range(len(genes)):
                for key, value in runDict.items(): # loop through all runs
                    bdf = read_as_df(value[0]) # files from baseline
                    tdf = read_as_df(value[1]) # files from tfmodel
                    if bdf.shape[0] > 0 and tdf.shape[0] > 0:
                        if genes[i] in bdf['gene'].values.tolist() and genes[i] in tdf['gene'].values.tolist():
                            #print( '    {:}/{:} outlier(s), gene={:}'.format(i+1, len(genes), genes[i]) )
                            bdf.set_index('gene', inplace=True)
                            tdf.set_index('gene', inplace=True)
                            # get R2 from each run
                            bR2 = bdf.loc[genes[i], 'R2']
                            tR2 = tdf.loc[genes[i], 'R2']
                            baselineR2.append(bR2)
                            tfR2.append(tR2)
                # filter by threshold to obtain qualified outlier(s)
                if len(baselineR2) > 0 or len(tfR2) >0:
                    stats, pvalue = scistats.wilcoxon(baselineR2, tfR2)
                    if pvalue < threshold:
                        qualifiedGeneList.append(genes[i])
                        wilcoxinPvalue.append(pvalue)
                        print( '    gene={:}, #bR2={:}, #tR2={:}, wilcoxin test: p-value={:}'.format(genes[i], len(baselineR2), len(tfR2), pvalue) )
    # combine back to outliers table
    df = sig_outliers.loc[sig_outliers['gene'].isin(qualifiedGeneList)]
    df.loc[:, 'wilcoxon p-value'] = wilcoxinPvalue
    print( 'threshold={:}, #pass={:}'.format(threshold, df.shape[0]))
    # save to file
    df.to_csv(prefix+'.'+model+'.outliers.pass.randomFoldSelection.txt', header=True, index=False, sep="\t")
