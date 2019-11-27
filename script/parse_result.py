"""
"""
import os
import glob 
import argparse
import pandas as pd

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def concate_files(fileList):
    all_dfs = []
    ttl_rows = 0
    for f in fileList:
        df = read_as_df(f)
        all_dfs.append(df)
        ttl_rows+=df.shape[0]
    combinedDf = pd.concat(all_dfs)

    if combinedDf.shape[0] != ttl_rows:
        print( 'ERROR: number of rows does not match: should be{:}, got{:}'.format(ttl_rows, combinedDf.shape[0]) )

    return combinedDf

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse model result. Given the path of result folders, return a concatenated result file of chr1 to chr22")
    parser.add_argument("-f", "--result_path",
                        required = True,
                        help = "Folder name string representing output folders.")
    parser.add_argument("-m", "--model_name",
                        choices = ['baseline', 'tf-binding', 'tf-both', 'tf-regulation'],
                        required = True,
                        help = "Choices of model name: [ baseline | tf-binding | tf-both | tf-regulation ]")
    parser.add_argument("-p", "--out_prefix",
                        nargs = "?",
                        default = "ModelResult",
                        help = "Prefix of output file.")
    args = parser.parse_args()


    # initialize
    resultPath = args.result_path
    modelName = args.model_name
    prefix = args.out_prefix

    # parse results and concate them together 
    resultFileList = glob.glob(resultPath+'/working*.txt') 
    df = concate_files(resultFileList)
    print( 'prefix={:}, model={:}, # result files={:}'.format(prefix, modelName, len(resultFileList)) )
    
    # save to file
    df.to_csv(prefix+'.'+modelName+'.result.txt', header=True, index=True, sep="\t")
