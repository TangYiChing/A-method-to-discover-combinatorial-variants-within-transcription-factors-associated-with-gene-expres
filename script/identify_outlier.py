"""
"""
import sys
import argparse
import pandas as pd
import subprocess as sp

RPATH = '/usr/bin/Rscript'
SCRIPTPATH = '/repo4/ytang4/GTEx_v6/geuvadis_data/script/'

def read_as_df(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0, sep="\t")
    return df

def has_headers(df, headerList):
    notIn = []
    for header in headerList:
        if not header in df.columns:
            notIn.append(header)
    if len(notIn) != 0:
        print( 'ERROR: must have headers={:}, missing={:}'.format(headerList, notIn) )
        sys.exit(1)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Generate y table for outlierTest and scatter plot.")
    parser.add_argument("-tf", "--tf_resultFile",
                        required = True,
                        help = "modeling result from tf model. Must have headers=[]")
    parser.add_argument("-bl", "--bl_resultFile",
                        required = True,
                        help = "modeling result from baseline model. Must have headers=[]")
    parser.add_argument("-p", "--out_prefix",
                        required = True,
                        help = "Prefix of output file.")
    args = parser.parse_args()

    # initialize
    tfFile = args.tf_resultFile
    blFile = args.bl_resultFile
    prefix = args.out_prefix

    # load data
    tf = read_as_df(tfFile)
    bl = read_as_df(blFile)

    # check headers
    has_headers(tf, ['R2', 'genename', 'w_TF'])
    has_headers(bl, ['R2', 'genename'])

    # exclude genes that has zero TFs (i.e. w_TF = No) 
    wTF = tf['w_TF']=='yes'
    useTF = tf['pct_TF']>0
    tf_wTF = tf[wTF & useTF] #tf.loc[tf['w_TF']=='yes']
    
    # subsetting by finding shared genes
    sharedGenes = sorted(list(set(tf_wTF.index) & set(bl.index)))
    sub_tf = tf_wTF.loc[sharedGenes][['R2', 'genename']]
    sub_bl = bl.loc[sharedGenes][['R2']]
    sub_tf.columns = ['TF_R2', 'genename']
    sub_bl.columns = ['baseline_R2']

    # join two tabels (i.e. ytable)
    combinedDf = pd.concat([sub_bl, sub_tf], axis=1, join='inner')
    combinedDf.to_csv(prefix+'.ytable.txt', header=True, index=True, sep="\t")
    #print( 'tf={:}, bl={:}, combinedDf={:}'.format(sub_tf.shape[0], sub_bl.shape[0], combinedDf.shape[0]) )
            
    # find outliers
    cmd = RPATH + ' ' + SCRIPTPATH + 'find_outlier.R --ytable ' + prefix+'.ytable.txt' + ' --out_prefix ' + prefix 
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    #print(cmd)
    process.wait()
 
