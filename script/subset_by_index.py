#!/data/anaconda3/bin/python3.6

"""
generate a subset of file

"""

import os
import sys
import random
import collections
import numpy as np
import pandas as pd
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(filename)s %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S')
logger = logging.getLogger(__name__)


if len(sys.argv) < 4:
    print("usage: generate_expression_genotype.py data_file.txt data_index_pos wanted_str_file.txt row_col.int output_path.str")
    sys.exit(1)
else:
    #get args
    prefix_file = os.path.basename(sys.argv[1]).split(".txt")[0]
    data_df = pd.read_csv(sys.argv[1], sep="\t", header=0)
    data_index = int(sys.argv[2])
    str_df = pd.read_csv(sys.argv[3], sep="\t", header=0)
    row_col_int = int(sys.argv[4])
    outPath_str = sys.argv[5]
   

# get wanted list from wanted_str_file.txt
n_cols = len(str_df.columns.tolist())
if n_cols > 1:
    print("found more than 1 column, only the 1st col will be used.")  
col_name = str_df.columns.tolist()[0]
wanted_list = str_df[col_name].tolist()

# set data index
data_index_label = data_df.columns.tolist()[data_index]


if row_col_int == 0:
    #this means the list of str_df can be found at rows of data_df
    #subset_df = data_df.copy().loc[wanted_list]
    subset_df = data_df.loc[data_df[data_index_label].isin(wanted_list)]
 
    # logging
    num_cols_data_df = len(data_df.index.tolist())
    num_cols_subset_df = len(subset_df.index.tolist())
    logging.debug("# of rows={:} in data file".format(num_cols_data_df))
    logging.debug("# of rows={:} in str file".format(num_cols_subset_df))


elif row_col_int == 1:
    #this means the list of str_df can be found at cols of data_df
    subset_df = data_df[wanted_list]

    # logging
    num_cols_data_df = len(data_df.columns.tolist())
    num_cols_subset_df = len(subset_df.columns.tolist())
    logging.debug("# of cols={:} in data file".format(num_cols_data_df))
    logging.debug("# of cols={:} in str file".format(num_cols_subset_df))


else:
    print("index col should be 1 or 0, given={:}".format(row_col_int))
    sys.exit(1)

# reset index
subset_df.to_csv(outPath_str+"/"+prefix_file+".subset.txt", sep="\t", header=True, index=False)
logging.debug("found output file at={:}".format(prefix_file+".subset.txt"))


