#!/data/anaconda3/bin/python3.6

"""
generate expression (gene_id by sample_id) and genotype (Id by sample_id)

"""

import os
import sys
import random
import collections
import numpy as np
import pandas as pd
import allel
import h5py
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(filename)s %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S')
logger = logging.getLogger(__name__)


if len(sys.argv) < 3:
    print("usage: generate_expression_genotype.py expression.txt genotype.vcf outPath.str")
    sys.exit(1)
else:
    prefix_expression = os.path.basename(sys.argv[1]).split(".txt")[0]
    prefix_genotype = os.path.basename(sys.argv[2]).split(".vcf")[0]
    outPath = sys.argv[3]

    expression = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0) # rows=gene_id, cols=sample_id
    allel.vcf_to_hdf5(sys.argv[2], prefix_genotype+'.h5', fields=['CHROM', 'samples', 'ID', 'GT'], overwrite=True)
    genotype_data = h5py.File(prefix_genotype+'.h5', mode='r')

# find shared ids  
id_dict = collections.OrderedDict()
exp_ids = expression.columns.tolist()
gen_ids = list(genotype_data['samples'])
num_exp_ids = len(exp_ids)
num_gen_ids = len(gen_ids)
for gid in gen_ids:
    for eid in exp_ids:
        if eid in gid:
            if not gid in id_dict:
                id_dict[gid] = eid
            else:
                print("found duplicates", gid, eid)

logging.debug("expression file: {:}".format(sys.argv[1]))
logging.debug("genotype file: {:}".format(sys.argv[2]))
logging.debug("# of ids (samples) in expression: {:}".format(num_exp_ids))
logging.debug("# of ids (samples) in genotype: {:}".format(num_gen_ids))
logging.debug("# of ids (samples) in both: {:}".format(len(list(id_dict.keys()))))

# extract expression to include samples in shared ids
shared_exp_ids = list(id_dict.values())
expression = expression[shared_exp_ids]
expression.to_csv(outPath+"/"+prefix_expression+".sharedsamples.txt", sep="\t", header=True, index=True)

# extract genotype to include samples in shared ids
gt = allel.GenotypeChunkedArray(genotype_data['calldata/GT'])
gt_array = gt.to_n_alt(fill=-1) # missing genotypes are denoted as -1
SnpSample_df = pd.DataFrame(data=list(gt_array),
                            index=genotype_data['variants/ID'],
                            columns=genotype_data['samples'])
shared_gen_ids = list(id_dict.keys())
genotype = SnpSample_df[shared_gen_ids]
genotype.columns = [gid.split('-')[0]+"-"+gid.split('-')[1] for gid in list(genotype.columns)] # modify GTExID
genotype.index.name = "Id"
genotype.to_csv(outPath+"/"+prefix_genotype+".sharedsamples.txt", sep="\t", header=True, index=True)

#logging.debug("find output at: {:}, {:}".format(prefix_expression+"sharedsamples.txt", prefix_genotype+"sharedsamples.txt"))
logging.debug("expression: {:}".format(expression.shape))
logging.debug("genotype: {:}".format(genotype.shape))
