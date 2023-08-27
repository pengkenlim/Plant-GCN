#setting sys.path for importing modules
import os
import sys
parent_module = "/home/ken/Plant-GCN/src/" # remove
sys.path.insert(0, parent_module) # remove
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import numpy as np
import pandas as pd

from data_processing import read_write , expression_matrix

expmat_path_original = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/QC_expression_data/expression_matrix.tsv"
expmat_path_5k = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5k/QC_expression_data/expression_matrix.tsv"
expmat_path_5h = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/QC_expression_data/expression_matrix.tsv"

expmat_o =  expression_matrix.load(expmat_path_original, expmatsep = "\t")
sample_list = list(expmat_o.columns)
sample_5k = list(np.random.choice(sample_list, size = 5000, replace=False))
expmat_5k = expression_matrix.subset(expmat_o, sample_5k)
expression_matrix.write(expmat_5k, expmat_path_5k)

expmat_o =  expression_matrix.load(expmat_path_original, expmatsep = "\t")
sample_list = list(expmat_o.columns)
sample_5h = list(np.random.choice(sample_list, size = 500, replace=False))
expmat_5h = expression_matrix.subset(expmat_o, sample_5h)
expression_matrix.write(expmat_5h, expmat_path_5h)

performance_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5k/Optimize_k/PCC/1_K/performance_dict.pkl")
performance_dict.keys()
