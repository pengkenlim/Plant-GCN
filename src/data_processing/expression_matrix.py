#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import pandas as pd
from scipy.stats import iqr
import numpy as np

def load(expmatfile, expmatsep = "\t", indexcolname="target_id"):
    #expmat_df = pd.read_csv(expmatfile, sep = expmatsep).set_index(indexcolname)
    expmat_df = pd.read_csv(expmatfile, sep = expmatsep, index_col=0)
    return expmat_df

def load_transposed(expmatfile, expmatsep = "\t", indexcolname="accession"):
    """"For expression matrix that is outputed by ken's pipeline where cols = genes and each row correspond to on sample"""
    #expmat_df = pd.read_csv( expmatfile, sep=expmatsep).set_index(indexcolname)
    expmat_df = pd.read_csv(expmatfile, sep = expmatsep, index_col=0)
    expmat_df = expmat_df[~expmat_df.index.duplicated(keep='first')].transpose() #get rid of duplicated sample rows then transpose
    return expmat_df

def subset(expmat_df, samples_to_keep):
     samples_to_drop = set(list(expmat_df.columns)) - set(samples_to_keep)
     samples_to_drop = list(samples_to_drop)
     expmat_df.drop(labels = samples_to_drop, inplace=True, axis = 'columns')
     return expmat_df

def write(expmat_df, path, expmatsep = "\t"):
     expmat_df.to_csv(path , sep = expmatsep)


def lowerfence_iqr_cutoff(mappingvalues):
    '''replaces kdecutoff (more reliable and is a well established statistical method to look for outliers). Basically returns the lower fence (Q1 - 1.5*IQR) using the interquartile range (IQR)
    if lower fence is lower than 20%, return 10 % instead'''
    return max([np.percentile(mappingvalues, 25) - 1.5*iqr(mappingvalues), 20])

def thresholder(maprate_dict, cutoff):
    """wrapper function for lowerfence_iqr_cutoff(). takes in the dictionary of {sample:PS%,...} in log file.
    upacks the dictionary and return lists of total, failed and passed sample as well"""
    if cutoff==0:
        cutoff= lowerfence_iqr_cutoff(list(maprate_dict.values()))
    total= list(maprate_dict.keys())
    failed= [sample for sample, maprate in maprate_dict.items() if maprate < cutoff or maprate == cutoff]
    passed= [sample for sample, maprate in maprate_dict.items() if maprate > cutoff ]
    return total, failed, passed, cutoff

def load_qc(path):
    maprate_dict = {}
    with open(path, "r") as f:
         for line in f:
              if line != "":
                    line_contents = line.split("\n")[0].split("\t")
                    sample, p_PS ,n_PS= line_contents
                    try:
                        p_PS, n_PS = float(p_PS), float(n_PS)
                        if sample not in maprate_dict.keys():
                            maprate_dict[sample] = p_PS
                    except:
                        pass
    return maprate_dict

