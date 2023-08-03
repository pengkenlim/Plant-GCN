#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

import pandas as pd

def load_exp_mat(expmatfile, expmatsep = "\t", genecolname="target_id"):
    expmat_df = pd.read_csv(expmatfile, sep = expmatsep).set_index(genecolname)
    return expmat_df
