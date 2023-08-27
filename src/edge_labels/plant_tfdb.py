#setting sys.path for importing modules
import os
import sys
parent_module = "/home/ken/Plant-GCN/src/" # remove
sys.path.insert(0, parent_module) # remove
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import pandas as pd
from data_processing import read_write

def parse_ara_annot(annot_path, All_genes, reg_type = ["A"]):
    TF_edges_set  = set()
    TF_annot_df = pd.read_csv(annot_path)
    for idx, line in TF_annot_df.iterrows():
        if line['Activate/Repress'] in reg_type:
            source = line['TF ID']
            target = line["Target ID"]
            if source in All_genes and target in All_genes:
                TF_edges_set.add("--".join(sorted([source, target])))
    return TF_edges_set

if __name__ == "__main__":
    Tid2Gid_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_BS/Remove_isoforms/Tid2Gid_dict.pkl")
    All_genes = list(Tid2Gid_dict.values())
    annot_path = "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/Regulations_in_ATRM.csv"
    positive_TF_edges = list(parse_ara_annot(annot_path, All_genes, reg_type = ["A"]))
    read_write.to_pickle(positive_TF_edges, "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_TF_edges.pkl")

    