#setting sys.path for importing modules
import os
import sys
parent_module = "/home/ken/Plant-GCN/src/"
sys.path.insert(0, parent_module)
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from data_processing import read_write

def parse_ara_annot(annot_path, All_genes, EVD_retain = ["IDA"], GO_cat = "BP"):
    annot_dict = {}
    if GO_cat == "BP":
        target_aspect = "P"
    elif GO_cat == "MF":
        target_aspect = "F" # add in         
    with open(annot_path ,  "r") as f:
        if EVD_retain != None:
            for line_no , line in enumerate(f):
                if line != "" and line_no >3 :
                    line = line.split("\t")
                    aspect = line[7]
                    EVD = line[9]
                    if aspect == target_aspect and EVD in EVD_retain:
                        GO = line[5]
                        gene = line[0].upper().split(".")[0]
                        if gene in All_genes:
                            annot_dict[GO] = annot_dict.get(GO, []) + [gene]
        else:
            for line_no , line in enumerate(f):
                if line != "" and line_no >3 :
                    line = line.split("\t")
                    aspect = line[7]
                    EVD = line[9]
                    if aspect == target_aspect:
                        GO = line[5]
                        gene = line[0].upper().split(".")[0]
                        if gene in All_genes:
                            annot_dict[GO] = annot_dict.get(GO, []) + [gene]

        return annot_dict

def annot_dict_to_edges(annot_dict):
    edges_set = set()
    GO_idx = 0
    for GO, genes in annot_dict.items():
        genes = list(set(genes)) # get rid of duplicates
        for idx, source in enumerate(genes):
            for target in genes[idx+1:]:
                edges_set.add("--".join(sorted([source, target])))
        GO_idx += 1
        if GO_idx % 100 ==0:
            print(f"Edges extracted for {GO_idx}/{len(annot_dict)} GO terms")
    return edges_set

if __name__ == "__main__":
    Tid2Gid_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_BS/Remove_isoforms/Tid2Gid_dict.pkl")
    All_genes = list(Tid2Gid_dict.values())
    positive_TF_edges = read_write.load_pickle( "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_TF_edges.pkl")
    positive_met_edges = read_write.load_pickle( "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_met_edges.pkl")
    annot_path = "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/ATH_GO_GOSLIM.txt"
    
    BP_IDA_annot_dict = parse_ara_annot(annot_path, All_genes, EVD_retain = ["IDA"], GO_cat = "BP")
    BP_IDA_GO_edges_set = annot_dict_to_edges(BP_IDA_annot_dict)
    
    MF_ALL_annot_dict = parse_ara_annot(annot_path, All_genes, EVD_retain = None, GO_cat = "MF")
    MF_ALL_GO_edges_set = annot_dict_to_edges(MF_ALL_annot_dict)
    positive_GO_edges = list(BP_IDA_GO_edges_set - MF_ALL_GO_edges_set - set(positive_TF_edges) - set(positive_met_edges))
    
    read_write.to_pickle(positive_GO_edges , "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_GO_edges.pkl")