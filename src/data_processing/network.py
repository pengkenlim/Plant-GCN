#setting sys.path for importing modules
import os
import sys


if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

def load_full():
    network ={}
    return network


def load_edges(positive_edges_corr_path, negative_edges_corr_path, negative_edges, score_types = ['Max','Avg','RAvg','RWA','RRWA']):
    positive_edges_cor_dict = {}
    negative_edges_cor_dict_temp = {}
    negative_edges_cor_dict = {}

    with open(positive_edges_corr_path) as f:
        for line_no, line in enumerate(f):
            if line_no != 0 and line != "": # skip first line of headers and last line
                edge, *_, scores = line.split("\t")
                score_dict  = {type : float(score) for type , score in zip(score_types, scores.split("\n")[0].split(","))}
                positive_edges_cor_dict[edge] = score_dict
    
    with open(negative_edges_corr_path) as f:
        for line_no, line in enumerate(f):
            if line_no != 0 and line != "": # skip first line of headers and last line
                edge, *_, scores = line.split("\t")
                score_dict  = {type : float(score) for type , score in zip(score_types, scores.split("\n")[0].split(","))}
                negative_edges_cor_dict_temp[edge] = score_dict
    
    for ds, edges in negative_edges.items():
        negative_edges_cor_dict[ds] = {}
        for edge in edges:
            negative_edges_cor_dict[ds][edge] = negative_edges_cor_dict_temp[edge]
    return positive_edges_cor_dict, negative_edges_cor_dict , score_types



if __name__ == "__main__": #for testing functions
    sys.path.insert(0,"/home/ken/Plant-GCN/src")
    from data_processing import read_write
    negative_met_edges = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/taxid3880/Label_edges/negative_met_edges.pkl")
    networkdir = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/taxid3880/Optimize_k/bicor/65_K"
    positive_edges_cor_dict, negative_edges_cor_dict = load_edges(networkdir, negative_met_edges)