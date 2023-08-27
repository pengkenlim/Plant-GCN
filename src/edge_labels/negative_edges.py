#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from data_processing import read_write

import numpy as np
def generate_dict(positive_edges, All_genes, iterations=100):
    negative_edges = {}
    positive_edges_set = set(positive_edges)
    for i in range(iterations):
        random_gene_array = np.random.choice(All_genes, size = (len(positive_edges_set)*2), replace=True) # draw x number of genes from all genes. x = 2*y where y is the number of edges
        negative_edges_set = set()
        
        #pair up adjacent genes within random gene list to form negative edges using index slicing.
        for i2 in range(len(positive_edges_set)):
            source_index = i2*2
            source = random_gene_array[source_index]
            target = random_gene_array[source_index + 1]
            if source == target: # unlikely as it might be, remove self pairing by using the gene further adjacent to the target as target
                target = random_gene_array[source_index + 2]
            negative_edges_set.add("--".join(sorted([source, target])))

        negative_edges_set = negative_edges_set - positive_edges_set #make sure negative and positive edges are mutually exclusive
        
        #compensate for edges removed.
        while len(negative_edges_set) < len(positive_edges_set):
             source_target = np.random.choice(All_genes, size = 2, replace=False)
             negative_edges_set.add("--".join(sorted(list(source_target))))

        negative_edges[i] = list(negative_edges_set)
        if i % 10 == 0:
            print(i,"/",iterations, "negative sample datasets generated")
    return negative_edges

def generate_dict_via_subset(negative_All_edges, positive_edges, iterations = 100):
    negative_edges = {}
    for i in range(iterations):
        negative_edges[i] = list(np.random.choice(negative_All_edges[i], size = len(positive_edges), replace=False))
    return negative_edges

def edge_editor(dataset_edges , toreplace, replacement):
    if type(dataset_edges) == list:
        new_dataset_edges = [edge.replace(toreplace, replacement) for edge in dataset_edges]
        return new_dataset_edges
    elif type(dataset_edges) == dict:
        new_dataset_edges_dict = {}
        for key , edge_list in dataset_edges.items():
            new_dataset_edges_dict[key] = [edge.replace(toreplace, replacement) for edge in edge_list]
            return new_dataset_edges_dict

if __name__ == "__main__":
    Tid2Gid_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_BS/Remove_isoforms/Tid2Gid_dict.pkl")
    All_genes = list(Tid2Gid_dict.values())
    positive_TF_edges = read_write.load_pickle( "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_TF_edges.pkl")
    positive_met_edges = read_write.load_pickle( "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_met_edges.pkl")
    positive_GO_edges = read_write.load_pickle( "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_GO_edges.pkl")
    
    positive_All_edges  = list(set(positive_met_edges + positive_TF_edges + positive_GO_edges ))
    read_write.to_pickle(positive_All_edges,"/home/ken/Plant-GCN/test_data/Arabidopsis_edges/positive_All_edges.pkl" )
    
    negative_All_edges = generate_dict(positive_All_edges, All_genes, iterations=100)
    
    negative_met_edges = generate_dict_via_subset(negative_All_edges, positive_met_edges ,iterations=100)
    negative_TF_edges = generate_dict_via_subset(negative_All_edges, positive_TF_edges ,iterations=100)
    negative_GO_edges = generate_dict_via_subset(negative_All_edges, positive_GO_edges ,iterations=100)
    
    read_write.to_pickle(negative_All_edges,"/home/ken/Plant-GCN/test_data/Arabidopsis_edges/negative_All_edges.pkl" )
    read_write.to_pickle(negative_met_edges,"/home/ken/Plant-GCN/test_data/Arabidopsis_edges/negative_met_edges.pkl" )
    read_write.to_pickle(negative_GO_edges,"/home/ken/Plant-GCN/test_data/Arabidopsis_edges/negative_GO_edges.pkl" )
    read_write.to_pickle(negative_TF_edges,"/home/ken/Plant-GCN/test_data/Arabidopsis_edges/negative_TF_edges.pkl" )
    
    for edge in list(positive_All_edges):
        source, target = edge.split("--")