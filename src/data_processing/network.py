#setting sys.path for importing modules
import os
import sys
import scipy
import math
import numpy as np

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

def calculate_edge_index(s, t, num_nodes):
    if s > t:
        s, t = t, s  # Ensure s < t
        return False, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1
    else:
        return True, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1

def init_flat_half_adj(num_edges):
    edge_values = np.empty(num_edges)
    edge_values[:] = math.nan
    return edge_values

def init_adj(num_genes):
    edge_values = np.empty( (num_genes, num_genes) )
    edge_values[:] = math.nan
    return edge_values


def extract_ranks(network_dir):
    genes = os.listdir(network_dir)
    gene_dict = {gene:i for i, gene in enumerate(genes)}
    num_genes = len(genes)
    num_edges = num_genes * (num_genes - 1) // 2
    upper_half = init_flat_half_adj(num_edges)
    lower_half = init_flat_half_adj(num_edges)
    for idx, source in enumerate(genes):
        with open(os.path.join(network_dir, source), "r") as f:
            for line_no, line in enumerate(f):
                if line_no != 0 and line != "": #skip first and last line
                    target, cor, Rank = line.split("\n")[0].split("\t")
                    if target != source:
                        source_idx = gene_dict[source]
                        target_idx = gene_dict[target]
                        upper, edge_index= calculate_edge_index(source_idx, target_idx, num_genes)
                        if upper:
                            upper_half[edge_index] = Rank
                        else:
                            lower_half[edge_index] = Rank
        if idx%1000 == 0:
            print(idx,"genes loaded")
    return genes , upper_half, lower_half, gene_dict



def write_MR_HRR(gene_dict, HRR_array, MR_array, old_network_dir, new_network_dir, genes):
    for idx, source in enumerate(genes):
        with open(os.path.join(old_network_dir, source), "r") as fin:
            with open(os.path.join(new_network_dir, source), "w") as fout:
                for line_no, line in enumerate(fin):
                    if line_no == 0:
                        fout.write(line.split("\n")[0] + "\tHRR\tMR\n")
                    elif line != "":
                        target, cor, Rank = line.split("\n")[0].split("\t")
                        cor = np.round(float(cor), 5)
                        if target == source:
                            MR , HRR = 1.0 , 1.0
                        else:
                            source_idx = gene_dict[source]
                            target_idx = gene_dict[target]
                            _, edge_index= calculate_edge_index(source_idx, target_idx, len(genes))
                            MR , HRR = MR_array[edge_index] , HRR_array[edge_index]
                            MR , HRR = np.round(MR, 3) , np.round(HRR, 3)
                        fout.write(f"{target}\t{cor}\t{Rank}\t{HRR}\t{MR}\n")
        if idx%1000 == 0:
            print(idx,"genes done")

########################################################
def load_full_network(network_dir,half=True):
    genes = os.listdir(network_dir)
    gene_dict = {gene:i for i, gene in enumerate(genes)}
    num_genes = len(genes)
    if half:
        num_edges = num_genes * (num_genes - 1) // 2
        network = {}
        network["RAW"] = init_flat_half_adj(num_edges)
        network["HRR"] = init_flat_half_adj(num_edges)
        network["MR"] = init_flat_half_adj(num_edges)

        for idx, source in enumerate(genes):
            with open(os.path.join(network_dir, source), "r") as f:
                for line_no, line in enumerate(f):
                    if line_no != 0 and line != "": #skip first and last line
                        target, cor, Rank , HRR, MR = line.split("\n")[0].split("\t")
                        if target != source:
                            source_idx = gene_dict[source]
                            target_idx = gene_dict[target]
                            upper, edge_index= calculate_edge_index(source_idx, target_idx, num_genes)
                            
                            if upper:
                                network["RAW"][edge_index] = float(cor)
                                network["HRR"][edge_index] = float(HRR)
                                network["MR"][edge_index] = float(MR)

            if idx%1000 == 0:
                print(idx,"genes loaded")
    else:
        network = {}
        network["RAW"] = init_adj(num_genes)
        network["HRR"] = init_adj(num_genes)
        network["MR"] = init_adj(num_genes)
        for idx, source in enumerate(genes):
            with open(os.path.join(network_dir, source), "r") as f:
                for line_no, line in enumerate(f):
                    if line_no != 0 and line != "": #skip first and last line
                        target, cor, Rank , HRR, MR = line.split("\n")[0].split("\t")

                        source_idx = gene_dict[source]
                        target_idx = gene_dict[target]

                        network["RAW"][target_idx, source_idx] = float(cor)
                        network["RAW"][source_idx , target_idx ] = float(cor)
                        network["HRR"][target_idx, source_idx] = float(HRR)
                        network["HRR"][source_idx, target_idx] = float(HRR)
                        network["MR"][target_idx, source_idx] = float(MR)
                        network["MR"][source_idx, target_idx] = float(MR)
            if idx%1000 == 0:
                print(idx,"genes loaded")
                            
    return genes , network, gene_dict

def load_full_network_ATTED(network_dir, ATTED_convert_dict, half=True):
    ATTED_genes = os.listdir(network_dir)
    gene_dict = {}
    for index, AT_gene in enumerate(ATTED_genes):
        try:
            gene = ATTED_convert_dict[AT_gene]
            gene_dict[gene] = index
        except:
            pass
    num_genes = len(gene_dict)
    if half:
        num_edges = num_genes * (num_genes - 1) // 2
        network = {}
        network["RAW"] = init_flat_half_adj(num_edges)
        for idx, source_ATTED in enumerate(ATTED_genes):
            #if idx < 1000:
            if True:
                try:
                    source = ATTED_convert_dict[source_ATTED]
                    with open(os.path.join(network_dir, source_ATTED), "r") as f:
                        for line_no, line in enumerate(f):
                            if line!="":
                                target_ATTED, cor = line.split("\n")[0].split("\t")
                                target = ATTED_convert_dict.get( target_ATTED,source)
                                if target != source:
                                    source_idx = gene_dict[source]
                                    target_idx = gene_dict[target]
                                    upper, edge_index= calculate_edge_index(source_idx, target_idx, num_genes)
                                    if upper:
                                        network["RAW"][edge_index] = float(cor)
                except:
                    pass
                
                if idx%1000 == 0:
                    print(idx,"genes loaded")
    else:
        network = {}
    return list(gene_dict.keys()) , network, gene_dict
########################################################

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
    #pass
    sys.path.insert(0,"/home/ken/Plant-GCN/src")
    from data_processing import read_write
    print("Loading network")
    genes , network_object , gene_dict = load_full_network("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Add_ranks/PCC_444k_RWA", half=False)
    read_write.establish_dir("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/H_clust/PCC_444k_RWA" , isdir=True)
    print("Saving network")
    read_write.to_pickle(gene_dict, "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/H_clust/PCC_444k_RWA/gene_dict.pkl")
    read_write.to_pickle(network_object, "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/H_clust/PCC_444k_RWA/network_object.pkl")
    #network_object = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/PCC_8k_Max/network_object.pkl")
    #gene_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/PCC_8k_Max/gene_dict.pkl")
    #performance_dict = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/PCC_8k_Max/performance_dict.pkl")
    #negative_met_edges = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/taxid3880/Label_edges/negative_met_edges.pkl")
    #networkdir = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/taxid3880/Optimize_k/bicor/65_K"
    #positive_edges_cor_dict, negative_edges_cor_dict = load_edges(networkdir, negative_met_edges)