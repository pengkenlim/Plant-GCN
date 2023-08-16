#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import numpy as np
import math
import concurrent.futures as cf
import multiprocessing as mp
from coexpression import ensemble



def calc_job(source_array, target_array, norm_weights_dict, cluster):
    np.seterr(all="ignore") # to remove this pesky runtime warning
    cor_values = np.einsum("ijk, ijk -> ij", norm_weights_dict[cluster][:,source_array,:], norm_weights_dict[cluster][:,target_array, :])
    cor_means = np.nanmean(cor_values, axis=0)
    return cor_means


def get_norm_weights(cluster, all_values, assignment):
    cluster_values = all_values[assignment == cluster]
    subvalues_D = np.array([cluster_values]) # shape is (1, len(cluster_values)) # so that it can handle 2D input in the future
    values_minus_med_D = subvalues_D - np.nanmedian(subvalues_D, axis=1).reshape(-1,1)
    UI_D= (values_minus_med_D)/(9*np.nanmedian(np.abs(values_minus_med_D), axis=1).reshape(-1,1))
    Identity_D = np.where(abs(UI_D) < 1, 1,0)
    nom_D = (values_minus_med_D)*((1-UI_D**2)**2)* Identity_D
    norm_weights = nom_D/np.sqrt(np.sum(nom_D**2, axis =1)).reshape(-1,1)
    return cluster , norm_weights
                

def precalc_job(idx, line, delimiter, k,Tid2Gid_dict, assignment):
    norm_weights_gene_dict = {cluster:[] for cluster in range(k)}
    parts = line.rstrip().split(delimiter)
    all_values = [float(i) for i in parts[1:]]
    Transcript_id = parts[0]
    if Transcript_id in Tid2Gid_dict.keys():
        gene = Tid2Gid_dict[Transcript_id]
        #genes.append(gene)
        all_values=np.array(all_values)
        for cluster in range(k):
            cluster , norm_weights = get_norm_weights(cluster, all_values, assignment)
            norm_weights_gene_dict[cluster] = norm_weights
        if idx % 5000 ==0:
                    print("k=", k,":", idx , "genes prepared.")
        return gene, norm_weights_gene_dict

def precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter="\t", workers=2):
    assignment= k_cluster_assignment_dict[k]
    genes = []
    norm_weights_dict = {cluster:[] for cluster in range(k)}
    with open(expmat_path, 'r') as fin:
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            results = [executor.submit(precalc_job, idx, line, delimiter, k,Tid2Gid_dict, assignment) for idx, line in enumerate(fin) if idx != 0]
            for f in cf.as_completed(results):
                gene, norm_weights_gene_dict = f.result()
                for cluster in range(k):
                    norm_weights_dict[cluster].append(norm_weights_gene_dict[cluster])
                genes.append(gene)
    print("k=", k,":","Transposing normaized weights...")
    for cluster in range(k):
        norm_weights_dict[cluster]= np.transpose(np.array(norm_weights_dict[cluster]), (1,0,2))
        print(f"Cluster {cluster} transposed.")
    gene_dict={}
    for idx , gene in enumerate(genes, start=0):
        gene_dict[gene] = idx
    return genes, gene_dict , norm_weights_dict



def calc_targeted(k, path, edges , genes, gene_dict, norm_weights_dict, workers = 2):
    manager = mp.Manager()
    shared_norm_weights_dict = manager.dict(norm_weights_dict)
    #add header if creating new file
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write("Edge\tcluster_cor\tMax,Avg,RAvg,RWA,RRWA\n")
    source_array , target_array= [], []
    calculated_edges =[]
    failed_edges = []
    for edge in list(edges):
        source, target = edge.split("-")
        if source in genes and target in genes:
            source_array.append(gene_dict[source])
            target_array.append(gene_dict[target])
            calculated_edges.append(edge)
        else:
            failed_edges.append(edge)
    ALL_cor_means = []
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        results = [executor.submit(calc_job, source_array, target_array, shared_norm_weights_dict, cluster) for cluster in range(k)]
        for f in cf.as_completed(results):
            cor_means = f.result()
            ALL_cor_means.append(cor_means)
    ALL_cor_means = np.array(ALL_cor_means)
    
    batch_ensemble_scores = []
    for mode in ["Max", "Avg", "RAvg", "RWA", "RRWA"]:
        batch_ensemble_scores.append(ensemble.aggregate(ALL_cor_means, mode, axis = 0))
    batch_ensemble_scores = np.array(batch_ensemble_scores)
    with open(path, "a") as f:
        for idx, edge in enumerate(calculated_edges):
            cluster_cor = ALL_cor_means[:,idx]
            ensemble_scores = batch_ensemble_scores[:,idx]
            cluster_cor = ",".join([str(i) for i in cluster_cor])
            ensemble_scores = ",".join([str(i) for i in ensemble_scores])
            f.write(f"{edge}\t{cluster_cor}\t{ensemble_scores}\n")
        
        for edge in failed_edges:
            cluster_cor = np.array([math.nan for i in range(k)])
            ensemble_scores = np.array([math.nan for i in ["Max", "Avg", "RAvg", "RWA", "RRWA"]])
            cluster_cor = ",".join([str(i) for i in cluster_cor])
            ensemble_scores = ",".join([str(i) for i in ensemble_scores])
            f.write(f"{edge}\t{cluster_cor}\t{ensemble_scores}\n")






def calc_untargeted():
    pass
        
def optimize_k(k, positive_met_edges_cor_path, negative_met_edges_cor_path ,expmat_path, Tid2Gid_dict,  k_cluster_assignment_dict, delim, workers, positive_met_edges, negative_met_edges_unpacked):
       genes, gene_dict , norm_weights_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
       print("Calculating and writing correlations of positive edges...")

       calc_targeted(k, positive_met_edges_cor_path, positive_met_edges , genes, gene_dict, norm_weights_dict, workers = workers)
       print("Calculating and writing correlations of negative edges...")
       calc_targeted(k, negative_met_edges_cor_path, negative_met_edges_unpacked , genes, gene_dict, norm_weights_dict, workers = workers)



