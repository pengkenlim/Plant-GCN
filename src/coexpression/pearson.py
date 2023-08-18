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
import warnings


def precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter="\t", workers=2):
    assignment= k_cluster_assignment_dict[k]
    genes = []
    nominators_dict, denominators_dict, = {cluster:[] for cluster in range(k)}, {cluster:[] for cluster in range(k)}
    with open(expmat_path, 'r') as fin:
        for idx, line in enumerate(fin):
            if idx != 0:
                parts = line.rstrip().split(delimiter)
                all_values = np.array([float(i) for i in parts[1:]])
                Transcript_id = parts[0]
                if Transcript_id in Tid2Gid_dict.keys():
                    gene = Tid2Gid_dict[Transcript_id]
                    genes.append(gene)               
                    for cluster in range(k): # precalc job
                        cluster_values = all_values[assignment == cluster]
                        nomi = cluster_values - np.array([(np.sum(cluster_values)/len(cluster_values))])
                        denomi = np.sqrt(np.sum(nomi**2))
                        nominators_dict[cluster].append(nomi)
                        denominators_dict[cluster].append(denomi)
            if idx % 5000 == 0:
                print("k=", k,":", idx , "genes prepared.")
    gene_dict={}
    for idx , gene in enumerate(genes, start=0):
        gene_dict[gene] = idx
    for cluster in range(k):
        nominators_dict[cluster] = np.array(nominators_dict[cluster])
        denominators_dict[cluster] = np.array(denominators_dict[cluster])
    
    return genes, gene_dict, nominators_dict, denominators_dict



def calc_job_k(source_array, target_array, shared_nominators_dict, shared_denominators_dict, cluster):
    cor_values = np.sum(shared_nominators_dict[cluster][source_array] * shared_nominators_dict[cluster][target_array], axis=1)/(shared_denominators_dict[cluster][source_array]* shared_denominators_dict[cluster][target_array])
    #cor_values = np.sum(shared_nominators_dict[source_array] * shared_nominators_dict[target_array], axis=1)/(shared_denominators_dict[source_array]* shared_denominators_dict[target_array])
    return cor_values


def calc_targeted(k, path, edges , genes, gene_dict, nominators_dict, denominators_dict , workers = 2):
    manager = mp.Manager()
    shared_nominators_dict = manager.dict(nominators_dict)
    shared_denominators_dict = manager.dict(denominators_dict)
    #add header if creating new file
    with open(path, "w") as f:
        f.write("Edge\tcluster_cor\tMax,Avg,RAvg,RWA,RRWA\n")
    source_array , target_array= [], []

    for edge in list(edges):
        source, target = edge.split("-")
        source_array.append(gene_dict[source])
        target_array.append(gene_dict[target])

    ALL_cor_means = []
    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        results = [executor.submit(calc_job_k, source_array, target_array, shared_nominators_dict, shared_denominators_dict, cluster) for cluster in range(k)]
        #results = [executor.submit(calc_job_k, source_array, target_array, nominators_dict[cluster], denominators_dict[cluster], cluster) for cluster in range(k)]
        for f in cf.as_completed(results):
            cor_means = f.result()
            ALL_cor_means.append(cor_means)
    ALL_cor_means = np.array(ALL_cor_means)
    batch_ensemble_scores = []
    for mode in ["Max", "Avg", "RAvg", "RWA", "RRWA"]:
        batch_ensemble_scores.append(ensemble.aggregate(ALL_cor_means, mode, axis = 0))
    batch_ensemble_scores = np.array(batch_ensemble_scores)
    with open(path, "a") as f:
        for idx, edge in enumerate(edges):
            cluster_cor = ALL_cor_means[:,idx]
            ensemble_scores = batch_ensemble_scores[:,idx]
            cluster_cor = ",".join([str(i) for i in cluster_cor])
            ensemble_scores = ",".join([str(i) for i in ensemble_scores])
            f.write(f"{edge}\t{cluster_cor}\t{ensemble_scores}\n")


def calc_untargeted():
    pass
        
def optimize_k(k, positive_met_edges_cor_path, negative_met_edges_cor_path ,expmat_path, Tid2Gid_dict,  k_cluster_assignment_dict, delim, workers, positive_met_edges, negative_met_edges_unpacked):
       genes, gene_dict, nominators_dict, denominators_dict = precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)

       print("Calculating and writing correlations of positive edges...")

       calc_targeted(k, positive_met_edges_cor_path, positive_met_edges , genes, gene_dict, nominators_dict, denominators_dict ,workers = workers)
       print("Calculating and writing correlations of negative edges...")
       calc_targeted(k, negative_met_edges_cor_path, negative_met_edges_unpacked , genes, gene_dict, nominators_dict,denominators_dict , workers = workers)