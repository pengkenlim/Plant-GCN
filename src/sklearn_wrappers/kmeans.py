#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pickle

#User functions

def iterate_over_krange(data, k_list, sizemax=0, randomstate=42):
    """Run kmeans clustering over a range of K"""
    kmeans_kwargs = {"init": "k-means++", "n_init": 10,"max_iter": 10000,"random_state": randomstate} #remove random
    silhouette_coefficients = []
    k_cluster_assignment_dict={}
    centroids_dict={}
    for k in k_list:
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(data)
        silhouette_coefficients.append(silhouette_score(data, kmeans.labels_))
        centroids_dict[k]= kmeans.cluster_centers_
        k_cluster_assignment_dict[k]= kmeans.labels_
        print(f"K-means iteration at k={k} complete.SC:{silhouette_score(data, kmeans.labels_)}\n")
    return k_cluster_assignment_dict , silhouette_coefficients, centroids_dict

def select_k(silhouette_coefficients,k_cluster_assignment_dict):
    """"select Ks where silhoette coefficients peaks"""
    selected_list = []
    score = 0
    k_list = list(k_cluster_assignment_dict.keys())
    for k, sc in zip(k_list, silhouette_coefficients):
        if sc < score:
            selected_list.append(prev_k)
        else:
            prev_k = k
        score = sc
    selected_list.append(k)
    selected_list = sorted(list(set(selected_list)))
    return selected_list
     
#if __name__ == "__main__":
#    silhouette_coefficients= {0:0.1,
#                              1:0.2,
#                              2:0.3,
#                              3:0.2,
#                              4:0.3}
#    print(select_k(silhouette_coefficients)) #testing is the function works
