#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import pickle

#User functions

def kmeans_kwalk(data, k_min, k_max, sizemax=0, randomstate=42):
    """Run kmeans clustering over a range of K"""
    kmeans_kwargs = {"init": "k-means++", "n_init": 10,"max_iter": 10000,"random_state": randomstate} #remove random
    silhouette_coefficients = []
    k_cluster_assignment_dict={}
    centroids_dict={}
    for k in range(k_min,k_max):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(data)
        silhouette_coefficients.append(silhouette_score(data, kmeans.labels_))
        centroids_dict[k]= kmeans.cluster_centers_
        k_cluster_assignment_dict[k]= kmeans.labels_
        print(f"K-means iteration at k={k} complete.SC:{silhouette_score(data, kmeans.labels_)}\n")
    return k_cluster_assignment_dict , silhouette_coefficients, centroids_dict