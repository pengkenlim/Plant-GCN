#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse
from coexpression import bicor , pearson , spearman #for now only bicor and PCC is online
from data_processing import read_write, network
from analyses import network_performance
import numpy as np
import pandas as pd

if __name__ == "__main__":
        
        parser= argparse.ArgumentParser(description="Build_ensemble_GCN.py.")
        
        parser.add_argument("-w", "--workers", type=int, metavar="", default=4,
        help = "Number of workers for parallelization." )

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Directory to output. Same as Label_edges.py and Partition_expression_data.py" )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )

        parser.add_argument("-cc","--correlation_coefficient", type=str, default = "bicor", choices = ["PCC","SCC","bicor"] ,
        help = "Correlation coefficient to optimize for. \'bicor\' by default.")

        parser.add_argument("-am","--aggregation_method", type=str, required = True, choices = ["Max","Avg","RAvg","RWA","RRWA"] ,
        help = "Aggregation method to generate ensemble network.")

        parser.add_argument("-k","--k_clusters", type=int, default = 0,
        help = "Number of clusters to partition the expression data. For k=0, will use best k as determined in Optimize_k.py step")

        args=parser.parse_args()
        
        workers=args.workers
        output_dir=args.output_dir
        delimiter=args.delimiter
        correlation_coefficient = args.correlation_coefficient
        aggregation_method = args.aggregation_method
        k_clusters = args.k_clusters

        if delimiter == "t":
            delim = "\t"
        else:
            delim = ","

        if k_clusters ==0:
            with open(os.path.join( output_dir, "Optimize_k",correlation_coefficient, "Best_worst_k.csv"), "r") as f:
                for line in f:
                    if aggregation_method in line:
                        max_score , max_k = line.split("\t")[1:3]
                
            print(f"k=0 specified. Will use best k of {max_k} as determined by Optimize_K.py")
            k= int(max_k)
        else:
            k = k_clusters
        
        #establish output sub directories
        sub_outdir = os.path.join(output_dir , "Build_ensemble_GCN")
        cc_sub_outdir = os.path.join(sub_outdir, correlation_coefficient)
        
        

        network_dir_name = "_".join([correlation_coefficient,
                                    f"{k}k",
                                    aggregation_method])
        
        network_path = os.path.join(cc_sub_outdir, network_dir_name)

        read_write.establish_dir(network_path , isdir=True)
        print(f"Building ensemble network with these specifications:\n\n\
              Correlation coefficient= {correlation_coefficient}\n\
              k= {k}\n\
              Aggregation method= {aggregation_method}\n\n\
              written to {network_path}")
        
        #loading required data
        Tid2Gid_dict = read_write.load_pickle(os.path.join(output_dir , "Remove_isoforms" , "Tid2Gid_dict.pkl"))
        k_cluster_assignment_dict =  read_write.load_pickle(os.path.join(output_dir, "Partition_expression_data", "k_cluster_assignment_dict.pkl"))
        expmat_path = os.path.join(output_dir, "QC_expression_data",  f"expression_matrix.{delimiter}sv")

        #start
        if correlation_coefficient == "bicor":
            build_ensemble_GCN = bicor.build_ensemble_GCN
        elif correlation_coefficient== "PCC":
            build_ensemble_GCN = pearson.build_ensemble_GCN
        elif correlation_coefficient == "SCC":
            build_ensemble_GCN = spearman.build_ensemble_GCN

        build_ensemble_GCN( Tid2Gid_dict, k_cluster_assignment_dict, expmat_path, k, network_path, aggregation_method, delim, workers)

