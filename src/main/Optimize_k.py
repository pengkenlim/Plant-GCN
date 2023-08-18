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
np.seterr(all="ignore")

if __name__ == "__main__":
        
        parser= argparse.ArgumentParser(description="Optimize_k.py.\n\
                                        Evaluate the performance of ensemble GCN across selected range of k.\
                                        Only edges in our dataset will have their ensemble weight edges calculated.")
        
        parser.add_argument("-w", "--workers", type=int, metavar="", default=4,
        help = "Number of workers for parallelization." )

        parser.add_argument("-o", "--ouput_dir", type=str, metavar="", required = True,
        help = "Directory to output. Same as Label_edges.py and Partition_expression_data.py" )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )

        parser.add_argument("-cc","--correlation_coefficient", type=str, default = "bicor", choices = ["PCC","SCC","bicor"] ,
        help = "Correlation coefficient to optimize for. \'bicor\' by default. NOTE: Only \'bicor\' is supported for now.")

        args=parser.parse_args()
        
        workers=args.workers
        ouput_dir=args.ouput_dir
        delimiter=args.delimiter
        correlation_coefficient = args.correlation_coefficient

        
        sub_outdir = os.path.join(ouput_dir , "Optimize_k")
        if delimiter == "t":
                delim = "\t"
        else:
                delim = ","
        #load data needed from workdirectory
        Tid2Gid_dict = read_write.load_pickle(os.path.join(ouput_dir , "Remove_isoforms" , "Tid2Gid_dict.pkl"))
        k_cluster_assignment_dict =  read_write.load_pickle(os.path.join(ouput_dir, "Partition_expression_data", "k_cluster_assignment_dict.pkl"))
        selected_k =  read_write.load_pickle(os.path.join(ouput_dir, "Partition_expression_data", "selected_k.pkl"))
        expmat_path = os.path.join(ouput_dir, "QC_expression_data",  f"expression_matrix.{delimiter}sv")
        positive_met_edges = read_write.load_pickle(os.path.join(ouput_dir, "Label_edges", "positive_met_edges.pkl"))
        negative_met_edges = read_write.load_pickle(os.path.join(ouput_dir, "Label_edges", "negative_met_edges.pkl"))

        negative_met_edges_unpacked = []
        for edges in negative_met_edges.values():
                negative_met_edges_unpacked.extend(edges)
        negative_met_edges_unpacked = list(set(negative_met_edges_unpacked))
        
        
        for cc in [correlation_coefficient]:
                #establish dataframes to store and write performance results
                AVG_full_df, AUC_ROC_full_df, AUC_PRC_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_summ_df, AUC_ROC_summ_df, AUC_PRC_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                cc_sub_outdir = os.path.join(sub_outdir ,cc)

                print(f"Calculating ensemble scores of edges across k range for correlation coefficient: \'{cc}\'")
                k_cluster_assignment_dict[1]= np.array([ 0 for i in range(len(k_cluster_assignment_dict[selected_k[0]]))])
                selected_k = [1]+ selected_k
                

                for k in selected_k:
                        k_sub_outdir = os.path.join(cc_sub_outdir, f"{k}_K")
                        read_write.establish_dir(k_sub_outdir, isdir =True)
                        
                        positive_met_edges_cor_path = os.path.join(k_sub_outdir, "positive_met_edges_cor.tsv")
                        negative_met_edges_cor_path = os.path.join(k_sub_outdir, "negative_met_edges_cor.tsv")
                        
                        if cc == "bicor":
                                optimize_k = bicor.optimize_k
                        elif cc == "PCC":
                                optimize_k = pearson.optimize_k
                        elif cc == "SCC":
                                optimize_k = spearman.optimize_k
                                
                        
                        optimize_k(k, positive_met_edges_cor_path, negative_met_edges_cor_path,expmat_path, Tid2Gid_dict,  
                                        k_cluster_assignment_dict, delim, workers, 
                                        positive_met_edges, negative_met_edges_unpacked)
                        
                        score_types = ['Max','Avg','RAvg','RWA','RRWA']
                        positive_edges_cor_dict, negative_edges_cor_dict , score_types = network.load_edges(positive_met_edges_cor_path,
                                                                                                        negative_met_edges_cor_path, 
                                                                                                        negative_met_edges, 
                                                                                                        score_types = score_types)
                        performance_dict = network_performance.evaluate(positive_edges_cor_dict, 
                                                                        negative_edges_cor_dict, 
                                                                        positive_met_edges , 
                                                                        negative_met_edges, 
                                                                        score_types)
                        read_write.to_pickle(performance_dict,  
                                        os.path.join(k_sub_outdir, "performance_dict.pkl"))
                        ##performance_dict = read_write.load_pickle(os.path.join(k_sub_outdir, "performance_dict.pkl"))
                        
                        network_performance.summarize_and_out(k_sub_outdir, score_types, performance_dict)

                        #store results in dataframe
                        
                        AVG_full_df = network_performance.cat_k_to_df(AVG_full_df, score_types, performance_dict, k, "AVG" , full = True)
                        AUC_ROC_full_df = network_performance.cat_k_to_df(AUC_ROC_full_df, score_types, performance_dict, k, "AUC_ROC" , full = True)
                        AUC_PRC_full_df = network_performance.cat_k_to_df(AUC_PRC_full_df, score_types, performance_dict, k, "AUC_PRC" , full = True)

                        AVG_summ_df = network_performance.cat_k_to_df(AVG_summ_df, score_types, performance_dict, k, "AVG" , full = False)
                        AUC_ROC_summ_df = network_performance.cat_k_to_df(AUC_ROC_summ_df, score_types, performance_dict, k, "AUC_ROC" , full = False)
                        AUC_PRC_summ_df = network_performance.cat_k_to_df(AUC_PRC_summ_df, score_types, performance_dict, k, "AUC_PRC" , full = False)

                        print(f"k= {k} completed.")
                
                print(f"Getting best k and worst k based on AVG(AUC_ROC,AUC_PRC) scores...")
                AVG_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_full_scores.csv"))
                AUC_ROC_full_df.to_csv(os.path.join(cc_sub_outdir, "AUC_ROC_full_scores.csv"))
                AUC_PRC_full_df.to_csv(os.path.join(cc_sub_outdir, "AUC_PRC_full_scores.csv"))
                
                AVG_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_summ_scores.csv"))
                AUC_ROC_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUC_ROC_summ_scores.csv"))
                AUC_PRC_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUC_ROC_summ_scores.csv"))

                result_string = network_performance.best_worst_k(AVG_summ_df, score_types)
                print(result_string)
                with open(os.path.join(cc_sub_outdir, "Best_worst_k.csv"), "w") as f:
                        f.write(result_string)

