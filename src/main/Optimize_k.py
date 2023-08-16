#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse
from coexpression import bicor #for now only bicor is online
from data_processing import read_write, network
from analyses import network_performance
import numpy as np
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

        parser.add_argument("-cc","--correlation_coefficient", type=str, default = "bicor", choices = ["PCC","SCC","bicor","all"] ,
        help = "Correlation coefficient to optimize for. \'bicor\' by default. \'all\' to optimize for all correlation coeffiecients. NOTE: Only \'bicor\' is supported for now.")

        args=parser.parse_args()
        
        workers=args.workers
        ouput_dir=args.ouput_dir
        delimiter=args.delimiter
        correlation_coefficient = args.correlation_coefficient
        if correlation_coefficient == "all":
                correlation_coefficients = ["PCC","SCC","bicor"]
        else:
                correlation_coefficients = [correlation_coefficient]
        
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
        
        
        for cc in correlation_coefficients:
                print(f"Calculating ensemble scores of edges across k range for correlation coefficient: \'{cc}\'")
                k_cluster_assignment_dict[1]= np.array([ 0 for i in range(len(k_cluster_assignment_dict[selected_k[0]]))])
                selected_k = [1]+ selected_k
                for k in selected_k:
                        k_sub_outdir = os.path.join(sub_outdir ,cc, f"{k}_K")
                        read_write.establish_dir(k_sub_outdir, isdir =True)
                        
                        
                        bicor.optimize_k(k, k_sub_outdir, expmat_path, Tid2Gid_dict,  
                                        k_cluster_assignment_dict, delim, workers, 
                                        positive_met_edges, negative_met_edges_unpacked)
                        
                        positive_edges_cor_dict, negative_edges_cor_dict , score_types = network.load_edges(positive_met_edges_cor_path,
                                                                                                        negative_met_edges_cor_path, 
                                                                                                        negative_met_edges, 
                                                                                                        score_types = ['Max','Avg','RAvg','RWA','RRWA'])
                        performance_dict = network_performance.evaluate(positive_edges_cor_dict, 
                                                                        negative_edges_cor_dict, 
                                                                        positive_met_edges , 
                                                                        negative_met_edges, 
                                                                        score_types)
                        read_write.to_pickle(performance_dict,  
                                        os.path.join(k_sub_outdir, "performance_dict.pkl"))
                        
                        network_performance.summarize_and_out(k_sub_outdir, score_types, performance_dict)
                        print(f"k= {k} completed.")

                        #genes, gene_dict , norm_weights_dict = bicor.precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
                        
                        #print("Calculating and writing correlations of positive edges...")
                        #positive_met_edges_cor_path = os.path.join(k_sub_outdir, "positive_met_edges_cor.tsv")
                        #bicor.calc_targeted(k, positive_met_edges_cor_path, positive_met_edges , genes, gene_dict, norm_weights_dict, workers = workers)
                        
                        #print("Calculating and writing correlations of negative edges...")
                        #negative_met_edges_cor_path = os.path.join(k_sub_outdir, "negative_met_edges_cor.tsv")
                        #bicor.calc_targeted(k, negative_met_edges_cor_path, negative_met_edges_unpacked , genes, gene_dict, norm_weights_dict, workers = workers)
                        
                        #positive_edges_cor_dict, negative_edges_cor_dict , score_types = network.load_edges(positive_met_edges_cor_path, negative_met_edges_cor_path, negative_met_edges, score_types = ['Max','Avg','RAvg','RWA','RRWA'])
                        #performance_dict = network_performance.evaluate(positive_edges_cor_dict, negative_edges_cor_dict, positive_met_edges , negative_met_edges, score_types)
                        
                        #read_write.to_pickle(performance_dict,  os.path.join(k_sub_outdir, "performance_dict.pkl"))
                        #with open(os.path.join(k_sub_outdir, "performance_summary.tsv"), "w") as f:
                        #        f.write(f"score_type\tAVG(AUC_ROC,AUC_PRC)\tAUC_ROC\tAUC_PRC\n")
                        #        for score_type in score_types:
                        #                f.write(f"{score_type}\t{performance_dict[score_type]['Quartiles']['AVG'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_ROC'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_PRC'][1]}\n")
                
                #print(f"Calculating ensemble scores of edges across k range for correlation coefficient: \'{cc}\'")
                #for k in selected_k:
                #        k_sub_outdir = os.path.join(sub_outdir ,cc, f"{k}_K")
                #        read_write.establish_dir(k_sub_outdir, isdir =True)
                #        genes, gene_dict , norm_weights_dict = bicor.precalc(expmat_path, Tid2Gid_dict, k_cluster_assignment_dict, k, delimiter=delim, workers=workers)
                #        
                #        print("Calculating and writing correlations of positive edges...")
                #        positive_met_edges_cor_path = os.path.join(k_sub_outdir, "positive_met_edges_cor.tsv")
                #        bicor.calc_targeted(k, positive_met_edges_cor_path, positive_met_edges , genes, gene_dict, norm_weights_dict, workers = workers)
                #       
                #        print("Calculating and writing correlations of negative edges...")
                #        negative_met_edges_cor_path = os.path.join(k_sub_outdir, "negative_met_edges_cor.tsv")
                #        bicor.calc_targeted(k, negative_met_edges_cor_path, negative_met_edges_unpacked , genes, gene_dict, norm_weights_dict, workers = workers)
                
                #print(f"Evaluating performance across k range for correlation coefficient: \'{cc}\'")
                #for k in selected_k:
                #        k_sub_outdir = os.path.join(sub_outdir ,cc, f"{k}_K")
                #        positive_met_edges_cor_path = os.path.join(k_sub_outdir, "positive_met_edges_cor.tsv")
                #        negative_met_edges_cor_path = os.path.join(k_sub_outdir, "negative_met_edges_cor.tsv")
                #        positive_edges_cor_dict, negative_edges_cor_dict , score_types = network.load_edges(positive_met_edges_cor_path, negative_met_edges_cor_path, negative_met_edges, score_types = ['Max','Avg','RAvg','RWA','RRWA'])
                #        performance_dict = network_performance.evaluate(positive_edges_cor_dict, negative_edges_cor_dict, positive_met_edges , negative_met_edges, score_types)
                #        read_write.to_pickle(performance_dict,  os.path.join(k_sub_outdir, "performance_dict.pkl"))
                #        with open(os.path.join(k_sub_outdir, "performance_summary.tsv"), "w") as f:
                #                f.write(f"score_type\tAVG(AUC_ROC,AUC_PRC)\tAUC_ROC\tAUC_PRC\n")
                #                for score_type in score_types:
                #                        f.write(f"{score_type}\t{performance_dict[score_type]['Quartiles']['AVG'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_ROC'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_PRC'][1]}\n")
                #        print(f"k= {k} completed.")


                      


                        



        



