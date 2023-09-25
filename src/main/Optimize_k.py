#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse
# import argparse to parse thread information so that we can set thread environment variable before importing numpy and sklearn modules

if __name__ == "__main__":
        
        parser= argparse.ArgumentParser(description="Optimize_k.py.\n\
                                        Evaluate the performance of ensemble GCN across selected range of k.\
                                        Only edges in our dataset will have their ensemble weight edges calculated.")
        
        parser.add_argument("-w", "--workers", type=int, metavar="", default=4,
        help = "Number of workers for parallelization." )

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Directory to output. Same as Label_edges.py and Partition_expression_data.py" )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )

        parser.add_argument("-cc","--correlation_coefficient", type=str, default = "bicor", choices = ["PCC","SCC","bicor"] ,
        help = "Correlation coefficient to optimize for. \'bicor\' by default.")

        parser.add_argument("-ara","--arabidopsis", action="store_true",
        help = "Performance evaluation for Arabidopsis. As seen in the publication.")

        parser.add_argument("-ara_dir","--arabidopsis_data_dir",type=str, default = "",
        help = "Directory for arabidopsis edges.")
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=12,
        help = "Number of threads for numpy operations." )

        args=parser.parse_args()
        
        workers=args.workers
        output_dir=args.output_dir
        delimiter=args.delimiter
        cc = args.correlation_coefficient
        arabidopsis=args.arabidopsis
        arabidopsis_data_dir = args.arabidopsis_data_dir
        threads = args.threads

        #set threads and then import
        os.environ["MKL_NUM_THREADS"] = str(threads)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
        os.environ["OMP_NUM_THREADS"] = str(threads)

        from coexpression import bicor , pearson , spearman #for now only bicor and PCC is online
        from data_processing import read_write, network
        from edge_labels import negative_edges
        from analyses import network_performance
        import numpy as np
        import pandas as pd

        
        sub_outdir = os.path.join(output_dir , "Optimize_k")
        if delimiter == "t":
                delim = "\t"
        else:
                delim = ","
        #load data needed from workdirectory
        Tid2Gid_dict = read_write.load_pickle(os.path.join(output_dir , "Remove_isoforms" , "Tid2Gid_dict.pkl"))
        k_cluster_assignment_dict =  read_write.load_pickle(os.path.join(output_dir, "Partition_expression_data", "k_cluster_assignment_dict.pkl"))
        selected_k =  read_write.load_pickle(os.path.join(output_dir, "Partition_expression_data", "selected_k.pkl"))
        expmat_path = os.path.join(output_dir, "QC_expression_data",  f"expression_matrix.{delimiter}sv")
        
        if arabidopsis:

                positive_All_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_All_edges.pkl"))
                negative_All_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_All_edges.pkl"))

                positive_met_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_met_edges.pkl"))
                negative_met_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_met_edges.pkl"))
                
                positive_GO_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_GO_edges.pkl"))
                negative_GO_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_GO_edges.pkl"))

                positive_TF_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_TF_edges.pkl"))
                negative_TF_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_TF_edges.pkl"))

                ##############################################################################################################

                #establish dataframes
                ##HM
                AVG_HM_full_df, AUC_ROC_HM_full_df, AUC_PRC_HM_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_HM_summ_df, AUC_ROC_HM_summ_df, AUC_PRC_HM_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()

                ##ALL
                AVG_All_full_df, AUC_ROC_All_full_df, AUC_PRC_All_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_All_summ_df, AUC_ROC_All_summ_df, AUC_PRC_All_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()

                ##met
                AVG_met_full_df, AUC_ROC_met_full_df, AUC_PRC_met_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_met_summ_df, AUC_ROC_met_summ_df, AUC_PRC_met_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()

                ##GO
                AVG_GO_full_df, AUC_ROC_GO_full_df, AUC_PRC_GO_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_GO_summ_df, AUC_ROC_GO_summ_df, AUC_PRC_GO_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()

                ##TF
                AVG_TF_full_df, AUC_ROC_TF_full_df, AUC_PRC_TF_full_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()
                AVG_TF_summ_df, AUC_ROC_TF_summ_df, AUC_PRC_TF_summ_df = pd.DataFrame() , pd.DataFrame(), pd.DataFrame()

                negative_All_edges_unpacked = []
                for edges in negative_All_edges.values():
                        negative_All_edges_unpacked.extend(edges)
                negative_All_edges_unpacked = list(set(negative_All_edges_unpacked))

                cc_sub_outdir = os.path.join(sub_outdir ,cc)
                
                print(f"Calculating ensemble scores of edges across k range for correlation coefficient: \'{cc}\'")
                k_cluster_assignment_dict[1]= np.array([ 0 for i in range(len(k_cluster_assignment_dict[selected_k[0]]))])
                selected_k = [1]+ selected_k

                for k in selected_k:
                        k_sub_outdir = os.path.join(cc_sub_outdir, f"{k}_K")
                        read_write.establish_dir(k_sub_outdir, isdir =True)
                        positive_All_edges_cor_path = os.path.join(k_sub_outdir, "positive_All_edges_cor.tsv")
                        negative_All_edges_cor_path = os.path.join(k_sub_outdir, "negative_All_edges_cor.tsv")
                        if cc == "bicor":
                                optimize_k = bicor.optimize_k
                        elif cc == "PCC":
                                optimize_k = pearson.optimize_k
                        elif cc == "SCC":
                                optimize_k = spearman.optimize_k
                        
                        #optimize_k(k, positive_All_edges_cor_path, negative_All_edges_cor_path, expmat_path, Tid2Gid_dict,  
                                        #k_cluster_assignment_dict, delim, workers, 
                                        #positive_All_edges, negative_All_edges_unpacked, threads)
                        score_types = ['Max','Avg','RAvg','RWA','RRWA']
                        
                        if k==1:
                                print("REMOVE THIS EDIT")
                                print(f"Loading network for k={k}")
                                positive_edges_cor_dict, negative_edges_cor_dict , score_types = network.load_edges(positive_All_edges_cor_path,
                                                                                                                negative_All_edges_cor_path, 
                                                                                                                negative_All_edges, 
                                                                                                        score_types = score_types)
                                print(f"Performance eval for k={k}")
                                performance_dict = network_performance.evaluate_ara( positive_edges_cor_dict,
                                                                negative_edges_cor_dict, 
                                                                positive_All_edges , 
                                                                negative_All_edges, 
                                                                positive_met_edges,
                                                                negative_met_edges,
                                                                positive_GO_edges,
                                                                negative_GO_edges,
                                                                positive_TF_edges,
                                                                negative_TF_edges,
                                                                score_types)
                                
                                read_write.to_pickle(performance_dict ,os.path.join(k_sub_outdir, "performance_dict_new.pkl"))
                                sys.exit()
                                
                                #performance_dict =  read_write.load_pickle(os.path.join(k_sub_outdir, "performance_dict_new.pkl"))
                                #HM
                                AVG_HM_full_df= network_performance.cat_k_to_df_ara(AVG_HM_full_df, score_types, performance_dict, k, "HM", "AVG" , full = True)
                                AUC_ROC_HM_full_df = network_performance.cat_k_to_df_ara(AUC_ROC_HM_full_df, score_types, performance_dict, k, "HM", "AUC_ROC" , full = True)
                                AUC_PRC_HM_full_df = network_performance.cat_k_to_df_ara(AUC_PRC_HM_full_df, score_types, performance_dict, k, "HM", "AUC_PRC" , full = True)

                                AVG_HM_summ_df= network_performance.cat_k_to_df_ara(AVG_HM_summ_df, score_types, performance_dict, k, "HM", "AVG" , full = False)
                                AUC_ROC_HM_summ_df = network_performance.cat_k_to_df_ara(AUC_ROC_HM_summ_df, score_types, performance_dict, k, "HM", "AUC_ROC" , full = False)
                                AUC_PRC_HM_summ_df = network_performance.cat_k_to_df_ara(AUC_PRC_HM_summ_df, score_types, performance_dict, k, "HM", "AUC_PRC" , full = False)

                                AVG_HM_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_HM_full_scores.csv"))
                                AUC_ROC_HM_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_HM_full_scores.csv"))
                                AUC_PRC_HM_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_HM_full_scores.csv"))
                                
                                AVG_HM_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_HM_summ_scores.csv"))
                                AUC_ROC_HM_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_HM_summ_scores.csv"))
                                AUC_PRC_HM_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_HM_summ_scores.csv"))


                                #All
                                AVG_All_full_df= network_performance.cat_k_to_df_ara(AVG_All_full_df, score_types, performance_dict, k, "All", "AVG" , full = True)
                                AUC_ROC_All_full_df = network_performance.cat_k_to_df_ara(AUC_ROC_All_full_df, score_types, performance_dict, k, "All", "AUC_ROC" , full = True)
                                AUC_PRC_All_full_df = network_performance.cat_k_to_df_ara(AUC_PRC_All_full_df, score_types, performance_dict, k, "All", "AUC_PRC" , full = True)

                                AVG_All_summ_df= network_performance.cat_k_to_df_ara(AVG_All_summ_df, score_types, performance_dict, k, "All", "AVG" , full = False)
                                AUC_ROC_All_summ_df = network_performance.cat_k_to_df_ara(AUC_ROC_All_summ_df, score_types, performance_dict, k, "All", "AUC_ROC" , full = False)
                                AUC_PRC_All_summ_df = network_performance.cat_k_to_df_ara(AUC_PRC_All_summ_df, score_types, performance_dict, k, "All", "AUC_PRC" , full = False)
                                
                                AVG_All_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_All_full_scores.csv"))
                                AUC_ROC_All_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_All_full_scores.csv"))
                                AUC_PRC_All_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_All_full_scores.csv"))
                                
                                AVG_All_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_All_summ_scores.csv"))
                                AUC_ROC_All_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_All_summ_scores.csv"))
                                AUC_PRC_All_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_All_summ_scores.csv"))

                                #Met
                                AVG_met_full_df= network_performance.cat_k_to_df_ara(AVG_met_full_df, score_types, performance_dict, k, "Met", "AVG" , full = True)
                                AUC_ROC_met_full_df = network_performance.cat_k_to_df_ara(AUC_ROC_met_full_df, score_types, performance_dict, k, "Met", "AUC_ROC" , full = True)
                                AUC_PRC_met_full_df = network_performance.cat_k_to_df_ara(AUC_PRC_met_full_df, score_types, performance_dict, k, "Met", "AUC_PRC" , full = True)

                                AVG_met_summ_df= network_performance.cat_k_to_df_ara(AVG_met_summ_df, score_types, performance_dict, k, "Met", "AVG" , full = False)
                                AUC_ROC_met_summ_df = network_performance.cat_k_to_df_ara(AUC_ROC_met_summ_df, score_types, performance_dict, k, "Met", "AUC_ROC" , full = False)
                                AUC_PRC_met_summ_df = network_performance.cat_k_to_df_ara(AUC_PRC_met_summ_df, score_types, performance_dict, k, "Met", "AUC_PRC" , full = False)

                                AVG_met_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_met_full_scores.csv"))
                                AUC_ROC_met_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_met_full_scores.csv"))
                                AUC_PRC_met_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_met_full_scores.csv"))
                                
                                AVG_met_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_met_summ_scores.csv"))
                                AUC_ROC_met_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_met_summ_scores.csv"))
                                AUC_PRC_met_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_met_summ_scores.csv"))
                                
                                #GO
                                AVG_GO_full_df= network_performance.cat_k_to_df_ara(AVG_GO_full_df, score_types, performance_dict, k, "GO", "AVG" , full = True)
                                AUC_ROC_GO_full_df = network_performance.cat_k_to_df_ara(AUC_ROC_GO_full_df, score_types, performance_dict, k, "GO", "AUC_ROC" , full = True)
                                AUC_PRC_GO_full_df = network_performance.cat_k_to_df_ara(AUC_PRC_GO_full_df, score_types, performance_dict, k, "GO", "AUC_PRC" , full = True)

                                AVG_GO_summ_df= network_performance.cat_k_to_df_ara(AVG_GO_summ_df, score_types, performance_dict, k, "GO", "AVG" , full = False)
                                AUC_ROC_GO_summ_df = network_performance.cat_k_to_df_ara(AUC_ROC_GO_summ_df, score_types, performance_dict, k, "GO", "AUC_ROC" , full = False)
                                AUC_PRC_GO_summ_df = network_performance.cat_k_to_df_ara(AUC_PRC_GO_summ_df, score_types, performance_dict, k, "GO", "AUC_PRC" , full = False)

                                
                                AVG_GO_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_GO_full_scores.csv"))
                                AUC_ROC_GO_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_GO_full_scores.csv"))
                                AUC_PRC_GO_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_GO_full_scores.csv"))
                                
                                AVG_GO_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_GO_summ_scores.csv"))
                                AUC_ROC_GO_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_GO_summ_scores.csv"))
                                AUC_PRC_GO_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_GO_summ_scores.csv"))

                                #TF
                                AVG_TF_full_df= network_performance.cat_k_to_df_ara(AVG_TF_full_df, score_types, performance_dict, k, "TF", "AVG" , full = True)
                                AUC_ROC_TF_full_df = network_performance.cat_k_to_df_ara(AUC_ROC_TF_full_df, score_types, performance_dict, k, "TF", "AUC_ROC" , full = True)
                                AUC_PRC_TF_full_df = network_performance.cat_k_to_df_ara(AUC_PRC_TF_full_df, score_types, performance_dict, k, "TF", "AUC_PRC" , full = True)

                                AVG_TF_summ_df= network_performance.cat_k_to_df_ara(AVG_TF_summ_df, score_types, performance_dict, k, "TF", "AVG" , full = False)
                                AUC_ROC_TF_summ_df = network_performance.cat_k_to_df_ara(AUC_ROC_TF_summ_df, score_types, performance_dict, k, "TF", "AUC_ROC" , full = False)
                                AUC_PRC_TF_summ_df = network_performance.cat_k_to_df_ara(AUC_PRC_TF_summ_df, score_types, performance_dict, k, "TF", "AUC_PRC" , full = False)

                                
                                AVG_TF_full_df.to_csv(os.path.join(cc_sub_outdir, "AVG_TF_full_scores.csv"))
                                AUC_ROC_TF_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_TF_full_scores.csv"))
                                AUC_PRC_TF_full_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_TF_full_scores.csv"))
                                
                                AVG_TF_summ_df.to_csv(os.path.join(cc_sub_outdir, "AVG_TF_summ_scores.csv"))
                                AUC_ROC_TF_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCROC_TF_summ_scores.csv"))
                                AUC_PRC_TF_summ_df.to_csv(os.path.join(cc_sub_outdir, "AUCPRC_TF_summ_scores.csv"))
                
                result_string = network_performance.best_worst_k(AVG_HM_summ_df, score_types)
                print(result_string)
                with open(os.path.join(cc_sub_outdir, "Best_worst_k.csv"), "w") as f:
                        f.write(result_string)
        else:
                positive_met_edges = read_write.load_pickle(os.path.join(output_dir, "Label_edges", "positive_met_edges.pkl"))
                negative_met_edges = read_write.load_pickle(os.path.join(output_dir, "Label_edges", "negative_met_edges.pkl"))

                negative_met_edges_unpacked = []
                for edges in negative_met_edges.values():
                        negative_met_edges_unpacked.extend(edges)
                negative_met_edges_unpacked = list(set(negative_met_edges_unpacked))
                
                
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
                                                                        positive_met_edges, 
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

