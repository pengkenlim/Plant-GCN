#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from scipy import stats, integrate
import numpy as np
from sklearn import metrics
import pandas as pd
import math

def tie_break(scores):
    scores = [score - (offset/100000000000) for offset, score in enumerate(scores)]
    return scores

def calc_AUC_ROC(pos_score, neg_score, return_thresholds =False):
    scores = pos_score + neg_score
    labels = [1 for i in pos_score] + [0 for i in neg_score]    
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, pos_label=1)
    AUC_ROC = metrics.auc(fpr, tpr)
    if return_thresholds:
        return fpr, tpr, thresholds, AUC_ROC
    else:
        return AUC_ROC

def calc_AUC_PRC(pos_score, neg_score, return_thresholds =False):
    scores = tie_break(pos_score) + tie_break(neg_score)
    labels = [1 for i in pos_score] + [0 for i in neg_score]    
    precision, recall, thresholds = metrics.precision_recall_curve(labels, scores, pos_label=1)
    AUC_PRC = np.abs(integrate.trapz(y=precision, x=recall))
    if return_thresholds:
        return precision, recall, thresholds, AUC_PRC
    else:
        return AUC_PRC


def extract_edge_score(edges, score_dict, score_type):
    """extract scores for a particcular score type"""
    scores = []
    for edge in edges:
        scores.append(score_dict[edge][score_type])
    return scores

def nan2value(scores, n=-1):
    """convert nan values to something else"""
    new_scores = [-1 if np.isnan(i) else i for i in scores]
    return new_scores

def calc_quartiles(performance_scores):
    return_array = list(stats.scoreatpercentile(performance_scores, [25,50,75], interpolation_method = "lower"))
    return return_array

def evaluate(positive_edges_cor_dict, negative_edges_cor_dict, positive_edges , negative_edges, score_types):
    """Main function."""
    performance_dict = {}
    for score_type in score_types:
        performance_dict[score_type] = {"AUC_ROC":{},"AUC_PRC":{}, "AVG":{}}
        positive_scores = nan2value(extract_edge_score( positive_edges, positive_edges_cor_dict, score_type))
        for ds , negative_edges_cor_dict_values  in negative_edges_cor_dict.items():
            negative_scores = nan2value(extract_edge_score(negative_edges[ds], negative_edges_cor_dict_values, score_type))
            performance_dict[score_type]["AUC_ROC"][ds] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =False)
            performance_dict[score_type]["AUC_PRC"][ds] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =False)
            performance_dict[score_type]["AVG"][ds] = np.mean([performance_dict[score_type]["AUC_ROC"][ds], performance_dict[score_type]["AUC_PRC"][ds]])
        performance_dict[score_type]["Quartiles"] = {}
        performance_dict[score_type]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type]["AUC_ROC"].values()))
        performance_dict[score_type]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type]["AUC_PRC"].values()))
        performance_dict[score_type]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type]["AVG"].values()))
    return performance_dict


def evaluate_ara(positive_edges_cor_dict, negative_edges_cor_dict, positive_All_edges , negative_All_edges, positive_met_edges, negative_met_edges,positive_GO_edges, negative_GO_edges, positive_TF_edges, negative_TF_edges, score_types):
    """main function."""
    performance_dict={}
    for score_type in score_types:
        performance_dict[score_type] = {}
        for edge_dataset in ["All", "Met", "GO", "TF"]:
            performance_dict[score_type][edge_dataset] = {"AUC_ROC":{},"AUC_PRC":{}, "AVG":{}}
            
            if edge_dataset == "All":
                positive_edges , negative_edges = positive_All_edges, negative_All_edges
            elif edge_dataset == "Met":
                positive_edges , negative_edges = positive_met_edges, negative_met_edges
            elif edge_dataset == "GO":
                positive_edges , negative_edges = positive_GO_edges, negative_GO_edges
            elif edge_dataset == "TF":
                positive_edges , negative_edges = positive_TF_edges, negative_TF_edges
            
            positive_scores = nan2value(extract_edge_score( positive_edges, positive_edges_cor_dict, score_type))
            positive_scores = [np.round(i, 5) for i in positive_scores]
            #if edge_dataset == "TF" and score_type == "Max":
            #        print("writing...")
            #        with open("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/PCC_8k_Max/TF_edges_kmeans.tsv", "w") as f:
            #            for edge , score in zip(positive_edges, positive_scores):
            #                f.write(f"{edge}\t{score}\n")
            for ds , negative_edges_cor_dict_values  in negative_edges_cor_dict.items():
                negative_scores = nan2value(extract_edge_score(negative_edges[ds], negative_edges_cor_dict_values, score_type))
                negative_scores = [np.round(i, 5) for i in negative_scores]
                performance_dict[score_type][edge_dataset]["AUC_ROC"][ds] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =False)
                performance_dict[score_type][edge_dataset]["AUC_PRC"][ds] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =False)
                performance_dict[score_type][edge_dataset]["AVG"][ds] = np.mean([performance_dict[score_type][edge_dataset]["AUC_ROC"][ds], 
                                                                                 performance_dict[score_type][edge_dataset]["AUC_PRC"][ds]])
            performance_dict[score_type][edge_dataset]["Quartiles"] = {}
            performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AUC_ROC"].values()))
            performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AUC_PRC"].values()))
            performance_dict[score_type][edge_dataset]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AVG"].values()))
            
            med_idx = list(performance_dict[score_type][edge_dataset]["AUC_ROC"].values()).index(performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_ROC"][1])
            negative_scores = nan2value(extract_edge_score(negative_edges[med_idx], negative_edges_cor_dict[med_idx], score_type))
            performance_dict[score_type][edge_dataset]["Thresholds"]={}
            performance_dict[score_type][edge_dataset]["Thresholds"]["AUC_ROC"] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =True)
            
            #AUC_PRC
            med_idx = list(performance_dict[score_type][edge_dataset]["AUC_PRC"].values()).index(performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_PRC"][1])
            negative_scores = nan2value(extract_edge_score(negative_edges[med_idx], negative_edges_cor_dict[med_idx], score_type))
            performance_dict[score_type][edge_dataset]["Thresholds"]["AUC_PRC"] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =True)

        #harmonic mean of scores
        performance_dict[score_type]["HM"] = {"AUC_ROC":{},"AUC_PRC":{}, "AVG":{}}
        for ds in negative_edges_cor_dict.keys():
            performance_dict[score_type]["HM"]["AUC_ROC"][ds] = stats.hmean([performance_dict[score_type][edge_dataset]["AUC_ROC"][ds] for edge_dataset in ["Met", "GO", "TF"]])
            performance_dict[score_type]["HM"]["AUC_PRC"][ds] = stats.hmean([performance_dict[score_type][edge_dataset]["AUC_PRC"][ds] for edge_dataset in ["Met", "GO", "TF"]])
            performance_dict[score_type]["HM"]["AVG"][ds] = np.mean([performance_dict[score_type]["HM"]["AUC_ROC"][ds], 
                                                                     performance_dict[score_type]["HM"]["AUC_PRC"][ds]])
            
            performance_dict[score_type]["HM"]["Quartiles"] = {}
            performance_dict[score_type]["HM"]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AUC_ROC"].values()))
            performance_dict[score_type]["HM"]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AUC_PRC"].values()))
            performance_dict[score_type]["HM"]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AVG"].values()))
            
    return performance_dict



def summarize_and_out(k_sub_outdir, score_types, performance_dict):
        with open(os.path.join(k_sub_outdir, "performance_summary.tsv"), "w") as f:
                f.write(f"score_type\tAVG(AUC_ROC,AUC_PRC)\tAUC_ROC\tAUC_PRC\n")
                for score_type in score_types:
                    f.write(f"{score_type}\t{performance_dict[score_type]['Quartiles']['AVG'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_ROC'][1]}\t{performance_dict[score_type]['Quartiles']['AUC_PRC'][1]}\n")

def summarize_and_out_ara(k_sub_outdir, score_types, performance_dict):
        with open(os.path.join(k_sub_outdir, "performance_summary.tsv"), "w") as f:
                format = ["score_type" , "AVG(AUCROC_HM, AUCPRC_HM)", "AVG(AUCROC_All, AUCPRC_All)" , "AVG(AUCROC_Met, AUCPRC_Met)", 
                          "AVG(AUCROC_GO, AUCPRC_GO)", "AVG(AUCROC_TF, AUCPRC_TF)"]
                f.write("\t".join(format) + "\n")
                for score_type in score_types:    
                    scores = [score_type, performance_dict[score_type]["HM"]["Quartiles"]["AVG"][1],
                              performance_dict[score_type]["All"]["Quartiles"]["AVG"][1],
                              performance_dict[score_type]["Met"]["Quartiles"]["AVG"][1],
                              performance_dict[score_type]["GO"]["Quartiles"]["AVG"][1],
                              performance_dict[score_type]["TF"]["Quartiles"]["AVG"][1]]
                    f.write("\t".join(scores) + "\n")


def cat_k_to_df(performance_df, score_types, performance_dict, k, metric , full = True):
    try:
        # check whether performance_df is empty
        out = performance_df[0]
        del out
    except:
        index_col = []
        index2_col = []
        if full:
            out = [[index_col.append(score_type) for i in range(3)] for score_type in score_types]
            out = [[index2_col.append(i) for i in ["Q1","Med","Q3"]] for score_type in score_types]
        else:
            out = [[index_col.append(score_type) for i in range(1)] for score_type in score_types]
            out = [index2_col.append("Med") for score_type in score_types]
        performance_df.index = index_col
        performance_df["Stat"] = index2_col
    
    col_values = []
    for score_type in score_types:
        if full:
            for i in range(3):
                col_values.append(performance_dict[score_type]['Quartiles'][metric][i])
        else:
            i = 1
            col_values.append(performance_dict[score_type]['Quartiles'][metric][i])

    performance_df[k] = col_values

    return performance_df


def cat_k_to_df_ara(performance_df, score_types, performance_dict, k, edge_dataset, metric , full = True):
    try:
        # check whether performance_df is empty
        out = performance_df[0]
        del out
    except:
        index_col = []
        index2_col = []
        if full:
            out = [[index_col.append(score_type) for i in range(3)] for score_type in score_types]
            out = [[index2_col.append(i) for i in ["Q1","Med","Q3"]] for score_type in score_types]
        else:
            out = [[index_col.append(score_type) for i in range(1)] for score_type in score_types]
            out = [index2_col.append("Med") for score_type in score_types]
        performance_df.index = index_col
        performance_df["Stat"] = index2_col
    
    col_values = []
    for score_type in score_types:
        if full:
            for i in range(3):
                col_values.append(performance_dict[score_type][edge_dataset]['Quartiles'][metric][i])
        else:
            i = 1
            col_values.append(performance_dict[score_type][edge_dataset]['Quartiles'][metric][i])

    performance_df[k] = col_values

    return performance_df


def best_worst_k(performance_df, score_types):
    result_string = ["Score_type\tmax_score\tmax_k\tmin_score\tmin_k\n"]
    for score_type in score_types:
        scores = list(performance_df.loc[score_type])[1:]
        k_list = list(performance_df.columns)[1:]
        
        max_score = np.max(scores )
        max_k = k_list[scores.index( max_score)]

        min_score = np.min(scores)
        min_k = k_list[scores.index( min_score)]

        result_string.append(f"{score_type}\t{max_score}\t{max_k}\t{min_score}\t{min_k}\n")
    result_string = "".join(result_string)
    return result_string

##############################################
def calculate_edge_index(s, t, num_nodes):
    if s > t:
        s, t = t, s  # Ensure s < t
        return False, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1
    else:
        return True, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1

def extract_edge_score_FULL(edges, network_object, gene_dict, score_type):
    """extract scores for a particcular score type"""
    scores = []
    for edge in edges:
        source, target = edge.split("--")
        try:
            source_index = gene_dict[source]
            target_index = gene_dict[target]
            _, edge_index = calculate_edge_index(source_index, target_index, len(gene_dict))       
            scores.append(network_object[score_type][edge_index])
        except:
            scores.append(math.nan)
    return scores


def evaluate_ara_FULL(gene_dict, network_object, positive_All_edges , negative_All_edges, positive_met_edges, negative_met_edges,positive_GO_edges, negative_GO_edges, positive_TF_edges, negative_TF_edges, score_types):
    """main function. For evaluation of FULL networks"""
    performance_dict={}
    for score_type in score_types:
        performance_dict[score_type] = {}
        for edge_dataset in ["All", "Met", "GO", "TF"]:
            performance_dict[score_type][edge_dataset] = {"AUC_ROC":{},"AUC_PRC":{}, "AVG":{}}
            if edge_dataset == "All":
                positive_edges , negative_edges = positive_All_edges, negative_All_edges
            elif edge_dataset == "Met":
                positive_edges , negative_edges = positive_met_edges, negative_met_edges
            elif edge_dataset == "GO":
                positive_edges , negative_edges = positive_GO_edges, negative_GO_edges
            elif edge_dataset == "TF":
                positive_edges , negative_edges = positive_TF_edges, negative_TF_edges
            
           
            if score_type in ["MR", "HRR"]:
                positive_scores = nan2value(extract_edge_score_FULL(positive_edges, network_object, gene_dict, score_type), len(gene_dict))
                positive_scores =  np.array(positive_scores) * -1  #so that positive samples have cor > than negative samples
                positive_scores = list(positive_scores)
            else:
                 positive_scores = nan2value(extract_edge_score_FULL(positive_edges, network_object, gene_dict, score_type))
                 #if edge_dataset == "TF":
                 #   print("writing...")
                 #   with open("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/PCC_8k_Max/TF_edges.tsv", "w") as f:
                 #       for edge , score in zip(positive_edges, positive_scores):
                 #           f.write(f"{edge}\t{score}\n")
            for ds , neg_edges  in negative_edges.items():
                
                if score_type in ["MR", "HRR"]:
                    negative_scores = nan2value(extract_edge_score_FULL(neg_edges, network_object, gene_dict, score_type), len(gene_dict))
                    negative_scores = np.array(negative_scores) * -1 
                    negative_scores = list(negative_scores)
                else:
                    negative_scores = nan2value(extract_edge_score_FULL(neg_edges, network_object, gene_dict, score_type))
                performance_dict[score_type][edge_dataset]["AUC_ROC"][ds] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =False)
                performance_dict[score_type][edge_dataset]["AUC_PRC"][ds] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =False)
                performance_dict[score_type][edge_dataset]["AVG"][ds] = np.mean([performance_dict[score_type][edge_dataset]["AUC_ROC"][ds], 
                                                                                 performance_dict[score_type][edge_dataset]["AUC_PRC"][ds]])
            performance_dict[score_type][edge_dataset]["Quartiles"] = {}
            performance_dict[score_type][edge_dataset]["Thresholds"] = {}
            #AUC_ROC
            performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AUC_ROC"].values()))            
            med_idx = list(performance_dict[score_type][edge_dataset]["AUC_ROC"].values()).index(performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_ROC"][1])
            if score_type in ["MR", "HRR"]:
                negative_scores = nan2value(extract_edge_score_FULL(negative_edges[med_idx], network_object, gene_dict, score_type), len(gene_dict))
                negative_scores = np.array(negative_scores) * -1 
                negative_scores = list(negative_scores)
            else:
                negative_scores = nan2value(extract_edge_score_FULL(negative_edges[med_idx], network_object, gene_dict, score_type))
            performance_dict[score_type][edge_dataset]["Thresholds"]["AUC_ROC"] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =True)
            #AUC_PRC
            performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AUC_PRC"].values()))
            med_idx = list(performance_dict[score_type][edge_dataset]["AUC_PRC"].values()).index(performance_dict[score_type][edge_dataset]["Quartiles"]["AUC_PRC"][1])
            if score_type in ["MR", "HRR"]:
                negative_scores = nan2value(extract_edge_score_FULL(negative_edges[med_idx], network_object, gene_dict, score_type), len(gene_dict))
                negative_scores = np.array(negative_scores) * -1 
                negative_scores = list(negative_scores)
            else:
                negative_scores = nan2value(extract_edge_score_FULL(negative_edges[med_idx], network_object, gene_dict, score_type))
            performance_dict[score_type][edge_dataset]["Thresholds"]["AUC_PRC"] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =True)
            #AVG
            performance_dict[score_type][edge_dataset]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type][edge_dataset]["AVG"].values()))
        #harmonic mean of scores
        performance_dict[score_type]["HM"] = {"AUC_ROC":{},"AUC_PRC":{}, "AVG":{}}
        for ds in negative_edges.keys():
            performance_dict[score_type]["HM"]["AUC_ROC"][ds] = stats.hmean([performance_dict[score_type][edge_dataset]["AUC_ROC"][ds] for edge_dataset in ["Met", "GO", "TF"]])
            performance_dict[score_type]["HM"]["AUC_PRC"][ds] = stats.hmean([performance_dict[score_type][edge_dataset]["AUC_PRC"][ds] for edge_dataset in ["Met", "GO", "TF"]])
            performance_dict[score_type]["HM"]["AVG"][ds] = np.mean([performance_dict[score_type]["HM"]["AUC_ROC"][ds], 
                                                                     performance_dict[score_type]["HM"]["AUC_PRC"][ds]])           
            performance_dict[score_type]["HM"]["Quartiles"] = {}
            performance_dict[score_type]["HM"]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AUC_ROC"].values()))
            performance_dict[score_type]["HM"]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AUC_PRC"].values()))
            performance_dict[score_type]["HM"]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type]["HM"]["AVG"].values()))          
    return performance_dict 
     
