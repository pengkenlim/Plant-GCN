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
import math

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
    scores = pos_score + neg_score
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
    new_scores = [-1 if np.isnan(i) else -i for i in scores]
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
        negative_scores = nan2value(extract_edge_score(negative_edges, negative_edges_cor_dict_values, score_type))
        
        performance_dict[score_type]["AUC_ROC"][ds] = calc_AUC_ROC(positive_scores, negative_scores, return_thresholds =False)
        performance_dict[score_type]["AUC_PRC"][ds] = calc_AUC_PRC(positive_scores, negative_scores, return_thresholds =False)
        performance_dict[score_type]["AVG"][ds] = np.mean(performance_dict[score_type]["AUC_ROC"][ds], performance_dict[score_type]["AUC_PRC"][ds])
    performance_dict[score_type]["Quartiles"] = {}
    performance_dict[score_type]["Quartiles"]["AUC_ROC"] = calc_quartiles(list(performance_dict[score_type]["AUC_ROC"].values()))
    performance_dict[score_type]["Quartiles"]["AUC_PRC"] = calc_quartiles(list(performance_dict[score_type]["AUC_PRC"].values()))
    performance_dict[score_type]["Quartiles"]["AVG"] = calc_quartiles(list(performance_dict[score_type]["AVG"].values()))
    return performance_dict
