#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import requests
from collections import Counter

def get_pathways(PMNCODE):
    met_annot_dict = {}
    exp_pathways_URL=f"https://pmn.plantcyc.org/{PMNCODE}/search-query?type=PATHWAY"
    page_string = requests.get(exp_pathways_URL).text
    
    pathway_codes = page_string.split("type=PATHWAY&object=")[1:]
    pathway_codes = [chunk.split("</A")[0] for chunk in pathway_codes]
    for PWY_NAME in pathway_codes:
        met_annot_dict[PWY_NAME.split("\">")[0]] = {"NAME" : PWY_NAME.split("\">")[1]}

    return met_annot_dict


def mine_info_generate_edges(PMNCODE, met_annot_dict, All_genes, EXP=False, criterion_2_cutoff = 1):
    for P_idx, PWY in enumerate(met_annot_dict.keys(), start = 1):
        URL=f"https://pmn.plantcyc.org/{PMNCODE}/pathway-genes?object=" + PWY
        page_string = requests.get(URL).text.split("\n")
        NAME= page_string[0]
        Genes={}
        Edges={"All":[], "Cri_1": [], "Cri_2": []}
        for line in page_string[3:]:
            if line != "":
                geneID = line.split("\t")[1].split(".")[0].upper()
                if EXP:
                    if line.split("\t")[-1] == "EV-EXP" and geneID in All_genes:
                    #if line.split("\t")[-1] == "EV-EXP":                
                        Genes[geneID]= line.split("\t")[3]
                else:
                    if geneID in All_genes:
                    #if True:
                        Genes[geneID]= line.split("\t")[3]
        met_annot_dict[PWY]["Genes"] = Genes

        #for cri_2
        RXN_counts = Counter(list(Genes.values()))

        for idx , source in enumerate(list(Genes.keys())):
            S_RXN = Genes[source]
            for target in list(Genes.keys())[idx+1:]:
                if target != source:
                    Edges["All"].append( "-".join(sorted([source, target])) )
                    if Genes[target] != S_RXN:
                        Edges["Cri_1"].append( "-".join(sorted([source, target])) )
                        if RXN_counts[Genes[target]] <= criterion_2_cutoff and RXN_counts[S_RXN] <= criterion_2_cutoff:
                            Edges["Cri_2"].append( "-".join(sorted([source, target])) )            

        met_annot_dict[PWY]["Edges"] = Edges

        met_annot_dict[PWY]["Name"] = NAME

        if P_idx % 10 == 0:
            print(P_idx, "pathways mined")
    return met_annot_dict

def edge_dump(met_annot_dict, path ,type="Cri_1"):
    string = []
    for PWY , info in met_annot_dict.items():
        Name = info["Name"]
        for edge in info["Edges"][type]:
            string.append(f"{edge}\t{PWY}\t{Name}\n")
    string = "".join(string)
    with open(path, "w") as f:
        f.write(string)

def extract_edges(met_annot_dict,type="Cri_1"):
    positive_met_edges = []
    for PWY , info in met_annot_dict.items():
        for edge in info["Edges"][type]:
            positive_met_edges.append(edge)
    return positive_met_edges


if __name__ == "__main__": # testing the code
    PMNCODE = "MTRUNCATULA"
    All_genes = []
    met_annot_dict_2 = get_pathways(PMNCODE)
    met_annot_dict_2 = mine_info_generate_edges(PMNCODE, met_annot_dict_2, All_genes, EXP=False, criterion_2_cutoff = 5)
    path = "./met_edges.tsv"
    edge_dump(met_annot_dict_2, path, type = "Cri_2")