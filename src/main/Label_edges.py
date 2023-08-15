#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from edge_labels import negative_edges, plantcyc
from data_processing import read_write

import argparse


if __name__ == "__main__":
        parser= argparse.ArgumentParser(description="Label_edges.py\n\
                                        Mine metabolic edges from Plant-cyc and construct sample-balanced labeled datasets for performance evaluation.\
                                        Metabolic edges are labeled positive samples while negative edges are randomly chosen.\
                                        Edges connnect enzymes within the same pathway that catalyses different reaction step.")


        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Working directory containing required data of which to output data. Requires Tid2Gid_dict.pkl from Remove_isoforms.py to get GeneIDs." )

        parser.add_argument("-sc", "--species_code", type= str, metavar="", required = True,
        help = "species code for querying plantcyc.")
        
        parser.add_argument("-cri", "--criterion", type = int, default = 2,  choices = [1,2],
        help = "Method of filtering edges. 1: Only retains edges connecting enzymes with experimentally verified reactions (EV code: EXP).\
                2: No EV requirements. RXNIDs with more than 5 enzymes assigned are considered as ambiguous and ignored.\
                Use 1 for model plants such as A. thaliana / Z. mays. Use 2 for non-model organisms.")
        
        parser.add_argument("-ds", "--dataset_size", type = int, default = 100, 
        help = "Number of datasets to generate. Increase more for better estimation of GCN performance and performance margin further down the line.")
        
        #arg parse
        args=parser.parse_args()
        output_dir = args.output_dir
        species_code = args.species_code
        criterion = args.criterion
        dataset_size = args.dataset_size

        
        #read Tid2Gid_dict.pkl to get all genes
        Tid2Gid_dict_path = os.path.join(output_dir, "Remove_isoforms", "Tid2Gid_dict.pkl")
        Tid2Gid_dict = read_write.load_pickle(Tid2Gid_dict_path)
        All_genes = list(Tid2Gid_dict.values())


        #mine Plantcyc using species code
        print(f"Getting pathways from plantcyc using species_code = {species_code} ...")
        met_annot_dict =  plantcyc.get_pathways(species_code)
        print(f"Getting pathways info and generating edges... Criterion = {criterion} ")
        if criterion == 1:
                met_annot_dict = plantcyc.mine_info_generate_edges(species_code, met_annot_dict, All_genes, EXP=True)
        elif criterion == 2:
                met_annot_dict = plantcyc.mine_info_generate_edges(species_code, met_annot_dict, All_genes, criterion_2_cutoff=5)

        print("Mining completed. Saving positive sample data to output dir...")
        #make subdirectory to write out
        sub_outdir = os.path.join(output_dir, "Label_edges")
        read_write.establish_dir(sub_outdir, isdir=True)
        
        #save tsv and pickle file
        read_write.to_pickle(met_annot_dict, os.path.join(sub_outdir, "met_annot_dict.pkl"))
        
        plantcyc.edge_dump(met_annot_dict, os.path.join(sub_outdir, "positive_met_edges.tsv") , type = f"Cri_{criterion}")
        
        positive_met_edges = plantcyc.extract_edges(met_annot_dict, type = f"Cri_{criterion}")
        read_write.to_pickle(positive_met_edges, os.path.join(sub_outdir, "positive_met_edges.pkl"))
        print(f"{len(positive_met_edges)} positive edges generated using Criterion = {criterion}.\n")
        
        print("Generating negative sample data...")
        negative_met_edges = negative_edges.generate_dict(positive_met_edges, All_genes ,iterations = dataset_size)
        read_write.to_pickle(negative_met_edges, os.path.join(sub_outdir, "negative_met_edges.pkl"))
        
        print("Label_edges.py completed")

        



