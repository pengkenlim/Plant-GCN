#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from data_processing import read_write , network
from analyses import network_performance
import argparse

if __name__ == "__main__":
    parser= argparse.ArgumentParser(description="Eval_ara_full_net_performance.py.py.\n\
                                            Evaluate predicitive performances of Arabidopsis GCNs produced by the pipline.")
    
    
    parser.add_argument("-n", "--net_name" , type=str, metavar = "" , default = "All,",
                        help = "name of network to evaluate. By default, all networks will be evaluated")
    parser.add_argument("-o", "--output_dir", type=str, required = True, help = "Working directory containing required data of which to output data." )
    
    
    parser.add_argument("-ara_dir","--arabidopsis_data_dir",type=str, default = "",
    help = "Directory for arabidopsis edges.")

    parser.add_argument("-ATTED","--ATTED_PATH",type=str, default = "",
    help = "Path to ATTED network directory.")

    parser.add_argument("-map_path","--mapping_path",type=str, default = "",
    help = "Path to pickled ATTED to AGi mapping dictionary.")
    
    args=parser.parse_args()
    output_dir = args.output_dir
    net_name = args.net_name
    arabidopsis_data_dir = args.arabidopsis_data_dir
    #arabidopsis_data_dir = "/home/ken/Plant-GCN/test_data/Arabidopsis_edges/" ##remove
    ATTED_PATH = args.ATTED_PATH
    mapping_path = args.mapping_path
    AR_sub_outdir = os.path.join(output_dir , "Add_ranks")
    sub_outdir = os.path.join(output_dir , "Evaluate_full_network")

    All_networks_paths =[]
    for netdir in os.listdir(AR_sub_outdir):
            if "k" in netdir:
                All_networks_paths.append(os.path.join( AR_sub_outdir, netdir))
    
    if ATTED_PATH =="":
    
        if net_name == "All":
            print(f"net_name={net_name} specified by user.\nProceeding to evaluate these networks:")
            out = [print(network_path) for network_path in All_networks_paths]
        else:
            print(f"net_name={net_name} specified by user.\nProceeding to evaluate:")
            All_networks_paths = [network_path for network_path in All_networks_paths if net_name.lower() in network_path.lower()]
            try:
                print(All_networks_paths[0])
            except:
                sys.exit(f"{net_name} network not found!")
    else:
        if os.path.exists(ATTED_PATH):
            All_networks_paths=[ATTED_PATH]
            ATTED_convert_dict = read_write.load_pickle(mapping_path)
            #print(ATTED_convert_dict)
            #print(len(ATTED_convert_dict))
            #sys.exit()
        else:
             sys.exit(f"{ATTED_PATH} network not found!")
    
    #loading edges
    print("Loading edge datasets...")
    positive_All_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_All_edges.pkl"))
    negative_All_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_All_edges.pkl"))

    positive_met_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_met_edges.pkl"))
    negative_met_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_met_edges.pkl"))
                
    positive_GO_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_GO_edges.pkl"))
    negative_GO_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_GO_edges.pkl"))

    positive_TF_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir, "positive_TF_edges.pkl"))
    negative_TF_edges = read_write.load_pickle(os.path.join(arabidopsis_data_dir , "negative_TF_edges.pkl"))
    
    #start main loop
    for network_path in All_networks_paths:
        sub_sub_outdir = os.path.join(sub_outdir, network_path.split("/")[-1])
        read_write.establish_dir(sub_sub_outdir, isdir = True)
        print(f"Loading {network_path}...")
        
        if ATTED_PATH =="":
            genes , network_object, gene_dict = network.load_full_network(network_path)
        else:
            genes, network_object, gene_dict = network.load_full_network_ATTED(network_path, ATTED_convert_dict)
            
        print(genes[:10])
        print(network_object["RAW"][:10])
        print( "len(gene_dict):", len(gene_dict))

        read_write.to_pickle(genes , os.path.join(sub_sub_outdir, "genes.pkl") )
        read_write.to_pickle(network_object ,os.path.join(sub_sub_outdir, "network_object.pkl") )
        read_write.to_pickle(gene_dict ,os.path.join(sub_sub_outdir, "gene_dict.pkl") )
        #genes = read_write.load_pickle( os.path.join(sub_sub_outdir, "genes.pkl"))
        #network_object = read_write.load_pickle("/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/Evaluate_full_network/networkobject.pkl")
        #network_object = read_write.load_pickle( os.path.join(sub_sub_outdir, "network_object.pkl"))
        #gene_dict = read_write.load_pickle( os.path.join(sub_sub_outdir, "gene_dict.pkl"))

        score_types = list(network_object.keys())
        print(f"Evaluating performance...")
        performance_dict = network_performance.evaluate_ara_FULL(gene_dict, network_object, positive_All_edges , negative_All_edges, positive_met_edges, negative_met_edges, positive_GO_edges, negative_GO_edges, positive_TF_edges, negative_TF_edges, score_types)
        read_write.to_pickle(performance_dict ,os.path.join(sub_sub_outdir, "performance_dict.pkl") )
        