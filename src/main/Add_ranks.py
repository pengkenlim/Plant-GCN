#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse
import numpy as np
import scipy

from data_processing import read_write , network

if __name__ == "__main__":
    parser= argparse.ArgumentParser(description="Add_ranks.py.\n\
                                            calculate Mutual Ranks (MR) and Highest reciprocal ranks (HRR) from raw correlation scores.")
    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
                        help = "Working directory containing required data of which to output data." )
    parser.add_argument("-n", "--net_name" , type=str, metavar = "" , default = "All,",
                        help = "name of network to calculate MR and HRR for. By default, all networks will be calculated")
    
    args=parser.parse_args()
    output_dir = args.output_dir
    net_name = args.net_name
    BGCN_sub_outdir = os.path.join(output_dir , "Build_ensemble_GCN")
    sub_outdir = os.path.join(output_dir , "Add_ranks")

    All_networks_paths =[]
    for cc in ["PCC", "SCC", "bicor" ]:
           netdirs = os.listdir(os.path.join(BGCN_sub_outdir, cc))
           for netdir in netdirs:
                if "k" in netdir:
                    All_networks_paths.append(os.path.join( BGCN_sub_outdir, cc, netdir))
    
    if net_name == "All":
        print(f"net_name={net_name} specified by user.\nProceeding to calculate MR and HRR for these networks:")
        out = [print(network_path) for network_path in All_networks_paths]
    else:
        print(f"net_name={net_name} specified by user.\nProceeding to calculate MR and HRR for:")
        All_networks_paths = [network_path for network_path in All_networks_paths if net_name.lower() in network_path.lower()]
        try:
            print(All_networks_paths[0])
        except:
            sys.exit(f"{net_name} network not found!")
    
    #start main loop
    for network_path in All_networks_paths:
        print(f"Loading network: {network_path}\n into adjacency matrix")
        genes , upper_half, lower_half, gene_dict = network.extract_ranks(network_path)
        
        print("Calculating HRR and MR of all edges...")
        HRR_array = np.max([upper_half , lower_half], axis = 0)
        print("HRR calculated")
        MR_array = np.sqrt(upper_half * lower_half)
        print("MR calculated")
        
        del upper_half, lower_half

        print("Writing ranks to new network at:")
        network_outpath = os.path.join( sub_outdir , network_path.split("/")[-1])
        read_write.establish_dir(network_outpath, isdir = True)
        print(network_outpath)
        network.write_MR_HRR(gene_dict, HRR_array, MR_array, network_path, network_outpath, genes)

        print("Add_ranks.py completed")
        
        


