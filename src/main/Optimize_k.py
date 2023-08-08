#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse

if __name__ == "__main__":
        parser= argparse.ArgumentParser(description="Optimize_k.py\n\
                                        ")
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=4,
        help = "Number of threads to use for PCA and Kmeans clustering. Handled by sklearn's parrelization." )

        parser.add_argument("-o", "--ouput_dir", type=str, metavar="", required = True,
        help = "Directory to output. Same as Label_edges.py and Partition_expression_data.py" )

        parser.add_argument("-i", "--input_path", type= str, metavar="", required = True,
        help = "Path to expression matrix.")

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )
        
        parser.add_argument("-pc", "--pc_no", type= int, metavar="", default = 100,
        help = "Number of Pincipal components to keep. Default = 100. It cannot be more than the number of samples and number of genes in in the expression matrix." )
        args=parser.parse_args()