#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

import argparse
# import argparse to parse thread information so that we can set thread environment variable before importing numpy and sklearn modules
if __name__ == "__main__":
        parser= argparse.ArgumentParser(description="Partition_expression_data.py\n\
                                        Load expression matrix. Perform sample standardization to correct batch effecct. Embed samples from gene-space into PC-space. Perform kmeans clustering across a range of K.\n\
                                        Outputs principal component embedings of samples (PCA data), mean silhoette coefficients from each clustering iteration, kmeans cluster assignments of samples in pickle format.")

        parser.add_argument("-ks", "--k_start", type=int, metavar="", default=2, 
        help = "Starting k for kmeans clustering.")
        
        parser.add_argument("-ke", "--k_end", type=int, metavar="", default=10, 
        help = "Ending k for kmeans clustering.")

        parser.add_argument("-st", "--step", type=int, metavar="", default=1,
        help = "Step to increase between each k. User is suggested to increase step if range is very big to decrease search space." )
        
        parser.add_argument("-t", "--threads", type=int, metavar="", default=4,
        help = "Number of threads to use for PCA and Kmeans clustering. Handled by sklearn's parrelization." )

        parser.add_argument("-o", "--ouput_dir", type=str, metavar="", required = True,
        help = "Directory to output." )

        parser.add_argument("-i", "--input_path", type= str, metavar="", required = True,
        help = "Path to expression matrix." )

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )

        args=parser.parse_args()
        threads = args.threads
        
        #set threads and then import
        os.environ["MKL_NUM_THREADS"] = str(threads)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        from sklearn_wrappers import kmeans , PCA

        k_start = args.k_start
        k_end = args.k_end
        step = args.step
        k_list = [i for i in range(k_start,k_end,step)]

        #check validity of input and output_dir
        




