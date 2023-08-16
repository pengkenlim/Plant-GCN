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

        parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
        help = "Working directory containing required data of which to output data." )

        #parser.add_argument("-i", "--input_path", type= str, metavar="", required = True,
        #help = "Path to expression matrix.")

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )
        
        parser.add_argument("-pc", "--pc_no", type= int, metavar="", default = 100,
        help = "Number of Pincipal components to keep. Default = 100. It cannot be more than the number of samples and number of genes in in the expression matrix." )

        args=parser.parse_args()
        threads = args.threads
        
        #set threads and then import
        os.environ["MKL_NUM_THREADS"] = str(threads)
        os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        from sklearn_wrappers import kmeans , PCA
        from data_processing import read_write, expression_matrix

        k_start = args.k_start
        k_end = args.k_end
        step = args.step
        k_list = [i for i in range(k_start,k_end+1,step)]

        #check validity of input then create outputdir if not already exists
        ouput_dir = args.output_dir
        #input_path = args.input_path
        sub_outdir =  os.path.join(ouput_dir, "Partition_expression_data")

        read_write.establish_dir(sub_outdir , isdir= True)

        #load exp_mat
        delimiter = args.delimiter
        if delimiter == "t":
                delim = "\t"
        else:
                delim= ","
        expmat_path =  os.path.join(ouput_dir, "QC_expression_data", f"expression_matrix.{delimiter}sv")
        pc_no = args.pc_no
        print("Loading expression matrix as pandas DataFrame...")
        expmat_df = expression_matrix.load(expmat_path , expmatsep = delim, indexcolname="target_id" )
        
        print(f"Performing standardization of expression matrix and embedding samples in {pc_no} principal components...")
        pca_data , pc_variances = PCA.standardize_transform(expmat_df, n_pcs=pc_no)
        del expmat_df
        print(f"PCA complete. {pc_no} PCs representing {sum(pc_variances)}% of variance retained. Writing PCA data to output..." )
        read_write.to_pickle(pca_data , os.path.join(sub_outdir,"PCA_data.pkl"))
        
        print(f"Performing K-means clustering....")
        k_cluster_assignment_dict , silhouette_coefficients, centroids_dict = kmeans.iterate_over_krange(pca_data, k_list ,randomstate=42)

        print(f"Writing k-means clustering data to output folder...")
        read_write.to_pickle(k_cluster_assignment_dict, os.path.join(sub_outdir,"k_cluster_assignment_dict.pkl"))
        read_write.to_pickle(silhouette_coefficients, os.path.join(sub_outdir, "silhouette_coefficients.pkl"))
        
        print("ks selected based on silhouette coefficient peaks:")
        selected_k = kmeans.select_k_window(silhouette_coefficients, k_cluster_assignment_dict)
        out = [print(f"k={k}")for k in selected_k]
        read_write.to_pickle(selected_k, os.path.join(sub_outdir,"selected_k.pkl"))

        print("Partition_expression_data.py complete")
    




