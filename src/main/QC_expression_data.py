#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

import argparse
from data_processing import expression_matrix , read_write



if __name__ == "__main__":
        parser= argparse.ArgumentParser(description="QC_expression_data.py\n\
                                        Load expression matrix. Perform QC via filtering out undesirable samples based on P_pseudoaligned mapping statistics from kallisto.")
      
        parser.add_argument("-ps", "--p_pseudoaligned", type=float, metavar="", default=20, 
        help = "Specifies p_pseudoaligned as outputed by kallisto as threshold for quality control. Expression data of accessions that do not meet this threshold will be removed from expression data.")

        parser.add_argument("-o", "--ouput_dir", type=str, metavar="", required = True,
        help = "Working directory containing required data of which to output data." )

        parser.add_argument("-i", "--input_path", type= str, metavar="", required = True,
        help = "Path to expression matrix.")

        parser.add_argument("-id", "--input_data_path", type= str, metavar="", required = True,
        help = "Path to kallisto pseudoalignment data. TSV file where each row is \'<SAMPLE>tabs<P_PSEUDOALIGNED>tabs<P_PSEUDOALIGNED>\'.")

        parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
        help = "Delimiter for input expression matrix.  -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default." )
        
        parser.add_argument("-t_n", "--tranpose_needed", action="store_true",
        help = "Use this if input expression matrix needs transposing (i.e. cols are genes, rows corresponsds to samples.")

        args=parser.parse_args()
        p_pseudoaligned= args.p_pseudoaligned
        ouput_dir= args.ouput_dir
        input_path= args.input_path
        input_data_path = args.input_data_path
        delimiter= args.delimiter
        tranpose_needed= args.tranpose_needed
        
        if delimiter == "t":
            delim = "\t"
        else:
            delim= ","

        #reading matrix as dataframe. Transpose if needed
        maprate_dict = expression_matrix.load_qc(input_data_path)
        if tranpose_needed:
            print("Reading expression matrix with transponse enabled...\n")
            expmat_df = expression_matrix.load_transposed(input_path, expmatsep = delim, indexcolname="accession")
        else:
            print("Reading expression matrix...\n")
            expmat_df = expression_matrix.load(input_path, expmatsep = delim, indexcolname="target_id")
        
        total, failed, passed, cutoff = expression_matrix.thresholder(maprate_dict, p_pseudoaligned)
        print(f"\n\nQC completed. Here are the stats:\nTotal samples: {len(total)}\nFailed: {len(failed)}\nPassed: {len(passed)}\nP_ps threshold: {cutoff}\n\n")
        
        sub_outdir =  os.path.join(ouput_dir, "QC_expression_data")
        #writing summarry.
        summary_path = os.path.join(sub_outdir, "qc_summary.txt")
        read_write.establish_dir(summary_path)
        with open(summary_path, "w") as f:
              f.write(f"QC completed. Here are the stats:\nTotal samples: {len(total)}\nFailed: {len(failed)}\nPassed: {len(passed)}\nP_ps threshold: {cutoff}")
        
        #removing failed samples form expression matrix and writting to file.
        print("Removing failed samples form expression matrix and writting to output...\n")
        #print(expmat_df)
        expmat_df = expression_matrix.subset(expmat_df, passed)
        expression_matrix.write(expmat_df, os.path.join(sub_outdir, f"expression_matrix.{delimiter}sv"), expmatsep = delim)

        print("QC_expression_data.py completed\n")


