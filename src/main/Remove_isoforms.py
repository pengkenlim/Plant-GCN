import os
import sys
import argparse
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)

from data_processing import fasta, read_write

if __name__ == "__main__":
        parser= argparse.ArgumentParser(description="Remove_isoforms.py\n\
                                        Parse cDNA fasta file(s). Retain only longest isoform for each gene.")
        ME_group_1 = parser.add_mutually_exclusive_group(required=True)
        ME_group_1.add_argument("-b", "--batch_mode", type=str, metavar="", help= "Process multiple fasta files in one run. Provide path to TSV file containing path to \
                                fasta files, output paths and parameters.")
        ME_group_1.add_argument("-i", "--input_fasta", type=str, metavar="", help= "Input path to fasta file to process")
        
        parser.add_argument("-o", "--output", type=str, metavar="", default="", 
        help = "Output path to write processed cDNA fasta file. Not required for --b / --batch_mode.")

        parser.add_argument("-de", "--delimiter", type=str, metavar="", default=".", 
        help = "Delimiter to use to find GeneID of transcripts. \".\" will be used by default.\n\
            For e.g. use --delimiter= \".\" if transcriptIDs look like  AT1GXXXXXX.1\n\
                Use --delimiter= \"_\" if transcriptIDs look like  AT1GXXXXXX_1")
        parser.add_argument("-l", "--level", type=int, metavar="", default=1, help= "Number of levels of delimiters to consider.\n\
                            Use level=1 if transcriptIDs look like  AT1GXXXXXX.1\n\
                            Use level=2 if transcriptIDs look like  AT1GXXXXXX.1.1.")
        
        args=parser.parse_args()
        if args.batch_mode == None:
                #proceed with single mode
                input_fasta = args.input_fasta
                output = args.output
                delimiter= args.delimiter
                level = args.level
                if output == "":
                    sys.exit("--output not specified")
                if not os.path.exists("/".join(output.split("/")[:-1])):
                       os.makedirs("/".join(output.split("/")[:-1]))

                fasta_contents_dict = fasta.load_fasta(input_fasta, sep = delimiter, level = level)
                print(f"Total number of transcripts in {input_fasta}: \n{len(fasta_contents_dict)-1}")
                gene_dict =  fasta.retain_longest_isoforms(fasta_contents_dict)
                print(f"Total number of genes: \n{len(gene_dict)}")
                fasta.fasta_dump(gene_dict,output)
                print(f"Primary transcripts fasta written to: \n{output}")
                Tid2Gid_dict = fasta.get_Tid_2_Gid_dict(gene_dict)
                read_write.to_pickle( Tid2Gid_dict, ".".join(output.split(".")[:-1])+"Tid2Gid_dict.pkl")
        else:
            pass
        
"/mnt/g/My\ Drive/Ken/projects/Kingdom/Kingdom_scripts/common_data/Annotations/ATTED-II_group/3711"