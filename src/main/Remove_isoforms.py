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
        
        parser.add_argument("-o", "--output_dir", type=str, metavar="", default="", 
        help = "Ouput working director to create subfolders in. Not required for --b / --batch_mode.")

        parser.add_argument("-de", "--delimiter", type=str, metavar="", default=".", 
        help = "Delimiter to use to find GeneID of transcripts. \".\" will be used by default.\n\
            For e.g. use --delimiter= \".\" if transcriptIDs look like  AT1GXXXXXX.1\n\
                Use --delimiter= \"_\" if transcriptIDs look like  AT1GXXXXXX_1")
        parser.add_argument("-l", "--level", type=int, metavar="", default=1, help= "Number of levels of delimiters to consider.\n\
                            Use level=1 if transcriptIDs look like  AT1GXXXXXX.1\n\
                            Use level=2 if transcriptIDs look like  AT1GXXXXXX.1.1.")
        
        parser.add_argument("-r", "--replace", type=str, metavar="", default="", help= "Substring within transcriptIDs to replace in order to get geneIDs")
        parser.add_argument("-rw", "--replace_with", type=str, metavar="", default= "", help= "replacement string")
        parser.add_argument("-ex", "--exclude", type=str, metavar="", default= "impossibleplaceholder", help= "exclude transcripts with these substrings")
        parser.add_argument("-des", "--description_mode", type=str, metavar="", default= "impossibleplaceholder", help= "use description to extract geneID. Start flank seq to split description ")
        parser.add_argument("-des_end", "--description_mode_end", type=str, metavar="", default= " ", help= "end flank seq to split description")
        
        args=parser.parse_args()
        if args.batch_mode == None:
                #proceed with single mode
                input_fasta = args.input_fasta
                output_dir = args.output_dir
                delimiter= args.delimiter
                level = args.level
                replace = args.replace
                replace_with = args.replace_with
                exclude = args.exclude
                description_mode = args.description_mode
                description_mode_end = args.description_mode_end
                if output_dir == "":
                    sys.exit("--output not specified")
                read_write.establish_dir(output_dir, isdir = True)
                sub_outdir = os.path.join(output_dir, "Remove_isoforms")
                read_write.establish_dir(sub_outdir, isdir = True)

                fasta_contents_dict = fasta.load_fasta(input_fasta, replace, replace_with, exclude, description_mode, description_mode_end, sep = delimiter, level = level)
                print(f"\n\nTotal number of transcripts in {input_fasta}: \n{len(fasta_contents_dict)}\n")
                gene_dict =  fasta.retain_longest_isoforms(fasta_contents_dict)
                print(f"Total number of genes: \n{len(gene_dict)}\n")
                fasta_outpath = os.path.join(sub_outdir, "primary_transcripts.fasta")
                fasta.fasta_dump(gene_dict, fasta_outpath)
                print(f"Primary transcripts fasta written to: \n{fasta_outpath}\n")
                Tid2Gid_dict = fasta.get_Tid_2_Gid_dict(fasta_contents_dict)
                eg_tid = list(Tid2Gid_dict.keys())[1]
                eg_gid = Tid2Gid_dict[eg_tid]
                print(f"Example of how Gene_ID is obtained from Transcript_ID:\n\
                      {eg_tid} --> {eg_gid}\n\n")

                read_write.to_pickle( Tid2Gid_dict, os.path.join(sub_outdir , "Tid2Gid_dict.pkl"))
        else:
            #batch mode is not yet implemented
            pass
        
        print("Remove_isoforms.py completed")