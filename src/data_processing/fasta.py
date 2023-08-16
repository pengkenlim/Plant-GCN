#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

#func definitions

def load_fasta(path, replace, replace_with, exclude, sep=".",  level =1):
    """use sep to specify delimiter to split TranscriptID into GeneID
    Default is '.' 
    e.g. GENE_1234.1 --> GENE_1234"""
    fasta_contents_dict ={}
    with open(path,"r" ) as f:
        contents = f.read()
    contents = contents.split(">")
    contents = [chunk for chunk in contents if chunk != ""]
    for chunk in contents:
        transcript_ID = chunk.split("\n")[0].split(" ")[0].split("\t")[0]
        if len(transcript_ID.split(sep)) > 1:
            Gene_ID = sep.join(transcript_ID.split(sep)[:-level])
        else:
            Gene_ID = transcript_ID
        Gene_ID = Gene_ID.upper().replace(replace, replace_with)
        length = len("".join(chunk.split("\n")[1:]))
        if exclude.upper() not in transcript_ID.upper():
            fasta_contents_dict[transcript_ID]= {"Gene_ID":Gene_ID, "Fasta_Chunk": chunk, "Length": length}
    return fasta_contents_dict

def retain_longest_isoforms(fasta_contents_dict):
    gene_dict ={}
    for transcript_ID , info in fasta_contents_dict.items():
        Gene_ID = info["Gene_ID"]
        if  info["Length"] > gene_dict.get(Gene_ID , {"Length" :0})["Length"]:
            gene_dict[Gene_ID] = info
    return gene_dict

def fasta_dump(gene_dict , path):
    with open(path, "w") as f:
        for key , info in gene_dict.items():
            f.write(">" + info["Fasta_Chunk"])

def get_Tid_2_Gid_dict(fasta_contents_dict):
    Tid_2_Gid_dict = {}
    for key , info in fasta_contents_dict.items():
        if key != "" and info["Gene_ID"] != "":
            Tid_2_Gid_dict[key] = info["Gene_ID"].upper()
    return Tid_2_Gid_dict