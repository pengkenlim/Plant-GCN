#LABEL EDGES
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc POPLAR -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3694/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc CHINESECABBAGE -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3711/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc ORYZA -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid39947/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 1 -sc ARA -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/

#QC EXP datapython 
#python ~/Plant-GCN/src/main/QC_expression_data.py -t_n --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-out/taxid3702_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-matrices/taxid3702_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702/
#python ~/Plant-GCN/src/main/QC_expression_data.py -t_n --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-out/taxid3694_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-matrices/taxid3694_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3694/
#python ~/Plant-GCN/src/main/QC_expression_data.py -t_n --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-out/taxid3711_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-matrices/taxid3711_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3711/
#python ~/Plant-GCN/src/main/QC_expression_data.py -t_n --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-out/taxid39947_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b2/postprocess/qc-matrices/taxid39947_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid39947/

#Partition expression data
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 1000 -pc 1000 -st 2 -ws 19 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702 -np -nk
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 384 -pc 384 -st 1 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3694
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 406 -pc 406 -st 1 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3711
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 991 -pc 991 -st 2 -ws 19 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid39947
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 100 -pc 100 -st 1 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 250 -pc 250 -st 1 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5k
python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 100 -pc 100 -st 1 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h
