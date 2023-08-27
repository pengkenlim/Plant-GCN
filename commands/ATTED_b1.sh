#LABEL EDGES

#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc GRAPE -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid29760/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc CORN -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4577/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc TOMATO -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4018/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc MMTRUNCATULA -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3880/
#python ~/Plant-GCN/src/main/Label_edges.py -ds 100 -cri 2 -sc SOY -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3847/

#QC-EXPRESSION DATA
#python ~/Plant-GCN/src/main/QC_expression_data.py  --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-out/taxid29760_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-matrices/taxid29760_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid29760/ -t_n
#python ~/Plant-GCN/src/main/QC_expression_data.py  --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-out/taxid4577_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-matrices/taxid4577_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4577 -t_n
#python ~/Plant-GCN/src/main/QC_expression_data.py  --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-out/taxid4081_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-matrices/taxid4081_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4081 -t_n
#python ~/Plant-GCN/src/main/QC_expression_data.py  --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-out/taxid3880_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-matrices/taxid3880_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3880 -t_n
#python ~/Plant-GCN/src/main/QC_expression_data.py  --input_data_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-out/taxid3847_qc_info.tsv --input_path /mnt/md2/ken/correlation_networks/Plant-GCN_data/From_wallace/ATTED_b1/postprocess/qc-matrices/taxid3847_expmat.tsv --ouput_dir /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3847 -t_n

#Partition expression data
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 100 -ks 2 -ke 202 -pc 202 -ws 10 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3880
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 100 -ks 2 -ke 601 -pc 601 -ws 11 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid29760
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 100 -ks 2 -ke 660 -pc 660 -ws 13 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid3847
#python ~/Plant-GCN/src/main/Partition_expression_data.py -t 100 -ks 2 -ke 806 -pc 806 -st 2 -ws 16 -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4081
python ~/Plant-GCN/src/main/Partition_expression_data.py -t 64 -ks 2 -ke 1000 -pc 1000 -st 5  -ws 19 -np -nk -o /mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b1/taxid4577
