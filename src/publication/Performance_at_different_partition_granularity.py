import pandas as pd
import matplotlib.pyplot as plt

#Plot AUCROC-HM and AUCROC-HM across different partition granularity


#AUCROC-HM
ymax , ymin = 0.75 , 0.64
##PCC
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/PCC/AUCROC_HM_full_scores.csv"


data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig(datapath.split(".csv")[0] +"_PCC.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.grid(axis='x', linewidth = 0.1)
plt.ylim(0.5, 0.75)
plt.savefig(datapath.split(".csv")[0] +"_sub_PCC.svg")
plt.clf()


##SCC
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/SCC/AUCROC_HM_full_scores.csv"

data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig( datapath.split(".csv")[0] +"_SCC.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.grid(axis='x', linewidth = 0.1)
plt.ylim(0.5, 0.75)
plt.savefig(datapath.split(".csv")[0] +"_sub_SCC.svg")
plt.clf()

##bicor
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/bicor/AUCROC_HM_full_scores.csv"

data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig( datapath.split(".csv")[0] +"_bicor.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.grid(axis='x', linewidth = 0.1)
plt.ylim(0.5, 0.75)
plt.savefig(datapath.split(".csv")[0] +"_sub_bicor.svg")
plt.clf()


#AUCPRC-HM
ymax , ymin = 0.74 , 0.64
##PCC
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/PCC/AUCPRC_HM_full_scores.csv"


data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig(datapath.split(".csv")[0] +"_PCC.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.ylim(0.5, 0.75)
plt.grid(axis='x', linewidth = 0.1)
plt.savefig(datapath.split(".csv")[0] +"_sub_PCC.svg")
plt.clf()


##SCC
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/SCC/AUCPRC_HM_full_scores.csv"

data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig( datapath.split(".csv")[0] +"_SCC.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.grid(axis='x', linewidth = 0.1)
plt.ylim(0.5, 0.75)
plt.savefig(datapath.split(".csv")[0] +"_sub_SCC.svg")
plt.clf()

##bicor
datapath = "/mnt/md2/ken/correlation_networks/Plant-GCN_data/ATTED_b2/taxid3702_5h/Optimize_k/bicor/AUCPRC_HM_full_scores.csv"

data = pd.read_csv(datapath, index_col= 0 )
methods = ["Max", "Avg", "RAvg", "RWA", "RRWA"]
k_list = list(data.columns)[2:]
k_int_list = [int(k) for k in k_list]

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.7)
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
    plt.grid(axis='x', linewidth = 0.1)
plt.ylim(ymin, ymax)
plt.savefig( datapath.split(".csv")[0] +"_bicor.svg")
plt.clf()

plt.figure(figsize = (8,8))
plt.rc('ytick', labelsize=12)
plt.rc('xtick', labelsize=12)

no_partition_score = data.iloc[1,1]
plt.plot([k_int_list[0], k_int_list[-1]],[no_partition_score, no_partition_score] ,  linestyle='dashed', linewidth=0.3, c="black",dashes=(10, 10))

for method in methods:
    sub_data = data[data.index == method]
    Q1 = list(sub_data[sub_data["Stat"] == "Q1"].values)[0][2:]
    Q3 = list(sub_data[sub_data["Stat"] == "Q3"].values)[0][2:]
    Med = list(sub_data[sub_data["Stat"] == "Med"].values)[0][2:]
    Lower_err = list(Med - Q1)
    Upper_err = list(Q3 - Med)
    yerr = [Lower_err, Upper_err]
    plt.xticks(ticks = [k_int_list[i] for i in range(len(k_int_list)) if i % 2 ==0][:-1] + [k_int_list[-1]])
    #plt.errorbar(k_int_list, Med, yerr=yerr, capsize=3, elinewidth=0.3, linewidth=0.7)
    plt.plot(k_int_list, Med, linewidth=0.35)
    
    plt.fill_between(k_int_list, list(Q3), list(Q1) , alpha=0.1)
plt.grid(axis='x', linewidth = 0.1)
plt.ylim(0.5, 0.75)
plt.savefig(datapath.split(".csv")[0] +"_sub_bicor.svg")
plt.clf()

