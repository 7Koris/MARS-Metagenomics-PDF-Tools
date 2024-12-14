import math
import argparse
import os
import random
import sys
import pandas as pd
import read_data as rd
import matplotlib.pyplot as plt
import alphas as alphas
import pickle
from matplotlib.backends.backend_pdf import PdfPages
from os import path, listdir
from scipy import stats as scistats

seed = 0

def generate_alpha_report(file_prefix, read_data: rd.ReadData):
    pd.options.mode.chained_assignment = None 
    tax_dict = read_data.get_tax_count_dict()
    reads = [read.get_assignment() for read in read_data.reads.values()]
    
    if not os.path.exists("../reports"):
        os.makedirs("../reports")
            
    output = PdfPages("../reports/" + file_prefix + ".alpha_report.pdf") 
    abundance_estimates = tax_dict["genus"]
    
    read_ids = [read.get_assignment() for read in read_data.reads.values()]
    read_ids_genus = []
    for read_id in read_ids:
        new_id = tax_dict.id_2_level_id(read_id, "genus")
        read_ids_genus.append(new_id)
    
    estimate_list = list(abundance_estimates.values())  
    count_sum = sum(estimate_list)
    richness = len(estimate_list)
    chao1 = alphas.chao1(estimate_list)
    shannons_alpha = alphas.shannons_alpha(estimate_list)
    bp = alphas.berger_parkers_alpha(estimate_list)
    simp = alphas.simpsons_alpha(estimate_list)
    in_simp = alphas.inverse_simpsons_alpha(estimate_list)
    
    evenness = shannons_alpha / math.log(richness)
    
    tax_table = tax_dict
    fungi_count = 0
    
    if "kingdom" in tax_table:
        if "Fungi" in tax_table["kingdom"]:
            fungi_count = tax_table["kingdom"]["Fungi"]
            
    fungus_species_ratio = 0
    if fungi_count != 0:
        fungus_species_ratio = (fungi_count / (count_sum - fungi_count))     
        
    fig, plots = plt.subplots(1, 2)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    species_dist = list(tax_table["species"].values())
    species_dist = sorted(species_dist, reverse=True)
    read_sum = sum(species_dist)

    plots[0].axis('off')
    plots[1].bar(range(len(species_dist)), species_dist, color='blue')
    plots[1].set_title("Species-Level Read Distribution of %i reads" %(round(read_sum, 0)))
    plots[1].set_xlabel("Species")
    plots[1].set_ylabel("Reads Assigned")
    plots[1].set_yscale('log')
    plt.axes().axis('off')
    h_offset = 0
    v_offset = 0.05
    plt.text(h_offset, 1 - v_offset, "Alpha Summary: %s" %(str(path.basename(file_prefix))))
    plt.text(h_offset, 0.95 - v_offset, "     Richness: %s" %(round(richness, 3)))
    plt.text(h_offset, 0.9 - v_offset, "     Chao1: %s" %(round(chao1, 3)))
    plt.text(h_offset, 0.85 - v_offset, "     Shannon's Alpha: %s" %(round(shannons_alpha, 3)))
    plt.text(h_offset, 0.8 - v_offset, "     Berger-Parker: %s" %(round(bp, 3)))
    plt.text(h_offset, 0.75 - v_offset, "     Simpson's Alpha: %s" %(round(simp, 3)))
    plt.text(h_offset, 0.7 - v_offset, "     Inverse Simpson's Alpha: %s" %(round(in_simp, 3)))
    plt.text(h_offset, 0.65 - v_offset, "     Evenness: %s" %(round(evenness, 3)))
    plt.text(h_offset, 0.6 - v_offset, "     Fungus-Species ratio: %s" %(round(fungus_species_ratio, 3)))
    plt.text(h_offset, 0.55 - v_offset, "     Fungi read count: %s" %(round(fungi_count, 3)))
    plt.text(h_offset, 0.5 - v_offset, "     Total filtered reads: %s" %(round(count_sum, 3)))
    plt.text(h_offset, 0.45 - v_offset, "     Total unfiltered reads: %s" %(round(read_data.processed_read_count, 3)))
    plt.tight_layout()
    output.savefig()
    plt.close()
       
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    
    # genus level heatmap of z-scores
    z_scores = scistats.zscore(list(tax_table["genus"].values()))
    labels = list(tax_table["genus"].keys())
    data_dict = {}
    for idx, label in enumerate(labels):
        data_dict[label] = z_scores[idx]
        
    data_dict = dict(sorted(data_dict.items(), key=lambda item: item[1], reverse=True))
    data = []
    data_labels = []
    
    for score in data_dict.values():
        data.append([score])
        data_labels.append(list(data_dict.keys())[list(data_dict.values()).index(score)])

    top_ten_genus = data[:10]
    top_ten_genus_labels = data_labels[:10]
    old_top_ten_labels = top_ten_genus_labels.copy()
    top_ten_counts = []
    for label in top_ten_genus_labels:
        top_ten_counts.append(tax_table["genus"][label])

    for idx, label in enumerate(top_ten_genus_labels):
        top_ten_genus_labels[idx] += " - "
        top_ten_genus_labels[idx] += str(top_ten_counts[idx])
        top_ten_genus_labels[idx] += " reads"
    
    # Bar plot 6 top ten genus bar plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    
    top_ten_counts = top_ten_counts[::-1]
    old_top_ten_labels = old_top_ten_labels[::-1]
    top_ten_genus_df = pd.DataFrame(top_ten_counts, columns=["Count"], index=old_top_ten_labels)
    sum_top = top_ten_genus_df["Count"].sum()
    other_count = len(reads) - sum_top
    top_ten_genus_df.loc["Other"] = other_count
    top_ten_genus_df.sort_values(by="Count", ascending=False, inplace=True)
    
    plot = top_ten_genus_df.T.plot.barh(stacked=True, legend=False, figsize=(12, 8), ax=ax1)
    plot.get_xaxis().set_ticks([])
    plot.get_yaxis().set_ticks([])
    labels = list(top_ten_genus_df.index)
    values = list(top_ten_genus_df["Count"].values)
    plt.legend(labels)
    handles, labels = plot.get_legend_handles_labels()
    
    for idx, label in enumerate(labels):
        labels[idx] += " - "
        labels[idx] += str(round(values[idx] / len(reads) * 100, 3)) + "%"
        
    plot.legend(handles[::-1], labels[::-1], title='Genus', loc='upper left', prop={'size': 6})
    plot.set_title("Top 10 Genus")    
    
    # Bar plot #1:  In relation to total reads: percentage of reads classified; percentage of reads unclassified (which would include reads filtered-out)
    classified_read_count = len(read_data.reads)
    unclassified_read_count = read_data.processed_read_count - classified_read_count
    classified_read_percentage = classified_read_count / read_data.processed_read_count
    unclassified_read_percentage = unclassified_read_count / read_data.processed_read_count
    read_count_data = [unclassified_read_count, classified_read_count]
    read_count_labels = ["Unclassified", "Classified"]
    read_percent_data = [unclassified_read_percentage, classified_read_percentage]
    read_count_df = pd.DataFrame(read_count_data, columns=["Count"], index=read_count_labels)
    plot = read_count_df.T.plot.barh(stacked=True, legend=False, figsize=(12, 8), ax=ax2)
    plot.get_xaxis().set_ticks([])
    plot.get_yaxis().set_ticks([])
    plt.legend(old_top_ten_labels)
    handles, labels = plot.get_legend_handles_labels()
    for idx, label in enumerate(labels):
        labels[idx] += " - "
        labels[idx] += str(round(read_percent_data[idx] * 100, 3)) + "%"
    plot.legend(handles[::-1], labels[::-1], title='Type', loc='upper left', prop={'size': 6})
    plot.set_title("Read Classification")
    
    # Bar plot #2:  In relation to reads classified:  percentage of reads bacterial; percentage of reads fungal; percentage of reads archaeal
    fungi_count = 0
    bacteria_count = 0
    archaea_count = 0
    other_eukaryota_count = 0
    if "kingdom" in tax_table:
        if "Fungi" in tax_table["kingdom"]:
            fungi_count = tax_table["kingdom"]["Fungi"]

    if "superkingdom" in tax_table:
        if "Bacteria" in tax_table["superkingdom"]:
            bacteria_count = tax_table["superkingdom"]["Bacteria"]
        if "Archaea" in tax_table["superkingdom"]:
            archaea_count = tax_table["superkingdom"]["Archaea"]
        if "Eukaryota" in tax_table["superkingdom"]:
            other_eukaryota_count = tax_table["superkingdom"]["Eukaryota"] - fungi_count
    k_data = [bacteria_count, fungi_count, archaea_count, other_eukaryota_count]
    k_labels = ["Bacteria", "Fungi", "Archaea", "Other Eukaryota"]
    k_df = pd.DataFrame(k_data, columns=["Count"], index=k_labels)
    k_plot = k_df.T.plot.barh(stacked=True, legend=False, figsize=(12, 8), ax=ax3)
    k_plot.get_xaxis().set_ticks([])
    k_plot.get_yaxis().set_ticks([])
    plt.legend(old_top_ten_labels)
    handles, labels = k_plot.get_legend_handles_labels()
    for idx, label in enumerate(labels):
        labels[idx] += " - "
        labels[idx] += str(round(k_data[idx] / len(reads) * 100, 3)) + "%"
    k_plot.legend(handles[::-1], labels[::-1], title='Kingdom', loc='upper left', prop={'size': 6})
    k_plot.set_title("Kingdom Classification")
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    
    # Bar plot #3:  In relation to reads classified as bacterial:  percentage breakdown of top-10 genera, with another added “other” category to bring-up total to 100%
    bacteria_read_data = rd.ReadData()
    for read in read_data.reads.values():
        bacteria_read_data.insert_read(read)
    bacteria_read_data.prune_reads_not_in_rank("superkingdom", "Bacteria")
    bacteria_tax_dict = bacteria_read_data.get_tax_count_dict()
    bacteria_genus_data = list(bacteria_tax_dict["genus"].values())
    bacteria_genus_labels = list(bacteria_tax_dict["genus"].keys())
    bacteria_genus_df = pd.DataFrame(bacteria_genus_data, columns=["Count"], index=bacteria_genus_labels)
    bacteria_genus_df.sort_values(by="Count", ascending=False, inplace=True)
    top_ten_genus = bacteria_genus_df[:10]
    remaining_genus = bacteria_genus_df[10:]
    remaining_genus_sum = remaining_genus.sum()
    other_count = remaining_genus_sum
    top_ten_genus.loc["Other"] = other_count
    total_sum = top_ten_genus["Count"].sum()
    top_ten_genus.sort_values(by="Count", ascending=False, inplace=True)
    df_values = top_ten_genus["Count"].values
    plot = top_ten_genus.T.plot.barh(stacked=True, legend=False,  figsize=(12, 8), ax=ax1)
    plot.get_xaxis().set_ticks([])
    plot.get_yaxis().set_ticks([])
    plt.legend(labels)
    handles, labels = plot.get_legend_handles_labels()
    
    for idx, label in enumerate(labels):
        labels[idx] += " - "
        labels[idx] += str(round(float(df_values[idx] / total_sum * 100), 3)) + "%"
        
    plot.legend(handles[::-1], labels[::-1], title='Genus', loc='upper left', prop={'size': 6})
    plot.set_title("Top 10 Genus of Kingdom Bacteria")    
    
    # Bar plot #4:  In relation to reads classified as fungal:  percentage breakdown of top-10 genera, with another added “other” category to bring-up total to 100%
    fungi_read_data = rd.ReadData()
    for read in read_data.reads.values():
        fungi_read_data.insert_read(read)
    fungi_read_data.prune_reads_not_in_rank("kingdom", "Fungi")
    if len(fungi_read_data.reads) != 0:
        fungi_tax_dict = fungi_read_data.get_tax_count_dict()
        fungi_genus_data = list(fungi_tax_dict["genus"].values())
        fungi_genus_labels = list(fungi_tax_dict["genus"].keys())
        fungi_genus_df = pd.DataFrame(fungi_genus_data, columns=["Count"], index=fungi_genus_labels)
        fungi_genus_df.sort_values(by="Count", ascending=False, inplace=True)
        top_ten_genus = fungi_genus_df[:10]
        remaining_genus = fungi_genus_df[10:]
        remaining_genus_sum = remaining_genus.sum()
        other_count = remaining_genus_sum
        total_sum = remaining_genus_sum + top_ten_genus["Count"].sum()
        top_ten_genus.loc["Other"] = other_count
        top_ten_genus.sort_values(by="Count", ascending=False, inplace=True)
        df_values = top_ten_genus["Count"].values
        plot = top_ten_genus.T.plot.barh(stacked=True, legend=False,  figsize=(12, 8), ax=ax2)
        plot.get_xaxis().set_ticks([])
        plot.get_yaxis().set_ticks([])
        plt.legend(labels)
        handles, labels = plot.get_legend_handles_labels()
        su_a = 0
        for idx, label in enumerate(labels):
            labels[idx] += " - "
            labels[idx] += str(round(float(df_values[idx] / total_sum * 100), 3)) + "%"
            su_a += round(df_values[idx] / total_sum * 100, 3)
        print(su_a)
        plot.legend(handles[::-1], labels[::-1], title='Genus', loc='upper left', prop={'size': 6})
        plot.set_title("Top 10 Genus of Kingdom Fungi")
    else:
        ax2.axis('off')
    
    # Bar plot #5:  In relation to reads classified as archaeal:  percentage breakdown of top-10 genera, with another added “other” category to bring-up total to 100%
    archael_read_data = rd.ReadData()
    for read in read_data.reads.values():
        archael_read_data.insert_read(read)
    archael_read_data.prune_reads_not_in_rank("superkingdom", "Archaea")
    
    if len(archael_read_data.reads) != 0:
        archael_tax_dict = archael_read_data.get_tax_count_dict()
        archael_genus_data = list(archael_tax_dict["genus"].values())
        archael_genus_labels = list(archael_tax_dict["genus"].keys())
        archael_genus_df = pd.DataFrame(archael_genus_data, columns=["Count"], index=archael_genus_labels)
        archael_genus_df.sort_values(by="Count", ascending=False, inplace=True)
        top_ten_genus = archael_genus_df[:10]
        remaining_genus = archael_genus_df[10:]
        remaining_genus_sum = remaining_genus.sum()
        other_count = remaining_genus_sum
        total_sum = remaining_genus_sum + top_ten_genus["Count"].sum()
        top_ten_genus.loc["Other"] = other_count
        top_ten_genus.sort_values(by="Count", ascending=False, inplace=True)
        print(top_ten_genus)
        df_values = top_ten_genus["Count"].values
        print(df_values)
        plot = top_ten_genus.T.plot.barh(stacked=True, legend=False,  figsize=(12, 8), ax=ax3)
        plot.get_xaxis().set_ticks([])
        plot.get_yaxis().set_ticks([])
        plt.legend(labels)
        handles, labels = plot.get_legend_handles_labels()
        su_a = 0
        for idx, label in enumerate(labels):
            labels[idx] += " - "
            labels[idx] += str(float(round(df_values[idx] / total_sum * 100, 3))) + "%"
            su_a += round(df_values[idx] / total_sum * 100, 3)
        print(su_a)
        plot.legend(handles[::-1], labels[::-1], title='Genus', loc='upper left', prop={'size': 6})
        plot.set_title("Top 10 Genus of Kingdom Archaea")
    else:
        ax3.axis('off')
    plt.tight_layout()
    output.savefig()
    plt.close()  
    
    filtered_reads = read_ids_genus.copy()    
    reads_processed = 0
    read_per_species = {}
    read_per_genus = {}
    num_species_rare = []
    num_genus_rare = []
    alphas_chao1_rare = []
    alphas_shannon_rare = []
    alphas_evenness_rare = []
   
    random.seed(seed)
    
    num_reads = 0
    while len(filtered_reads) > 0:
        num_reads += 1
        index = random.choice(range(len(filtered_reads)))
        (id) = filtered_reads.pop(index)
        reads_processed += 1
        

        if id not in read_per_species:
            read_per_species[id] = 0
            
        if id not in read_per_genus:
            read_per_genus[id] = 0
        
        read_per_species[id] += 1
        read_per_genus[id] += 1
        
        num_species = 0
        for id in read_per_species.keys():
            if read_per_species[id] > 0:
                num_species += 1
        num_species_rare.append(num_species)
        
        num_genus = 0
        for id in read_per_genus.keys():
            if read_per_genus[id] > 0:
                num_genus += 1
        num_genus_rare.append(num_genus)        
            
        estimate_list = list(read_per_species.values())  
        alphas_chao1_rare.append(alphas.chao1(estimate_list))
        if len(estimate_list) == 1:
            sa = 0
        else:
            sa = alphas.shannons_alpha(estimate_list)
        alphas_shannon_rare.append(sa)
        if math.log(num_species) == 0:
            alphas_evenness_rare.append(0)
        else:
            alphas_evenness_rare.append( sa / math.log(num_species))
        
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    
    # num taxa over reads
    ax1.set_title("Taxa Count")
    ax1.plot(range(1, reads_processed+1), num_species_rare)
    ax1.set_xlabel("Reads Processed")
    ax1.set_ylabel("Taxa Count")
    
    # chao1 over reads
    ax2.set_title("Chao1")
    ax2.plot(range(1, reads_processed+1), alphas_chao1_rare)
    ax2.set_xlabel("Reads Processed")
    ax2.set_ylabel("Chao1")
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    
    # Shannon's Alpha
    ax1.set_title("Shannon's Alpha over Reads Processed")
    ax1.plot(range(1, reads_processed+1), alphas_shannon_rare)
    ax1.set_xlabel("Reads Processed")
    ax1.set_ylabel("Shannon's Alpha")
    
    # num species over reads
    ax2.set_title("Evenness over Reads Processed")
    ax2.plot(range(1, reads_processed+1), alphas_evenness_rare, label="Number of Species")
    ax2.set_xlabel("Reads Processed")
    ax2.set_ylabel("Evenness")
    plt.tight_layout()
    output.savefig()
    plt.close()
    print("Alpha report written to " + output._file.fh.name)
    output.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alpha Diversity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file", type=str, help="Pickle file name", default=None)
    parser.add_argument("-s", "--seed", type=str, help="Seed for random number generator", default=0)

    
    args = parser.parse_args()
    config = vars(args)
    file = config["file"]
    seed = config["seed"]
    
    if file is None:
        files = listdir("../reads")
        for file in files:
            if path.splitext(file)[1] != ".p":
                continue
            file = path.splitext(file)[0] + ".p"
            break
        if file is None:
            sys.exit("No files found in reads ~/reads directory.") 
            
    data = pickle.load(open("../reads/" + file, "rb"))
    file = path.splitext(file)[0]
    read_data = rd.ReadData()
    read_data.load_data(data)
    
    generate_alpha_report(file, read_data)
    