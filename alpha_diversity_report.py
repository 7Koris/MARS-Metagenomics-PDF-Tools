import math
import argparse
import random
from os import path
import sys

from matplotlib import colors
import pandas as pd
import read_data as rd
import numpy as np
import composition_stats as cs
import numpy as np
from scipy import stats as scistats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# import rpy2.robjects as robjects
import alphas as alphas
import taxonomy as tax

def generate_alpha_report(file_prefix, read_data: rd.ReadData):
    tax_dict = read_data.get_tax_count_dict()
    reads = [read.get_assignment() for read in read_data.reads.values()]
    
    output = PdfPages(file_prefix + ".alpha_report.pdf")
    
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
    
    
    for key in tax_table.keys():
        if key != "phylum":
            continue
        
        rank_freq_dict = tax_table[key]

        for sub_key in rank_freq_dict.keys():
            rank_freq_dict[sub_key] = rank_freq_dict[sub_key] / count_sum
            
       
        wedge, text= plots[0].pie(rank_freq_dict.values())
        labels = list(rank_freq_dict.keys())
        values = list(rank_freq_dict.values())
        for idx, label in enumerate(labels):
            label += " "
            label += str(round(values[idx] * 100, 3))
            label += "%"
            labels[idx] = label
        plots[0].legend(wedge, labels, title=key,  bbox_to_anchor=(0.9, 0, 0.5, 1))
      
    plots[1].axis('off')
    h_offset = 0.45
    v_offset = 0.25
    plots[1].text(h_offset, 1 - v_offset, "Alpha Summary: %s" %(str(path.basename(file_prefix))))
    plots[1].text(h_offset, 0.95 - v_offset, "     Richness: %s" %(round(richness, 3)))
    plots[1].text(h_offset, 0.9 - v_offset, "     Chao1: %s" %(round(chao1, 3)))
    plots[1].text(h_offset, 0.85 - v_offset, "     Shannon's Alpha: %s" %(round(shannons_alpha, 3)))
    plots[1].text(h_offset, 0.8 - v_offset, "     Berger-Parker: %s" %(round(bp, 3)))
    plots[1].text(h_offset, 0.75 - v_offset, "     Simpson's Alpha: %s" %(round(simp, 3)))
    plots[1].text(h_offset, 0.7 - v_offset, "     Inverse Simpson's Alpha: %s" %(round(in_simp, 3)))
    plots[1].text(h_offset, 0.65 - v_offset, "     Evenness: %s" %(round(evenness, 3)))
    plots[1].text(h_offset, 0.6 - v_offset, "     Fungus-Species ratio: %s" %(round(fungus_species_ratio, 3)))
    plots[1].text(h_offset, 0.55 - v_offset, "     Fungi read count: %s" %(round(fungi_count, 3)))
    plots[1].text(h_offset, 0.5 - v_offset, "     Total reads: %s" %(round(count_sum, 3)))
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    fig, plots = plt.subplots(1, 4)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    plots[0].axis('off')
    plots[2].axis('off')
    
    # phylum level heatmap
    # z_scores = scistats.zscore(list(tax_table["phylum"].values()))
    # data = []
    # for score in z_scores:
    #     data.append([score])
    # data = np.array(data)
    # im = plots[1].imshow(data, aspect='auto', cmap='hot')
    # plots[1].set_title("Phylum Z-Scores")
    # plots[1].get_xaxis().set_visible(False)
    # plots[1].set_yticks(np.arange(len(list(tax_table["phylum"].keys()))), labels=list(tax_table["phylum"].keys()))
    # plots[1].figure.colorbar(im)
    
    # # superkingdom level heatmap
    # z_scores = scistats.zscore(list(tax_table["superkingdom"].values()))
    # data = []
    # for score in z_scores:
    #     data.append([score])
    # data = np.array(data)
    # im = plots[3].imshow(data, aspect='auto', cmap='hot')
    # plots[3].set_title("Superkingdom Z-Scores")
    # plots[3].get_xaxis().set_visible(False)
    # plots[3].set_yticks(np.arange(len(list(tax_table["superkingdom"].keys()))), labels=list(tax_table["superkingdom"].keys()))
    # plots[3].figure.colorbar(im)
    
    # output.savefig()
    # plt.close()
    
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    #plots.axis('off')
    
    # genus level heatmap
    z_scores = scistats.zscore(list(tax_table["genus"].values()))
    labels = list(tax_table["genus"].keys())
    data_dict = {}
    for idx, label in enumerate(labels):
        data_dict[label] = z_scores[idx]
        
    # sort dict
    data_dict = dict(sorted(data_dict.items(), key=lambda item: item[1], reverse=True))
    data = []
    data_labels = []
    
    for score in data_dict.values():
        data.append([score])
        data_labels.append(list(data_dict.keys())[list(data_dict.values()).index(score)])

 
    top_ten_genus = data[:10]
    top_ten_genus_labels = data_labels[:10]
    top_ten_counts = []
    for label in top_ten_genus_labels:
        top_ten_counts.append(tax_table["genus"][label])

    for idx, label in enumerate(top_ten_genus_labels):
        top_ten_genus_labels[idx] += " - "
        top_ten_genus_labels[idx] += str(top_ten_counts[idx])
        top_ten_genus_labels[idx] += " reads"

    data = np.array(top_ten_genus)
    im = plots.imshow(top_ten_genus, aspect='auto', cmap='hot', norm=colors.Normalize(vmin=min(z_scores), vmax=max(z_scores)))
    plots.set_title("Genus Z-Scores top 10")
    plots.get_xaxis().set_visible(False)
    plots.set_yticks(np.arange(len(top_ten_genus_labels)), labels=top_ten_genus_labels)
    plots.figure.colorbar(im)
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    # Graph rarefaction curves
    filtered_reads = read_ids_genus.copy()
    # for idx, level, unit, read_i, identitiy_score, length in mu_data.mapping_units.itertuples():
    #     if unit in mu_data.mapping_unit_2_tax_id: # if unit was not filtered
    #         id = mu_data.mapping_unit_2_tax_id[unit]
    #         filtered_reads.append((unit, id, identitiy_score, length))
        
    reads_processed = 0
    read_per_species = {}
    read_per_otu = {}
    read_per_genus = {}
    
    num_species_rare = []
    num_genus_rare = []
    alphas_chao1_rare = []
    alphas_shannon_rare = []
    alphas_evenness_rare = []
    
    while len(filtered_reads) > 0:
        index = random.choice(range(len(filtered_reads)))
        (id) = filtered_reads.pop(index)
        reads_processed += 1
        
        # Compute new values
        if id not in read_per_species:
            read_per_species[id] = 0
            
        if id not in read_per_genus:
            read_per_genus[id] = 0
        # if unit not in read_per_otu:
        #     read_per_otu[unit] = 0
        
        read_per_species[id] += 1
        read_per_genus[id] += 1
        # read_per_otu[unit] += 1

        
        # compute current number of species
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
        
        
        # compute current number of otus
        # num_otus = 0
        # for unit in read_per_otu.keys():
        #     if read_per_otu[unit] > 0:
        #         num_otus += 1
        # num_otus_rare.append(num_otus)
        
            
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
        
    # num species over reads
    plt.figure(figsize=(12, 8))
    plt.title("Genus Count")
    plt.plot(range(1, reads_processed+1), num_species_rare)
    plt.xlabel("Reads Processed")
    plt.ylabel("Number of Species")
    output.savefig()
    plt.close()
    
    # # num otus over reads
    # plt.figure(figsize=(12, 8))
    # plt.title("Number of Units")
    # plt.plot(range(1, reads_processed+1), num_otus_rare)
    # plt.xlabel("Reads Processed")
    # plt.ylabel("Number of Units")
    # output.savefig()
    # plt.close()
    
    # num species over reads
    plt.figure(figsize=(12, 8))
    plt.title("Chao1")
    plt.plot(range(1, reads_processed+1), alphas_chao1_rare)
    plt.xlabel("Reads Processed")
    plt.ylabel("Chao1")
    output.savefig()
    plt.close()
    
    # Shannon's Alpha
    plt.figure(figsize=(12, 8))
    plt.title("Shannon's Alpha")
    plt.plot(range(1, reads_processed+1), alphas_shannon_rare)
    plt.xlabel("Reads Processed")
    plt.ylabel("Shannon's Alpha")
    output.savefig()
    plt.close()
    
    # num species over reads
    plt.figure(figsize=(12, 8))
    plt.title("Evenness")
    plt.plot(range(1, reads_processed+1), alphas_evenness_rare, label="Number of Species")
    plt.xlabel("Reads Processed")
    plt.ylabel("Evenness")
    output.savefig()
    plt.close()
    print("Alpha report written to " + output._file.fh.name)
    output.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alpha Diversity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    #parser.add_argument("-M", "--mtsv", action="store_true", help="Classifcation file is MTSV output")
    parser.add_argument("-mtsv", "--mtsv-file",  type=str, help="MTSV File")
    parser.add_argument("-mtsvl", "--mtsv-lookup-file",  type=str, help="MTSV Lookup File")
    parser.add_argument("-meta", "--meta-maps-file",  type=str, help="MetaMaps File")
    parser.add_argument("-I", "--ignore-ids", action="store_true", help="Ignore ids in the file .ignoreids")
    parser.add_argument("-S", "--skip-coverage-filter", action="store_true", help="Skip the coverage filtering step. Saves time if you already have a .ignoreids file. Will not generate an outlier pdf.")
    parser.add_argument("-metaref", "--meta-maps-reference-file", type=str, help="MetaMaps file to filter MTSV with")
    parser.add_argument("-mtsvref", "--mtsv-reference-file", type=str, help="MTSV file to filter MetaMaps reads with")


    
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    #file_prefix = config["classification_file_prefix"]
    # max_outlier_coverage = config["max_avg_outlier_coverage"]
    # proportion = config["trim_proportion"]
    ignore_ids = config["ignore_ids"]
    skip_coverage_filter = config["skip_coverage_filter"]
    meta_file = config["meta_maps_file"]
    mtsv_file = config["mtsv_file"]
    mtsv_lookup_file = config["mtsv_lookup_file"]
    meta_ref_file = config["meta_maps_reference_file"]
    mtsv_ref_file = config["mtsv_reference_file"]
    file_prefix = None
       
    id_blacklist_dict = {}
    if ignore_ids:
        print("Ignoring ids from file .ignoreids")
        for id in open(".ignoreids", "r"):
            id = id.strip()
            if id not in id_blacklist_dict:
                id_blacklist_dict[id] = True

    read_data = None
    output_name = ""
    
    if meta_file:
        read_data = rd.ReadData()
        read_data.parse_metamaps_reads_2_taxon(meta_file)
        output_name += path.basename(meta_file)
        if mtsv_ref_file:
            read_data.parse_mtsv_reads(mtsv_ref_file, mtsv_lookup_file, True)
            read_data.prune_non_incidental_reads()
            read_data.prune_by_level("genus")
        
    if mtsv_file:
        read_data = rd.ReadData()
        read_data.parse_mtsv_reads(mtsv_file, mtsv_lookup_file)
        output_name += path.basename(mtsv_file)
        if meta_ref_file:
            read_data.parse_metamaps_reads_2_taxon(meta_ref_file, True)
            read_data.prune_non_incidental_reads()
            read_data.resolve_lca()
            read_data.prune_by_level("genus")
        
    if not mtsv_file and not meta_file:
        sys.exit("Error: Must provide either a metamaps or mtsv file")
    
    if mtsv_file and meta_file:
        sys.exit("Error: Cannot provide both a metamaps and mtsv file. Choose one. Other file must be used as a reference file with -metaref or -mtsvref")
        
    
    
    # if not skip_coverage_filter:
    #     mu_data.load_coverage()
    #     mu_data.filter_sig_bin_outliers()
    generate_alpha_report(output_name, read_data)
    