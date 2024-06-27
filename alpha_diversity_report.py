import math
import argparse
import random
from os import path
import MappingUnitData as mu
import numpy as np
import composition_stats as cs
import numpy as np
from scipy import stats as scistats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# import rpy2.robjects as robjects
import alphas as alphas
import taxonomy as tax

def generate_alpha_report(mu_data):
    output = PdfPages(file_prefix + ".rarefaction.pdf")
    abundance_estimates = mapping_units.get_abundance_estimates()
    estimate_list = list(abundance_estimates.values())  
    count_sum = sum(estimate_list)
    richness = len(estimate_list)
    chao1 = alphas.chao1(estimate_list)
    shannons_alpha = alphas.shannons_alpha(estimate_list)
    bp = alphas.berger_parkers_alpha(estimate_list)
    simp = alphas.simpsons_alpha(estimate_list)
    in_simp = alphas.inverse_simpsons_alpha(estimate_list)
    
    evenness = shannons_alpha / math.log(richness)
    
    tax_table = tax.get_lineage_counts(abundance_estimates, False)
    fungi_count = 0
    
    if "kingdom" in tax_table:
        if "Fungi" in tax_table["kingdom"]:
            fungi_count = tax_table["kingdom"]["Fungi"]
            
    fungus_species_ratio = 0
    if fungi_count != 0:
        fungus_species_ratio = (fungi_count / (richness-fungi_count))        
    
    #plt.figure(figsize=(12, 8))
   
    #plt.axis('off')
    
    # Fam level heatmap
    
    # Genus level heatmap
    
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
    z_scores = scistats.zscore(list(tax_table["phylum"].values()))
    data = []
    for score in z_scores:
        data.append([score])
    data = np.array(data)
    im = plots[1].imshow(data, aspect='auto', cmap='hot')
    plots[1].set_title("Phylum Z-Scores")
    plots[1].get_xaxis().set_visible(False)
    plots[1].set_yticks(np.arange(len(list(tax_table["phylum"].keys()))), labels=list(tax_table["phylum"].keys()))
    plots[1].figure.colorbar(im)
    
    # superkingdom level heatmap
    z_scores = scistats.zscore(list(tax_table["superkingdom"].values()))
    data = []
    for score in z_scores:
        data.append([score])
    data = np.array(data)
    im = plots[3].imshow(data, aspect='auto', cmap='hot')
    plots[3].set_title("Superkingdom Z-Scores")
    plots[3].get_xaxis().set_visible(False)
    plots[3].set_yticks(np.arange(len(list(tax_table["superkingdom"].keys()))), labels=list(tax_table["superkingdom"].keys()))
    plots[3].figure.colorbar(im)
    
    output.savefig()
    plt.close()
    
    # Graph rarefaction curves
    filtered_reads = []
    for idx, level, unit, read_i, identitiy_score, length in mu_data.mapping_units.itertuples():
        if unit in mu_data.mapping_unit_2_tax_id: # if unit was not filtered
            id = mu_data.mapping_unit_2_tax_id[unit]
            filtered_reads.append((unit, id, identitiy_score, length))
        
    reads_processed = 0
    read_per_species = {}
    read_per_otu = {}
    
    num_species_rare = []
    num_otus_rare = []
    alphas_chao1_rare = []
    alphas_shannon_rare = []
    alphas_evenness_rare = []
    
    while len(filtered_reads) > 0:
        index = random.choice(range(len(filtered_reads)))
        (unit, id, identitiy_score, length) = filtered_reads.pop(index)
        reads_processed += 1
        
        # Compute new values
        if id not in read_per_species:
            read_per_species[id] = 0
        if unit not in read_per_otu:
            read_per_otu[unit] = 0
        
        read_per_species[id] += 1
        read_per_otu[unit] += 1

        
        # compute current number of species
        num_species = 0
        for id in read_per_species.keys():
            if read_per_species[id] > 0:
                num_species += 1
        num_species_rare.append(num_species)
        
        # compute current number of otus
        num_otus = 0
        for unit in read_per_otu.keys():
            if read_per_otu[unit] > 0:
                num_otus += 1
        num_otus_rare.append(num_otus)
        
            
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
    plt.title("Number of Species")
    plt.plot(range(1, reads_processed+1), num_species_rare)
    plt.xlabel("Reads Processed")
    plt.ylabel("Number of Species")
    output.savefig()
    plt.close()
    
    # num otus over reads
    plt.figure(figsize=(12, 8))
    plt.title("Number of OTUs")
    plt.plot(range(1, reads_processed+1), num_otus_rare)
    plt.xlabel("Reads Processed")
    plt.ylabel("Number of OTUs")
    output.savefig()
    plt.close()
    
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
    print("Rarefaction curves written to " + file_prefix + ".rarefaction.pdf")
    output.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alpha Diversity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-m", "--max-avg-outlier-coverage", type=float, help="Maximum average outlier coverage", default=0.0)
    parser.add_argument("-t", "--trim-proportion", type=float, help="Proportion of coverages to trim out", default=0.003)
    parser.add_argument("-I", "--ignore-ids", action="store_true", help="Ignore ids in the file .ignoreids")
    parser.add_argument("-S", "--skip-coverage-filter", action="store_true", help="Skip the coverage filtering step. Saves time if you already have a .ignoreids file. Will not generate an outlier pdf.")
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    file_prefix = config["classification_file_prefix"]
    max_outlier_coverage = config["max_avg_outlier_coverage"]
    proportion = config["trim_proportion"]
    ignore_ids = config["ignore_ids"]
    skip_coverage_filter = config["skip_coverage_filter"]
    
    excluded_tax_ids = []
    if ignore_ids:
        for id in open(".ignoreids", "r"):
            id = id.strip()
            if id not in excluded_tax_ids:
                excluded_tax_ids.append(id)
    
    mapping_units = mu.MappingUnitData(file_prefix, min_freq, excluded_tax_ids)
    
    if not skip_coverage_filter:
        mapping_units.load_coverage()
        mapping_units.filter_coverage_tm_outliers(max_outlier_coverage, proportion)

        
    
    generate_alpha_report(mapping_units)
    