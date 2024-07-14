import math
import argparse
import random
from os import path

from matplotlib import colors
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
import betas as betas

def generate_beta_report(output_name, all_read_data_names, all_read_data):
    output = PdfPages(output_name + ".beta_report.pdf")
    
    all_fungi_ratios = []
    all_fungi_counts = []
    all_phylum_scores = []
    all_abundance_estimates = {}
    for idx, read_data in enumerate(all_read_data):
        current_name = all_read_data_names[idx]
        tax_dict = read_data.get_tax_dict()
        
        abundance_estimates = tax_dict["species"]
        reads = [read.get_assignment() for read in read_data.reads.values()]
    

        all_abundance_estimates[current_name] = abundance_estimates
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
            
        all_fungi_ratios.append(fungus_species_ratio)
        all_fungi_counts.append(fungi_count)  
            
        # phylum level z-scores
        z_scores = scistats.zscore(list(tax_table["phylum"].values()))
        labels = list(tax_table["phylum"].keys())
        z_dict = {}
        for i in range(len(labels)):
            z_dict[labels[i]] = z_scores[i]
        all_phylum_scores.append(z_dict)
   
    # Collect all phylums 
    all_phylums = {}
    for phylum_scores in all_phylum_scores:
        for phylum in phylum_scores.keys():
            if phylum not in all_phylums:
                all_phylums[phylum] = 1

    # Insert 0 values for missing entries in the z-scores
    for phylum_scores in all_phylum_scores:
        for phylum in all_phylums.keys():
            if phylum not in phylum_scores:
                phylum_scores[phylum] = 0
    
    # Collect all ids in the abundance estimates
    all_ids = {}
    for prefix in all_abundance_estimates:
        abundance_estimates = all_abundance_estimates[prefix]
        for id in abundance_estimates.keys():
            if id not in all_ids:
                all_ids[id] = 1

    # Insert 0 values for missing entries in the abundance estimates
    for prefix in all_abundance_estimates:
        abundance_estimates = all_abundance_estimates[prefix]
        for id in all_ids.keys():
            if id not in abundance_estimates:
                abundance_estimates[id] = 0
                
    all_unit_sets = all_read_data
    all_unit_names = all_read_data_names
    # Beta matrices
    bc_matrix = np.zeros((len(all_unit_sets), len(all_unit_sets)))
    jaccard_matrix = np.zeros((len(all_unit_sets), len(all_unit_sets)))
    city_block_matrix = np.zeros((len(all_unit_sets), len(all_unit_sets)))
    aitchison_matrix = np.zeros((len(all_unit_sets), len(all_unit_sets)))
    for current_name in all_unit_names:
        for other_name in all_unit_names:
            if current_name == other_name:
                bc_matrix[all_unit_names.index(current_name)][all_unit_names.index(other_name)] = 0 # same sample
                continue
            current_abundance_estimates = all_abundance_estimates[current_name]
            other_abundance_estimates = all_abundance_estimates[other_name]
            
            bc_matrix[all_unit_names.index(current_name)][all_unit_names.index(other_name)] = betas.bray_curtis_distance(current_abundance_estimates, other_abundance_estimates, all_ids)
            jaccard_matrix[all_unit_names.index(current_name)][all_unit_names.index(other_name)] = betas.jaccard_distance(current_abundance_estimates, other_abundance_estimates, all_ids)
            city_block_matrix[all_unit_names.index(current_name)][all_unit_names.index(other_name)] = betas.city_block_distance(current_abundance_estimates, other_abundance_estimates, all_ids)
            #aitchison_matrix[all_unit_names.index(current_name)][all_unit_names.index(other_name)] = betas.aitchison_distance(current_abundance_estimates, other_abundance_estimates, all_ids)
       

    # Phylum level heatmap comparison between samples
    data = []
    y_labels = []
    for phylum in all_phylums.keys():
        y_labels.append(phylum)
        row = []
        for prefix in all_unit_names:
            row.append(all_phylum_scores[all_unit_names.index(prefix)][phylum])
        data.append(row)
    x_labels = all_unit_names
    
    # Rarefaction curves
    curves_dict = {}
    for idx, current_read_data in enumerate(all_read_data):
        current_name = all_read_data_names[idx]
        
        filtered_reads = [read.get_assignment() for read in current_read_data.reads.values()]
            
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
            (id) = filtered_reads.pop(index)
            reads_processed += 1
            
            # Compute new values
            if id not in read_per_species:
                read_per_species[id] = 0
            # if unit not in read_per_otu:
            #     read_per_otu[unit] = 0
            
            read_per_species[id] += 1
            #read_per_otu[unit] += 1

            
            # compute current number of species
            num_species = 0
            for id in read_per_species.keys():
                if read_per_species[id] > 0:
                    num_species += 1
            num_species_rare.append(num_species)
            
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
        curves_dict[current_name] = {}
        curves_dict[current_name]["num_species"] = num_species_rare
        curves_dict[current_name]["num_otus"] = num_otus_rare
        curves_dict[current_name]["chao1"] = alphas_chao1_rare
        curves_dict[current_name]["shannon"] = alphas_shannon_rare
        curves_dict[current_name]["evenness"] = alphas_evenness_rare
    
    # Plot phylum level heatmap
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    im = plots.imshow(data, aspect='auto', cmap='twilight_shifted')
    plots.set_title("Phylum Z-Scores")
    plots.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45, size='small')
    plots.set_yticks(np.arange(len(y_labels)), labels=y_labels) 
    plots.figure.colorbar(im)
    #plt.show()
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    # Plot Bray-Curtis matrix heatmap
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    im = plots.imshow(bc_matrix, aspect='auto')
    plots.set_title("Bray-Curtis Distance")
    plots.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45, size='small')
    plots.set_yticks(np.arange(len(x_labels)), labels=x_labels, size='small')
    plots.figure.colorbar(im)
    mean = np.mean(bc_matrix)
    
    for i in range(len(bc_matrix)):
        for j in range(len(bc_matrix)):
            text = plots.text(j, i, bc_matrix[i, j], ha="center", va="center", color='white' if bc_matrix[i, j] < mean else 'black')
    
    plt.tight_layout()
    output.savefig()
    plt.close()
    
    # Plot Jaccard matrix heatmap
    # fig, plots = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    # im = plots.imshow(jaccard_matrix, aspect='auto')
    # plots.set_title("Jaccard Distance")
    # plots.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45, size='small')
    # plots.set_yticks(np.arange(len(x_labels)), labels=x_labels, size='small')
    # plots.figure.colorbar(im)
    # mean = np.mean(jaccard_matrix)
    
    # for i in range(len(jaccard_matrix)):
    #     for j in range(len(jaccard_matrix)):
    #         text = plots.text(j, i, jaccard_matrix[i, j], ha="center", va="center", color='white' if jaccard_matrix[i, j] < mean else 'black')
    
    # plt.tight_layout()
    # output.savefig()
    
    # Plot City Block matrix heatmap
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    im = plots.imshow(city_block_matrix, aspect='auto')
    plots.set_title("City Block Distance")
    plots.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45, size='small')
    plots.set_yticks(np.arange(len(x_labels)), labels=x_labels, size='small')
    plots.figure.colorbar(im)
    mean = np.mean(city_block_matrix)
    
    for i in range(len(city_block_matrix)):
        for j in range(len(city_block_matrix)):
            text = plots.text(j, i, city_block_matrix[i, j], ha="center", va="center", color='white' if city_block_matrix[i, j] < mean else 'black')
    
    plt.tight_layout()
    output.savefig()
    
    # # Plot Aitchison matrix heatmap
    # fig, plots = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    # im = plots.imshow(jaccard_matrix, aspect='auto')
    # plots.set_title("Aitchison Distance")
    # plots.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45, size='small')
    # plots.set_yticks(np.arange(len(x_labels)), labels=x_labels, size='small')
    # plots.figure.colorbar(im)
    # mean = np.mean(aitchison_matrix)
    
    # for i in range(len(aitchison_matrix)):
    #     for j in range(len(aitchison_matrix)):
    #         text = plots.text(j, i, aitchison_matrix[i, j], ha="center", va="center", color='white' if aitchison_matrix[i, j] < mean else 'black')
    
    # plt.tight_layout()
    # output.savefig()
    plt.close()
    
    # PCoA
    # https://medium.com/@conniezhou678/applied-machine-learning-part-12-principal-coordinate-analysis-pcoa-in-python-5acc2a3afe2d
    # fig, plots = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    # mds_matrix = MDS(n_components=2, dissimilarity='precomputed').fit_transform(bc_matrix)
    
    
    
    # #plots.scatter(mds_matrix[:, 0], mds_matrix[:, 1])
    # for i, txt in enumerate(all_unit_names):
    #     plots.annotate(txt, (mds_matrix[i, 0], mds_matrix[i, 1]))
    # plots.set_xlabel('PC1')
    # plots.set_ylabel('PC2')
    # plots.set_title('PCoA')
    # plt.tight_layout()
    
    # output.savefig()

    
    # PCA
    # fig, plots = plt.subplots(1, 1)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    # mds_matrix = MDS(n_components=2, dissimilarity='precomputed').fit_transform(bc_matrix)
    # plots.scatter(mds_matrix[:, 0], mds_matrix[:, 1])
    # for i, txt in enumerate(all_unit_names):
    #     plots.annotate(txt, (mds_matrix[i, 0], mds_matrix[i, 1]))
    # plots.set_xlabel('PC1')
    # plots.set_ylabel('PC2')
    # plots.set_title('PCA')
    # plt.tight_layout()
    
    # output.savefig()
    # output.close()
    
    # return 
    
    # Plot rarefaction curves    
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    for name in curves_dict.keys():
        plots.plot(curves_dict[name]["num_species"], label=name)
    plots.set_title("Number of Species")
    #plots.xlabel("Reads Processed")
    plt.legend(loc='upper left')
    output.savefig()
    plt.close()
    
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    for name in curves_dict.keys():
        plots.plot(curves_dict[name]["chao1"], label=name)
    plots.set_title("Chao1")
    #plots.xlabel("Reads Processed")
    plt.legend(loc='upper left')
    output.savefig()
    plt.close()
    
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    for name in curves_dict.keys():
        plots.plot(curves_dict[name]["shannon"], label=name)
    plots.set_title("Shannon's Alpha")
    #plots.xlabel("Reads Processed")
    plt.legend(loc='upper left')
    output.savefig()
    plt.close()
    
    fig, plots = plt.subplots(1, 1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    for name in curves_dict.keys():
        plots.plot(curves_dict[name]["evenness"], label=name)
    plots.set_title("Evenness")
    #.xlabel("Reads Processed")
    plt.legend(loc='upper left')
    output.savefig()
    plt.close()
    
    output.close()
        
  
        
    # fig, plots = plt.subplots(1, 4)
    # fig.set_figwidth(12)
    # fig.set_figheight(8)
    # plots[0].axis('off')
    # plots[2].axis('off')
    
    # # phylum level heatmap
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
        
       
    
    # # num species over reads
    # plt.figure(figsize=(12, 8))
    # plt.title("Evenness")
    # #plt.plot(range(1, reads_processed+1), alphas_evenness_rare, label="Number of Species")
    # plt.xlabel("Reads Processed")
    # plt.ylabel("Evenness")
    # output.savefig()
    # plt.close()
    # print("Rarefaction curves written to " + file_prefix + ".rarefaction.pdf")
    # output.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Beta Diversity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-I", "--ignore-ids", action="store_true", help="Ignore ids in the file .ignoreids")
    #parser.add_argument("-S", "--skip-coverage-filter", action="store_true", help="Skip the coverage filtering step. Saves time if you already have a .ignoreids file. Will not generate an outlier pdf.")
    parser.add_argument("-mtsv", "--mtsv-files", nargs='+', help="MTSV files")
    parser.add_argument("-mtsvl", "--mtsv-lookup-files", nargs='+', help="MTSV lookup files (IN ORDER OF MTSV FILES)")
    parser.add_argument("-meta", "--meta-maps-files", nargs='+', help="MetaMaps files")
    
    parser.add_argument("-metaref", "--meta-maps-reference-files", nargs='+', help="MetaMaps files to filter MTSV with (IN ORDER OF MTSV FILES)")
    parser.add_argument("-mtsvref", "--mtsv-reference-files", nargs='+', help="MTSV files to filter MetaMaps with (IN ORDER OF MTSV FILES)")
    parser.add_argument("-mtsvrefl", "--mtsvref-lookup-files", nargs='+', help="Lookup files for MTSV refs (IN ORDER OF MTSV FILES)")

    
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    mtsv_files = config["mtsv_files"]
    meta_files = config["meta_maps_files"]
    meta_ref_files = config["meta_maps_reference_files"]
    mtsv_ref_files = config["mtsv_reference_files"]
    mtsv_lookup_files = config["mtsv_lookup_files"]
    mtsv_ref_lookup_files = config["mtsvref_lookup_files"]
    ignore_ids = config["ignore_ids"]
    #skip_coverage_filter = config["skip_coverage_filter"]
    
    output_name = ""
    
    all_read_data = []
    
    all_meta_refs = []
    all_mtsv_refs = []
    all_mtsv_ref_lookup_files = []
    all_lookup_files = []
    all_read_data_names = []
    
    if not meta_ref_files is None:
        for file in meta_ref_files:
            all_meta_refs.append(file)
    
    if not mtsv_ref_files is None:
        for file in mtsv_ref_files:
            all_mtsv_refs.append(file)

    if not mtsv_lookup_files is None:
        for file in mtsv_lookup_files:
            all_lookup_files.append(file)

    if not mtsv_ref_lookup_files is None:
        for file in mtsv_ref_lookup_files:
            all_mtsv_ref_lookup_files.append(file)
    
    if not meta_files is None:
        for idx, file in enumerate(meta_files):
            if len(output_name) > 0:
                output_name += "_"
            output_name += path.basename(file)
            read_data = rd.ReadData()
            read_data.parse_metamaps_reads_2_taxon(file)
            
            if len(all_mtsv_refs) > 0:
                read_data.parse_mtsv_reads(all_mtsv_refs[idx], all_mtsv_ref_lookup_files[idx], True)
                read_data.prune_non_incidental_reads()
            
            all_read_data.append(read_data)
            all_read_data_names.append(path.basename(file))
    
    if not mtsv_files is None:
        for idx, file in enumerate(mtsv_files):
            if len(output_name) > 0:
                output_name += "_"
            output_name += path.basename(file)
            read_data = rd.ReadData()
            read_data.parse_mtsv_reads(file, all_lookup_files[idx])
            
                        
            if len(all_meta_refs) > 0:
                read_data.parse_metamaps_reads_2_taxon(all_meta_refs[0], True)
                print("non incidental reads", len(read_data._incidence_dict))
                read_data.prune_non_incidental_reads()
            read_data.resolve_lca()
            
            all_read_data.append(read_data)
            all_read_data_names.append(path.basename(file))
        
                
  

    generate_beta_report(output_name, all_read_data_names, all_read_data)
    