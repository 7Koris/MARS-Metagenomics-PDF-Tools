import sys
import re
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: Main.py <file prefix> <min plot frequency>")
        sys.exit(1)
        
    file_prefix = sys.argv[1]
    min_plot_freq = float(sys.argv[2])
    
    lengths_file = file_prefix + ".EM.lengthAndIdentitiesPerMappingUnit"
    coverage_file = file_prefix + ".EM.contigCoverage"

    lengths_and_ids = pd.read_csv(lengths_file, delimiter="\t")
    lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"] #in what cases would we not have an EqualCoverageUnit?
    data_coverage = pd.read_csv(coverage_file, delimiter="\t")

    taxon_id_counts_per_unit = lengths_and_ids["ID"].value_counts() #Count of each taxonomic ID's occurrences
    taxon_id_counts_per_unit.sort_values(ascending=False, inplace=True)
    taxon_id_freq_per_unit = taxon_id_counts_per_unit / taxon_id_counts_per_unit.sum() #Frequency of each taxonomic ID's occurrences
    
    ids_to_plot = {} #Track id's we want to plot
    taxon_id_2_mapping_units = {} #Appears to be a table for converting a taxonomic id to its corresponding label. Possible that a single id may have multiple labels?
    
    #For all taxon ID counts, check if we are going to plot it
    for i in range(1, len(taxon_id_counts_per_unit)):
        
        current_id_label = taxon_id_counts_per_unit.index[i] #Label from ID column
        current_freq = taxon_id_freq_per_unit.iloc[i]
        
        if (current_freq >= min_plot_freq): #if the current unit's frequency is too low, we will not plot it
            
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label) #Get tax ID from label
            if len(matches) != 1:
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
            
            taxon_id = matches[0]
            ids_to_plot[taxon_id] = 1 #Mark id for plotting
            
            if (taxon_id not in taxon_id_2_mapping_units):
                taxon_id_2_mapping_units[taxon_id] = {}
            taxon_id_2_mapping_units[taxon_id][current_id_label] = 1 #One is arbitrary. We can plug in taxon_id and get current_id_label in the future
            
    taxon_coverage_indices = {} #track data_coverage indices in which a given taxon id occurs
    for index, entry in enumerate(data_coverage['taxonID']):
        if (entry not in taxon_coverage_indices):
            taxon_coverage_indices[entry] = []
        taxon_coverage_indices[entry].append(index)
    
    all_reads_lengths = [] #Track all lengths (integers)
    all_reads_identities = [] #Track all identities (floating point numbers)
    all_windows_coverages = [] #Track all windows coverages (floating point numbers)
    all_identity_densities = [] #Track all identity densities (floating point numbers)
    all_identitity_min_not_0 = [] #Track all identity min not 0 (floating point numbers)
    
    iter = 0
    for taxon_id in ids_to_plot.keys():  
        iter += 1
        if iter >= 2:
            break
        
        if not len(taxon_coverage_indices[int(taxon_id)]) > 0:
                print("Error: Could not find " + taxon_id + " in " + coverage_file)
                sys.exit(1)
        reads_count = 0
        
        for mapping_unit in taxon_id_2_mapping_units[taxon_id].keys(): #For each mapping unit; ex: 'kraken:taxid|595496|NC_012759.1'
            reads_count = reads_count + int(taxon_id_counts_per_unit[mapping_unit])
            coverage_mapping_unit_indices = [index for index, entry in enumerate(data_coverage['contigID']) if entry == mapping_unit ]
            
            if not len(coverage_mapping_unit_indices) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
            all_windows_coverages.append(data_coverage['readCoverage'][coverage_mapping_unit_indices])
            
            lengths_and_ids_mapping_unit = lengths_and_ids.loc[ lengths_and_ids['ID'] == mapping_unit ]
            
            if not len(lengths_and_ids_mapping_unit) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            all_reads_lengths = all_reads_lengths + list(lengths_and_ids_mapping_unit['Length'])
            all_reads_identities = all_reads_identities + list(lengths_and_ids_mapping_unit['Identity'])

        if not len(all_reads_lengths) == int(reads_count):
            print("Error: taxon_id counts do not match number of read lengths")
            sys.exit(1)
        
        #Preprocessing for plotting the 
        vector_identities = [0] * 101
        vector_identities = pd.Series(vector_identities)
        all_reads_identities_100 = []
        
        for ri in all_reads_identities:
            all_reads_identities_100.append(float(round(ri * 100) + 1.0))
            
        identities_table = pd.DataFrame(all_reads_identities_100)
        identities_table = identities_table / identities_table.sum()
        
        print(identities_table)
        
        for identity in identities_table.keys():
            vector_identities.iloc[int(identity)] = identities_table.iloc[int(identity)]
        
        all_identity_densities.append(vector_identities)
        all_identitity_min_not_0.append(min(vector_identities[vector_identities > 0]))
    
    iter = 0
    #do plotting
    pdf_output = PdfPages(file_prefix + ".identitiesAndCoverage.pdf")
    for taxon_id in ids_to_plot.keys():
        iter += 1
        if iter >= 2:
            break
           
        taxon_label = str(data_coverage["equalCoverageUnitLabel"][taxon_coverage_indices[int(taxon_id)][0]])
        
        fig = plt.figure(figsize=(12, 8))
        fig.text(0.5, 0.95, "MetaMaps mapping summary for " + taxon_label + " (taxon ID " + str(taxon_id) + ")"  + " - " + str(reads_count) + " + mapped reads assigned", ha='center', va='center')
        
        #generate read length histogram
        read_length_histogram_plot = plt.subplot2grid(loc=(0, 0), rowspan=1, colspan=1,  shape=(2, 3))
        read_length_histogram_plot.hist(all_reads_lengths, bins=100, color='blue', edgecolor='black', linewidth=1.2)
        read_length_histogram_plot.set_xlim(0, max(lengths_and_ids['Length']))
        read_length_histogram_plot.set_title("Read Length Histogram")
                            
        #generate identities bar plot
        read_identities_plot = plt.subplot2grid(loc=(0, 1), rowspan=1, colspan=1,  shape=(2, 3))
        read_identities_plot.set_title("Read Identities")
        read_identities_plot.bar(range(101), vector_identities.iloc[int(min(all_identitity_min_not_0)):int(len(vector_identities))], color='blue', edgecolor='black', linewidth=1.2)
        #read_identities_plot.set_ylim(0, max(all_identity_densities))
        
        #generate genome window coverage histogram
        genome_window_coverage_plot = plt.subplot2grid(loc=(0, 2), rowspan=1, colspan=1,  shape=(2, 3))
        genome_window_coverage_plot.set_title("Genome Window Coverage")
        genome_window_coverage_plot.hist(all_windows_coverages, bins=100, color='blue', edgecolor='black', linewidth=1.2)
        
        #generate plot for all genome window coverages
        genome_wide_coverage_over_all_contigs_plot = plt.subplot2grid(loc=(1, 0), rowspan=1, colspan=3,  shape=(2, 3))
        genome_wide_coverage_over_all_contigs_plot.set_title("Genome Wide Coverage Over All Contigs")
        
        indices_coverage_taxon_id = data_coverage['taxonID'].loc[data_coverage['taxonID'] == int(taxon_id)].index
        if (not len(indices_coverage_taxon_id) > 0):
            print("Error: Could not find " + taxon_id + " in " + coverage_file)
            sys.exit(1)
        
        all_windows_coverages = pd.Series()
        # all_windows_coverages_colors = []
        
        mapping_units = taxon_id_2_mapping_units[taxon_id].keys() #Get the list of all unique labels that match our taxonomic id
        mapping_units = pd.unique(data_coverage['contigID'].loc[data_coverage['taxonID'] == int(taxon_id)])
        
        reads_count = 0
        for mapping_unit in mapping_units:
            reads_count = reads_count + int(taxon_id_counts_per_unit[mapping_unit])
            if (not len(coverage_mapping_unit_indices) > 0):
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            coverage_mapping_unit_indices = data_coverage['contigID'].loc[data_coverage['contigID'] == mapping_unit].index #Get all indices for entries that match the current mapping unit
            current_coverages = data_coverage['readCoverage'][coverage_mapping_unit_indices]
            all_windows_coverages = pd.concat([all_windows_coverages, current_coverages], axis=0, ignore_index=True)
            
            # this_mapping_unit_color = "blue"
            # if ((mapping_unit_i % 2) == 0):
            #     this_mapping_unit_color = "red"
            #all_windows_coverages_colors = all_windows_coverages_colors + [this_mapping_unit_color] * len(coverages_this_mapping_unit)
        genome_wide_coverage_over_all_contigs_plot.scatter(list(range(len(all_windows_coverages))), all_windows_coverages, c='r', s=1)
   

        fig.text(0.5, 0.01, "Coordinate concatenated genome (1000s)", ha='center', va='bottom')
        pdf_output.savefig()
        plt.close()
    pdf_output.close()

if __name__ == "__main__":
    main()