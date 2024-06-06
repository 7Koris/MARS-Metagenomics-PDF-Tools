import sys
import re
import pandas as pd
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
    lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"] # In what cases would we not have an EqualCoverageUnit?
    data_coverage = pd.read_csv(coverage_file, delimiter="\t")

    taxon_id_counts_per_unit = lengths_and_ids["ID"].value_counts() # Count of each taxonomic ID's occurrences
    taxon_id_counts_per_unit.sort_values(ascending=False, inplace=True)
    taxon_id_freq_per_unit = taxon_id_counts_per_unit / taxon_id_counts_per_unit.sum() # Frequency of each taxonomic ID's occurrences
    
    ids_to_plot = {} # Track id's we want to plot
    taxon_id_2_mapping_units = {} # We can plug in taxon_id and get an mapping unit label ('kraken:taxid|595496|NC_012759.1') in the future. Possibly needed because a taxon_id may have multiple mapping units?
    
    # For all taxon ID counts, check if we are going to plot it
    for i in range(1, len(taxon_id_counts_per_unit)):
        current_id_label = taxon_id_counts_per_unit.index[i] # Label from ID column
        current_freq = taxon_id_freq_per_unit.iloc[i]
        
        if (current_freq >= min_plot_freq): # if the current unit's frequency is too low, we will not plot it
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label) # Get tax ID from label
            if len(matches) != 1:
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
            
            taxon_id = matches[0]
            ids_to_plot[taxon_id] = 1 #Mark id for plotting
            
            if (taxon_id not in taxon_id_2_mapping_units.keys()):
                taxon_id_2_mapping_units[taxon_id] = {}
            taxon_id_2_mapping_units[taxon_id][current_id_label] = 1 
            
    taxon_coverage_indices = {} #track data_coverage indices in which a given taxon id occurs
    for index, entry in enumerate(data_coverage['taxonID']):
        if (entry not in taxon_coverage_indices.keys()):
            taxon_coverage_indices[entry] = []
        taxon_coverage_indices[entry].append(index)

    min_id_x = None #for id plot
    max_id_x = None #for id plot
    pdf_output = PdfPages(file_prefix + ".identitiesAndCoverage.pdf")
    plot_dict = {} #Dictionary that will store all data we process for plotting
    
    # Process data for plotting
    for taxon_id in ids_to_plot.keys():
        plot_dict[taxon_id] = {}
        
        taxon_label = str(data_coverage["equalCoverageUnitLabel"][taxon_coverage_indices[int(taxon_id)][0]])
        plot_dict[taxon_id]["taxon_label"] = taxon_label
        
        if not len(taxon_coverage_indices[int(taxon_id)]) > 0:
            print("Error: Could not find " + taxon_id + " in " + coverage_file)
            sys.exit(1)
        
        reads_count = 0
        all_reads_lengths = []
        all_reads_identities = []
        all_window_coverage_units = []

        #For the given taxon_id, get all mapping units assocaited with it
        for mapping_unit in taxon_id_2_mapping_units[taxon_id].keys(): # For each mapping unit; ex: 'kraken:taxid|595496|NC_012759.1'
            reads_count = reads_count + int(taxon_id_counts_per_unit[mapping_unit])
            coverage_mapping_unit_indices = [index for index, entry in enumerate(data_coverage['contigID']) if entry == mapping_unit]
            
            if not len(coverage_mapping_unit_indices) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            all_window_coverage_units.extend(data_coverage['readCoverage'][coverage_mapping_unit_indices])
            lengths_and_ids_mapping_unit = lengths_and_ids.loc[ lengths_and_ids['ID'] == mapping_unit ]
            
            if not len(lengths_and_ids_mapping_unit) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            all_reads_lengths.extend(lengths_and_ids_mapping_unit['Length'])
            all_reads_identities.extend(lengths_and_ids_mapping_unit['Identity'])
        
        if not len(all_reads_lengths) == int(reads_count):
            print("Error: taxon_id counts do not match number of read lengths")
            sys.exit(1)
        
        plot_dict[taxon_id]["all_reads_lengths"] = all_reads_lengths
        plot_dict[taxon_id]["all_reads_identities"] = all_reads_identities
        plot_dict[taxon_id]["all_window_coverage_units"] = all_window_coverage_units
        plot_dict[taxon_id]["reads_count"] = reads_count
        
        all_reads_identities_100 = []
        
        for ri in all_reads_identities:
            new_id = int(round(ri * 100) + 1)
            if (max_id_x == None or new_id > max_id_x):
                max_id_x = new_id
            if (min_id_x == None or new_id < min_id_x):
                min_id_x = new_id
            
            all_reads_identities_100.append(new_id)
        
        plot_dict[taxon_id]["id_freq_table"] = pd.Series(all_reads_identities_100).value_counts(normalize=True).sort_index() #Frequencies of read identity counts
    
    # Create plots
    for taxon_id in ids_to_plot.keys():
        fig = plt.figure(figsize=(12, 8))
        fig.text(0.5, 0.95, "MetaMaps mapping summary for " + plot_dict[taxon_id]["taxon_label"] + " (taxon ID " + str(taxon_id) + ")"  + " - " + str(plot_dict[taxon_id]["reads_count"]) + " mapped reads assigned", ha='center', va='center')
        
        #generate read length histogram
        read_length_histogram_plot = plt.subplot2grid(loc=(0, 0), rowspan=1, colspan=1,  shape=(2, 3))
        read_length_histogram_plot.hist(plot_dict[taxon_id]["all_reads_lengths"], bins='sturges', edgecolor='black', linewidth=1.2) #R bins (breaks) default is 30
        read_length_histogram_plot.set_xlim(0, lengths_and_ids['Length'].max())
        read_length_histogram_plot.set_title("Read Length Histogram", fontsize='small')
        read_length_histogram_plot.set_xlabel("Read Length", fontsize='small')
        
        #generate identities bar plot
        read_identities_plot = plt.subplot2grid(loc=(0, 1), rowspan=1, colspan=1,  shape=(2, 3))
        read_identities_plot.set_title("Read Identities", fontsize='small')
        read_identities_plot.bar(plot_dict[taxon_id]["id_freq_table"].keys(), plot_dict[taxon_id]["id_freq_table"].values, color='blue', edgecolor='black', linewidth=1.2)
        read_identities_plot.set_xlim(min_id_x, max_id_x)
        read_identities_plot.set_xlabel("Identity", fontsize='small')
        
        #read_identities_plot.set_ylim(0, identities_table.max())
        read_identities_plot.set_ylim(0, 1) #In the R script the y lim is set to the table max but the example PDF output does not reflect this. I will set it to 1 for now.
            
        #generate genome window coverage histogram
        genome_window_coverage_plot = plt.subplot2grid(loc=(0, 2), rowspan=1, colspan=1,  shape=(2, 3))
        genome_window_coverage_plot.set_title("Genome Window Coverage", fontsize='small')
        genome_window_coverage_plot.hist(plot_dict[taxon_id]["all_window_coverage_units"], bins='sturges', edgecolor='black', linewidth=1.2)
        genome_window_coverage_plot.set_xlabel("Coverage", fontsize='small')
        
        # generate plot for all genome window coverages
        genome_wide_coverage_over_all_contigs_plot = plt.subplot2grid(loc=(1, 0), rowspan=1, colspan=3,  shape=(2, 3))
        genome_wide_coverage_over_all_contigs_plot.set_title("Genome Wide Coverage Over All Contigs", fontsize='small')
        
        indices_coverage_taxon_id = data_coverage['taxonID'].loc[data_coverage['taxonID'] == int(taxon_id)].index
        if (not len(indices_coverage_taxon_id) > 0):
            print("Error: Could not find " + taxon_id + " in " + coverage_file)
            sys.exit(1)
         
        mapping_units = taxon_id_2_mapping_units[taxon_id].keys() #Get the list of all unique labels that match our taxonomic id
        
        reads_count = 0
        iterator = 1
        all_windows_coverages = []
        all_colors = []
        for mapping_unit in mapping_units:
            reads_count = reads_count + int(taxon_id_counts_per_unit[mapping_unit])
            if (not len(coverage_mapping_unit_indices) > 0):
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
            
            coverage_mapping_unit_indices = data_coverage['contigID'].loc[data_coverage['contigID'] == mapping_unit].index #Get all coverage indices for entries that match the current mapping unit
            current_coverages = data_coverage['readCoverage'][coverage_mapping_unit_indices].values
            all_windows_coverages.extend(current_coverages)

            current_color = 'blue'
            if ((iterator % 2) == 0):
                current_color = "red"
            
            for i in range(len(current_coverages)):
                all_colors.append(current_color)
            iterator += 1
            
        genome_wide_coverage_over_all_contigs_plot.scatter(list(range(len(all_windows_coverages))), all_windows_coverages, s=1, c=all_colors)
        genome_wide_coverage_over_all_contigs_plot.set_xlabel("Coordinate concatenated genome (1000s)", fontsize='small')
        # Cleanup for next plot
        pdf_output.savefig()
        plt.close()
    pdf_output.close()

if __name__ == "__main__":
    main()