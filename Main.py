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
    lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"] #in what cases would we not have an EqualCoverageUnit?
    data_coverage = pd.read_csv(coverage_file, delimiter="\t")

    taxon_id_counts_per_unit = lengths_and_ids["ID"].value_counts() #Count of each taxonomic ID's occurrences
    taxon_id_counts_per_unit.sort_values(ascending=False, inplace=True)
    taxon_id_freq_per_unit = taxon_id_counts_per_unit / taxon_id_counts_per_unit.sum() #Frequency of each taxonomic ID's occurrences
    
    ids_to_plot = {} #Track id's we want to plot
    taxon_id_2_mapping_units = {} #Not sure what this is for yet
    
    #For all taxon ID counts, check if we are going to plot it
    for i in range(1, len(taxon_id_counts_per_unit)):
        
        current_id_label = taxon_id_counts_per_unit.index[i] #Label from ID column
        current_freq = taxon_id_freq_per_unit.iloc[i]
        
        if (current_freq >= min_plot_freq): #if the current unit's frequency is too low, we will not plot it
            
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label) #Get tax ID from label
            if len(matches) != 1: #why is this == 2 in the original R script?
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
            
            taxon_id = matches[0]
            ids_to_plot[taxon_id] = 1 #Mark id for plotting
            
            if (taxon_id not in taxon_id_2_mapping_units):
                taxon_id_2_mapping_units[taxon_id] = {}
            taxon_id_2_mapping_units[taxon_id][current_id_label] = 1 #??
            
    
    
    all_identity_densities = [] #Not sure what this is for yet
    all_identitity_min_not_0 = [] #Not sure what this is for yet
    
    taxon_coverage_indices = {} #track data_coverage indices in which a given taxon id occurs
    for index, entry in enumerate(data_coverage['taxonID']):
        if (entry not in taxon_coverage_indices):
            taxon_coverage_indices[entry] = []
        taxon_coverage_indices[entry].append(index)
    
    all_reads_lengths = [] #Track all lengths (integers)
    all_reads_identities = [] #Track all identities (floating point numbers)
    all_windows_coverages = []
    
    
    for taxon_id in ids_to_plot.keys():  
        if not len(taxon_coverage_indices[int(taxon_id)]) > 0:
                print("Error: Could not find " + taxon_id + " in " + coverage_file)
                sys.exit(1)
        reads_count = 0
        for mapping_unit in taxon_id_2_mapping_units[taxon_id].keys(): #For each mapping unit; ex: 'kraken:taxid|595496|NC_012759.1'
            reads_count = reads_count + int(taxon_id_counts_per_unit[mapping_unit])
            coverage_indices_mapping_unit = [index for index, entry in enumerate(data_coverage['contigID']) if entry == mapping_unit ]
            
            if not len(coverage_indices_mapping_unit) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
            all_windows_coverages.append(data_coverage['readCoverage'][coverage_indices_mapping_unit])
            
            lengths_and_ids_mapping_unit = lengths_and_ids.loc[ lengths_and_ids['ID'] == mapping_unit ]
            
            if not len(lengths_and_ids_mapping_unit) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            all_reads_lengths = all_reads_lengths + list(lengths_and_ids_mapping_unit['Length'])
            all_reads_identities = all_reads_identities + list(lengths_and_ids_mapping_unit['Identity'])

        if not len(all_reads_lengths) == int(reads_count):
            print("Error: taxon_id counts do not match number of read lengths")
            sys.exit(1)
        
    #do plotting
    pdf_output = PdfPages(file_prefix + ".identitiesAndCoverage.pdf")
    for taxon_id in ids_to_plot.keys(): #Iterate through all ids that we will plot         
        taxon_label = str(data_coverage["equalCoverageUnitLabel"][taxon_coverage_indices[int(taxon_id)][0]])
        
        fig = plt.figure(figsize=(12, 8))
        fig.text(0.5, 0.95, "MetaMaps mapping summary for " + taxon_label + " (taxon ID " + str(taxon_id) + ")"  + " - " + str(reads_count) + " + mapped reads assigned", ha='center', va='center')
        
        #generate read length histogram
        read_length_histogram_plot = plt.subplot2grid(loc=(0, 0), rowspan=1, colspan=1,  shape=(2, 3))
        read_length_histogram_plot.hist(all_reads_lengths, bins=100, color='blue', edgecolor='black', linewidth=1.2)
        read_length_histogram_plot.set_xlim(0, max(lengths_and_ids['Length']))
        read_length_histogram_plot.set_title("Read Length Histogram")
            
        #generate identities bar plot
        # read_identities_plot = plt.subplot2grid(loc=(0, 1), rowspan=1, colspan=1,  shape=(2, 3))
        # read_identities_plot.set_title("Read Identities")
        
        #generate genome window coverage histogram
        # genome_window_coverage_plot = plt.subplot2grid(loc=(0, 2), rowspan=1, colspan=1,  shape=(2, 3))
        # genome_window_coverage_plot.set_title("Genome Window Coverage")
        
        #generate plot for all genome window coverages
        # genome_wide_coverage_over_all_contigs_plot = plt.subplot2grid(loc=(1, 0), rowspan=1, colspan=3,  shape=(2, 3))
        # genome_wide_coverage_over_all_contigs_plot.set_title("Genome Wide Coverage Over All Contigs")

        fig.text(0.5, 0.01, "Coordinate concatenated genome (1000s)", ha='center', va='bottom')
        pdf_output.savefig()
        plt.close()
    pdf_output.close()
                
                
                
                
                
                
                
            
    
        
            
            
            
            
    


if __name__ == "__main__":
    main()