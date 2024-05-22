import sys
import re
import pandas as pd
import matplotlib as plt
from matplotlib.backends.backend_pdf import PdfPages


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: Main.py <file prefix> <min plot frequency>")
        sys.exit(1)
        
    file_prefix = sys.argv[1]
    min_plot_freq = float(sys.argv[2])

    length_and_ids_per_mapping_unit = pd.read_csv(file_prefix + ".EM.lengthAndIdentitiesPerMappingUnit", delimiter="\t")
    length_and_ids_per_mapping_unit = length_and_ids_per_mapping_unit[length_and_ids_per_mapping_unit["AnalysisLevel"] == "EqualCoverageUnit"] #in what cases would we not have an EqualCoverageUnit?
    data_coverage = pd.read_csv(file_prefix + ".EM.contigCoverage", delimiter="\t")

    taxon_id_counts_per_unit = length_and_ids_per_mapping_unit["ID"].value_counts() #Count of each taxonomic ID's occurrences
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
            
    plot_labels = []
    identities_per_label = []
    lengths_per_label = []
  
    pdf_output = PdfPages(file_prefix + ".identitiesAndCoverage.pdf")
    pdf_output_coverage_plots = PdfPages(file_prefix + ".coveragePerContig_MetaMapComplete.pdf")
    #Todo: Configure pdf
    #See: https://matplotlib.org/stable/gallery/misc/multipage_pdf.html
    # par(mar = c(5, 3, 5, 3), oma = c(0,0,2,0))
    # m <- rbind(c(1,2,3), c(4, 4, 4))
    # layout(m)
    
    all_identity_densities = [] #Not sure what this is for yet
    all_identitity_min_not_0 = [] #Not sure what this is for yet
    
    #for 
    


if __name__ == "__main__":
    main()