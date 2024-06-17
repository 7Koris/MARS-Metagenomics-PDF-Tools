import sys
import re
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import argparse


def plot_identities() -> None:
    print("")
    print("Python-MetaMaps Identity Plotter")
    print("\tThis script is based on plotIdentities_EM.R by the MetaMaps authors at https://github.com/DiltheyLab/MetaMaps")
    print("\tFor help, use -h or --help")
    print("")
    
    parser = argparse.ArgumentParser(description="Plot MetaMaps Identity Results", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-t", "--min-trim-mean", type=float, help="Minimum coverage trim mean value to consider a taxon ID as an outlier. Outliers will be written to a separate PDF where trim_mean <= min_trim_mean", default=0.0)
    parser.add_argument("-p", "--trim-proportion", type=float, help="Proportion of sorted coverage data to trim from both ends for outlier detection", default=0.03)
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    args = parser.parse_args()
    config = vars(args)
    
    start_time = time.time()
    t0 = time.time()
    
    file_prefix = config["classification_file_prefix"]
    min_plot_freq = config["min_frequency"]
    min_tm = config["min_trim_mean"]
    trim_proportion = config["trim_proportion"]
    verbose = config["verbose"]
    
    lengths_file = file_prefix + ".EM.lengthAndIdentitiesPerMappingUnit"
    coverage_file = file_prefix + ".EM.contigCoverage"

    lengths_and_ids = pd.read_csv(lengths_file, delimiter="\t")
    lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"] # In what cases would we not have an EqualCoverageUnit?
    data_coverage = pd.read_csv(coverage_file, delimiter="\t")
    
    if verbose:
        print("Data read took", time.time() - t0, "seconds to run")
    
    counts_per_mapping_unit = pd.Series(lengths_and_ids["ID"].value_counts().sort_values(ascending=False)) # Count of each taxonomic ID's occurrences
    taxon_id_freq_per_unit = counts_per_mapping_unit / counts_per_mapping_unit.sum() # Frequency of each taxonomic ID's occurrences

    ids_to_plot = {} # Track id's we want to plot
    taxon_id_2_mapping_units = {} # We can plug in taxon_id and get an mapping unit label ('kraken:taxid|595496|NC_012759.1') in the future. Possibly needed because a taxon_id may have multiple mapping units?
    
    t0 = time.time()
    # For all taxon ID counts, check if we are going to plot it
    for i in range(0, len(counts_per_mapping_unit)):
        current_id_label = counts_per_mapping_unit.index[i] # Label from ID column
        current_freq = taxon_id_freq_per_unit.iloc[i]
  
        if (current_freq >= min_plot_freq): # if the current unit's frequency is too low, we will not plot it
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label) # Get tax ID from label
            if len(matches) != 1: #Only one part of the label should contain id
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
            
            taxon_id = matches[0]
            ids_to_plot[taxon_id] = 1 #Mark id for plotting
            
            if (taxon_id not in taxon_id_2_mapping_units.keys()):
                taxon_id_2_mapping_units[taxon_id] = {}
            taxon_id_2_mapping_units[taxon_id][current_id_label] = 1 
    if verbose:
        print("Total number of ID labels:", len(counts_per_mapping_unit))
        print("Number of labels that will be plotted after filtering:", len(ids_to_plot))
    
    if verbose:
        print("ID filtering took", time.time() - t0, "seconds to run")
        t0 = time.time()
    
    taxon_coverage_indices = {} #track data_coverage indices in which a given taxon id occurs
    for index, entry in enumerate(data_coverage['taxonID']):
        if (entry not in taxon_coverage_indices.keys()):
            taxon_coverage_indices[entry] = []
        taxon_coverage_indices[entry].append(index)
    
    if verbose:
        print("Taxon coverage indices took", time.time() - t0, "seconds to run")

    min_id_x = None # for id plot
    max_id_x = None # for id plot
    max_id_y = None # for id plot
    plot_dict = {} # Dictionary that will store all data we process for plotting
    
    # Process data for plotting and store in plot_dict
    num_outliers = 0
    t0 = time.time()
    for taxon_id in ids_to_plot.keys():
        plot_dict[taxon_id] = {}
        
        taxon_label = str(data_coverage["equalCoverageUnitLabel"][taxon_coverage_indices[int(taxon_id)][0]])
        plot_dict[taxon_id]["taxon_label"] = taxon_label
        
        if not len(taxon_coverage_indices[int(taxon_id)]) > 0:
            print("Error: Could not find " + taxon_id + " in " + coverage_file)
            sys.exit(1)
        
        reads_count = 0
        all_reads_lengths = [] #Sequence lengths associated with an ID
        all_reads_identities = [] #Sequence identities associated with an ID
        all_windows_coverages = [] #Coverage values associated with an ID

        #For the given taxon_id, get all mapping units assocaited with it
        for mapping_unit in taxon_id_2_mapping_units[taxon_id].keys(): # For each mapping unit; ex: 'kraken:taxid|595496|NC_012759.1'
            reads_count = reads_count + int(counts_per_mapping_unit[mapping_unit])
            coverage_mapping_unit_indices = [index for index, entry in enumerate(data_coverage['contigID']) if entry == mapping_unit]
            
            if not len(coverage_mapping_unit_indices) > 0:
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            all_windows_coverages.extend(data_coverage['readCoverage'][coverage_mapping_unit_indices])
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
        plot_dict[taxon_id]["all_windows_coverages"] = all_windows_coverages
        plot_dict[taxon_id]["reads_count"] = reads_count
        
        all_reads_identities_100 = []
        for ri in all_reads_identities:
            new_id = int(round(ri * 100) + 1)
            all_reads_identities_100.append(new_id)
            
            if (max_id_x == None or new_id > max_id_x):
                max_id_x = new_id
            if (min_id_x == None or new_id < min_id_x):
                min_id_x = new_id
            
        
        plot_dict[taxon_id]["id_freq_table"] = pd.Series(all_reads_identities_100).value_counts(normalize=True).sort_index() #Frequencies of read identity counts
        if max_id_y == None or plot_dict[taxon_id]["id_freq_table"].max() > max_id_y:
            max_id_y = plot_dict[taxon_id]["id_freq_table"].max()

        # mapping_units = taxon_id_2_mapping_units[taxon_id].keys() 
        # If we want to filter mapping_units for this step, uncomment the previous line and comment out the following line
        mapping_units = data_coverage['contigID'].loc[data_coverage['taxonID'] == int(taxon_id)].unique()
        reads_count = 0
        iterator = 1
        all_windows_coverages = []
        all_colors = []
        for mapping_unit in mapping_units:
            reads_count = reads_count + int(counts_per_mapping_unit[mapping_unit])
            coverage_mapping_unit_indices = data_coverage['contigID'].loc[data_coverage['contigID'] == mapping_unit].index #Get all coverage indices for entries that match the current mapping unit
           
            if (not len(coverage_mapping_unit_indices) > 0):
                print("Error: Could not find " + mapping_unit + " in " + coverage_file)
                sys.exit(1)
                
            current_coverages = data_coverage['readCoverage'][coverage_mapping_unit_indices].values
            all_windows_coverages.extend(current_coverages)

            current_color = 'blue'
            if ((iterator % 2) == 0):
                current_color = "red"
            
            for i in range(len(current_coverages)):
                all_colors.append(current_color)
            iterator += 1
        plot_dict[taxon_id]["coverage_all_contigs"] = all_windows_coverages
        plot_dict[taxon_id]["all_colors"] = all_colors
        plot_dict[taxon_id]["genome_wide_reads_count"] = reads_count
        
        tm = stats.trim_mean(plot_dict[taxon_id]["coverage_all_contigs"], trim_proportion)
        if (tm <= min_tm):
            num_outliers += 1
            plot_dict[taxon_id]["outlier"] = True
        else:
            plot_dict[taxon_id]["outlier"] = False
        plot_dict[taxon_id]["trim_mean"] = tm
    
    print("Number of outliers:", num_outliers)
    if verbose:   
        print("Data processing took", time.time() - t0, "seconds to run")
    
    outliers_output = PdfPages(file_prefix + ".outliers.identitiesAndCoverage.pdf")
    pdf_output = PdfPages(file_prefix + ".identitiesAndCoverage.pdf")
    # Create plots
    t0 = time.time()
    
    id_list = list(ids_to_plot.keys())
    id_list.sort(key=lambda x: plot_dict[x]["trim_mean"], reverse=True)
    for key in id_list:
        print(key, plot_dict[key]["trim_mean"])
    
    for taxon_id in ids_to_plot.keys():
        fig = plt.figure(figsize=(12, 8))
        fig.text(0.01, 0.01, "MetaMaps mapping summary for " + plot_dict[taxon_id]["taxon_label"] + " (taxon ID " + str(taxon_id) + ")"  + " - " + str(plot_dict[taxon_id]["reads_count"]) + " mapped reads assigned", ha='left', va='center')
        
        # generate read length histogram
        read_length_histogram_plot = plt.subplot2grid(loc=(0, 0), rowspan=1, colspan=1,  shape=(2, 3))
        read_length_histogram_plot.hist(plot_dict[taxon_id]["all_reads_lengths"], bins='sturges', edgecolor='black', linewidth=1.2) #R bins (breaks) default is 30
        read_length_histogram_plot.set_xlim(0, lengths_and_ids['Length'].max())
        read_length_histogram_plot.set_title("Read Length Histogram", fontsize='small')
        read_length_histogram_plot.set_xlabel("Read Length", fontsize='small')
        
        # generate identities pie chart
        read_identities_plot = plt.subplot2grid(loc=(0, 1), rowspan=1, colspan=1,  shape=(2, 3))
        read_identities_plot.set_title("Read Identities", fontsize='small')
        read_identities_plot.pie(plot_dict[taxon_id]["id_freq_table"].values, plot_dict[taxon_id]["id_freq_table"].values, labels=plot_dict[taxon_id]["id_freq_table"].index, autopct='%1.1f%%', labeldistance=1)
                    
        # generate genome window coverage histogram
        genome_window_coverage_plot = plt.subplot2grid(loc=(0, 2), rowspan=1, colspan=1,  shape=(2, 3))
        genome_window_coverage_plot.set_title("Genome Window Coverage", fontsize='small')
        genome_window_coverage_plot.hist(plot_dict[taxon_id]["all_windows_coverages"], bins='sturges', edgecolor='black', linewidth=1.2)
        genome_window_coverage_plot.set_xlabel("Coverage", fontsize='small')
        
        # generate heat map of genome coverage
        genome_wide_coverage_over_all_contigs_plot = plt.subplot2grid(loc=(1, 0), rowspan=1, colspan=3,  shape=(2, 3))        
        genome_wide_coverage_over_all_contigs_plot.set_title("Genome-wide coverage over all contigs for " + plot_dict[taxon_id]["taxon_label"] + " (taxon ID " + str(taxon_id) + ")" + " - " + str(plot_dict[taxon_id]["genome_wide_reads_count"]) + " mapped reads assigned", fontsize='small')
        genome_wide_coverage_over_all_contigs_plot.imshow([plot_dict[taxon_id]["coverage_all_contigs"]], aspect='auto', cmap='hot', interpolation='bicubic')
        genome_wide_coverage_over_all_contigs_plot.set_xlabel("Coordinate concatenated genome (1000s)", fontsize='small')
        genome_wide_coverage_over_all_contigs_plot.get_yaxis().set_visible(False)
        
        fig.tight_layout()
        
        if plot_dict[taxon_id]["outlier"]:
            outliers_output.savefig()
        else:
            pdf_output.savefig()
            
        plt.close()
        
    if verbose:
        print("Plot generations took", time.time() - t0, "seconds to run")
        print("Total script took", time.time() - start_time, "seconds to run")
        
    print(str(pdf_output.get_pagecount()), "page(s) written to", pdf_output._file.fh.name)
    pdf_output.close()
    if outliers_output.get_pagecount() > 0:
        print(str(outliers_output.get_pagecount()), "outlier page(s) written to", outliers_output._file.fh.name)
        outliers_output.close()


if __name__ == "__main__":
    plot_identities()