import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import argparse
import MappingUnitData as mu


def plot_identities(file_prefix, min_plot_freq, min_tm: float=0.0, trim_proportion: float=.003, output: str="", excluded_tax_ids: list[str]=[], skip_coverage_filter: bool=False, verbose: bool=False) -> None:
    start_time = time.time()
    t0 = time.time()
    
    mu_data = mu.MappingUnitData(file_prefix, excluded_taxon_ids=excluded_tax_ids)
    mu_data.filter_by_frequency(min_plot_freq)
    mu_data.load_coverage()
    if not skip_coverage_filter:
        mu_data.filter_coverage_tm_outliers(min_tm, trim_proportion, True)
    
    if verbose:
        print("Data read and filtering took", time.time() - t0, "seconds to run")
    
    min_id_x = None # for id plot
    max_id_x = None # for id plot
    max_id_y = None # for id plot
    max_read_length = None
    
    plot_dict = {}
    
    for idx, level, unit, read_i, identitiy_score, length in mu_data.mapping_units.itertuples():
        if unit in mu_data.mapping_unit_2_tax_id: # if unit was not filtered
            id = mu_data.mapping_unit_2_tax_id[unit]
            id = str(id)
            
            if id not in plot_dict:
                plot_dict[id] = {}
            if "lengths" not in plot_dict[id]:
                plot_dict[id]["lengths"] = []
            plot_dict[id]["lengths"].append(length)
            
            if max_read_length is None or length > max_read_length:
                max_read_length = length
                
            identitiy_score = int(round(identitiy_score * 100) + 1)
            
            if (max_id_x == None or identitiy_score > max_id_x):
                max_id_x = identitiy_score
            if (min_id_x == None or identitiy_score < min_id_x):
                min_id_x = identitiy_score
              
            if "identitiy_scores" not in plot_dict[id]:
                plot_dict[id]["identitiy_scores"] = []
            plot_dict[id]["identitiy_scores"].append(identitiy_score)

    for id in plot_dict.keys():
        if "all_reads_count" not in plot_dict[id]:
            plot_dict[id]["all_reads_count"] = 0
        if "reads_count" not in plot_dict[id]:
            plot_dict[id]["reads_count"] = 0
        
        plot_dict[id]["id_freq_table"] = pd.Series(plot_dict[id]["identitiy_scores"]).value_counts(normalize=True).sort_index()
        if max_id_y == None or plot_dict[id]["id_freq_table"].max() > max_id_y:
            max_id_y = plot_dict[id]["id_freq_table"].max()
        
        for contig in mu_data.tax_id_2_filtered_contigs[id].keys():
            if "coverages" not in plot_dict[id]:
                plot_dict[id]["coverages"] = []
            plot_dict[id]["coverages"].extend(mu_data.contig_2_coverages[contig])
        
        iter = 1
        for contig in mu_data.tax_id_2_all_contigs[id].keys():
            plot_dict[id]["all_reads_count"] += mu_data.counts_per_unit[contig]
            current_color = 'blue'
            if ((iter % 2) == 0):
                current_color = "red"
            iter += 1
                
            if "all_coverages" not in plot_dict[id]:
                plot_dict[id]["all_coverages"] = []
            if "all_colors" not in plot_dict[id]:
                plot_dict[id]["all_colors"] = []
            plot_dict[id]["all_coverages"].extend(mu_data.contig_2_coverages[contig])
            plot_dict[id]["all_colors"].extend([current_color] * len(mu_data.contig_2_coverages[contig]))
        
        for unit in mu_data.tax_id_2_mapping_units[id].keys():
            plot_dict[id]["reads_count"] += mu_data.counts_per_unit[unit] 
    
    if output == "":
        pdf_output = PdfPages(file_prefix + ".plots.pdf")
        outliers_output = PdfPages(file_prefix + ".plots.outliers.pdf")
        all_output = PdfPages(file_prefix + ".plots.all.pdf")
    else:
        pdf_output = PdfPages(output + "plots.pdf")
        outliers_output = PdfPages(output + ".outliers.pdf")
        all_output = PdfPages(output + ".all.pdf")
    
    for taxon_id in mu_data.filtered_tax_ids.keys():
        fig = plt.figure(figsize=(12, 8))
        fig.text(0.01, 0.01, "MetaMaps mapping summary for " + mu_data.tax_id_2_name[taxon_id] + " (taxon ID " + str(taxon_id) + ")"  + " - " + str(plot_dict[taxon_id]["reads_count"]) + " mapped reads assigned", ha='left', va='center')
        
        # generate read length histogram
        read_length_histogram_plot = plt.subplot2grid(loc=(0, 0), rowspan=1, colspan=1,  shape=(2, 3))
        read_length_histogram_plot.hist(plot_dict[taxon_id]["lengths"], bins='sturges', edgecolor='black', linewidth=1.2) #R bins (breaks) default is 30
        read_length_histogram_plot.set_xlim(0, max_read_length)
        read_length_histogram_plot.set_title("Read Length Histogram", fontsize='small')
        read_length_histogram_plot.set_xlabel("Read Length", fontsize='small')
        
        # generate identities bar plot
        read_identities_plot = plt.subplot2grid(loc=(0, 1), rowspan=1, colspan=1,  shape=(2, 3))
        read_identities_plot.set_title("Read Identities", fontsize='small')
        read_identities_plot.bar(plot_dict[taxon_id]["id_freq_table"].keys(), plot_dict[taxon_id]["id_freq_table"].values, color='blue', edgecolor='black', linewidth=1.2)
        read_identities_plot.set_xlim(min_id_x, max_id_x)
        read_identities_plot.set_xlabel("Identity", fontsize='small')
        read_identities_plot.set_ylim(0, max_id_y)
                    
        # generate genome window coverage histogram
        genome_window_coverage_plot = plt.subplot2grid(loc=(0, 2), rowspan=1, colspan=1,  shape=(2, 3))
        genome_window_coverage_plot.set_title("Genome Window Coverage", fontsize='small')
        genome_window_coverage_plot.hist(plot_dict[taxon_id]["coverages"], bins='sturges', edgecolor='black', linewidth=1.2)
        genome_window_coverage_plot.set_xlabel("Coverage", fontsize='small')
        
        # generate scatter plot for all genome window coverages
        genome_wide_coverage_over_all_contigs_plot = plt.subplot2grid(loc=(1, 0), rowspan=1, colspan=3,  shape=(2, 3))        
        genome_wide_coverage_over_all_contigs_plot.set_title("Genome-wide coverage over all contigs for " + mu_data.tax_id_2_name[taxon_id] + " (taxon ID " + str(taxon_id) + ")" + " - " + str(plot_dict[taxon_id]["all_reads_count"]) + " mapped reads assigned", fontsize='small')
        genome_wide_coverage_over_all_contigs_plot.scatter(list(range(len(plot_dict[taxon_id]["all_coverages"]))), plot_dict[taxon_id]["all_coverages"], s=1, c=plot_dict[taxon_id]["all_colors"])
        genome_wide_coverage_over_all_contigs_plot.set_xlabel("Coordinate concatenated genome (1000s)", fontsize='small')
        
        fig.tight_layout()
        
        if not skip_coverage_filter:
            if mu_data.tax_id_is_outlier[taxon_id]:
                outliers_output.savefig()
            else:
                pdf_output.savefig()
            all_output.savefig()
        else:
            pdf_output.savefig()
            
        plt.close()
        
    if verbose:
        print("Total script took", time.time() - start_time, "seconds to run")
        
    print(str(pdf_output.get_pagecount()), "page(s) written to", pdf_output._file.fh.name)
    pdf_output.close()
    if outliers_output.get_pagecount() > 0:
        print(str(outliers_output.get_pagecount()), "outlier page(s) written to", outliers_output._file.fh.name)
        print(str(all_output.get_pagecount()), "page(s) written to", all_output._file.fh.name)
        outliers_output.close()
        all_output.close()


if __name__ == "__main__":
    print("")
    print("Python-MetaMaps Identity Plotter")
    print("\tThis script is based on plotIdentities_EM.R by the MetaMaps authors at https://github.com/DiltheyLab/MetaMaps")
    print("\tFor help, use -h or --help")
    print("")
    
    parser = argparse.ArgumentParser(description="Plot MetaMaps Identity Results", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-t", "--min-trim-mean", type=float, help="Minimum coverage trim mean value to consider a taxon ID as an outlier. Outliers will be written to a separate PDF where trim_mean <= min_trim_mean", default=0.0)
    parser.add_argument("-p", "--trim-proportion", type=float, help="Proportion of sorted coverage data to trim from both ends for outlier detection", default=0.003)
    parser.add_argument("-o", "--output", type=str, help="Output file name", default="")
    parser.add_argument("-I", "--ignore-ids", action="store_true", help="Ignore ids in the file .ignoreids")
    parser.add_argument("-S", "--skip-coverage-filter", action="store_true", help="Skip the coverage filtering step. Saves time if you already have a .ignoreids file. Will not generate an outlier pdf.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    args = parser.parse_args()
    config = vars(args)
    
    file_prefix = config["classification_file_prefix"]
    min_plot_freq = config["min_frequency"]
    min_tm = config["min_trim_mean"]
    trim_proportion = config["trim_proportion"]
    output = config["output"]
    ignore_file = config["ignore_ids"]
    skip_coverage_filter = config["skip_coverage_filter"]
    verbose = config["verbose"]
    
    excluded_tax_ids = []
    if ignore_file:
        for id in open(".ignoreids", "r"):
            id = id.strip()
            if id not in excluded_tax_ids:
                excluded_tax_ids.append(id)
    
    plot_identities(file_prefix, min_plot_freq, min_tm, trim_proportion, output, excluded_tax_ids, skip_coverage_filter, verbose)