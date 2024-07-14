import os
import argparse
import read_data as rd

def create_coverage_outllier_table(file_prefix, min_freq):
    mapping_units = rd.ReadDataMM(file_prefix, min_freq)
   
    mapping_units.load_coverage()
    #mapping_units.filter_coverage_tm_outliers(max_outlier_coverage, proportion, True)
    mapping_units.filter_sig_bin_outliers(3, True)
    
    outliers = []
    for id in mapping_units.filtered_tax_ids.keys():
        if mapping_units.tax_id_is_outlier[id]:
            outliers.append(id)
            
    iter = 0
    # Create text file with outliers
    with open(".ignoreids", "w") as f:
        for outlier in outliers:
            iter += 1
            f.write(str(outlier) + "\n")
    print(iter, "outliers written to .ignoreids")
    
    
def get_coverage_outllier_list(file_prefix, min_freq) -> list:
    mapping_units = rd.ReadDataMM(file_prefix, min_freq)
   
    mapping_units.load_coverage()
    mapping_units.filter_sig_bin_outliers(3, True)
    
    outliers = []
    for id in mapping_units.filtered_tax_ids.keys():
        if mapping_units.tax_id_is_outlier[id]:
            outliers.append(id)
            
    return outliers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Outlier detection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-t", "--trim-proportion", type=float, help="Proportion of coverages to trim out", default=0.003)
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    file_prefix = config["classification_file_prefix"]
    proportion = config["trim_proportion"]
    
    create_coverage_outllier_table(file_prefix, min_freq)
    
    