import os
import argparse
import MappingUnitData as mu

def create_coverage_outllier_table(file_prefix, min_freq, max_outlier_coverage, proportion):
    mapping_units = mu.MappingUnitData(file_prefix, min_freq)
   
    mapping_units.load_coverage()
    #mapping_units.filter_coverage_tm_outliers(max_outlier_coverage, proportion, True)
    mapping_units.filter_z_score_outliers(50, True)
    
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Outlier detection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    parser.add_argument("-m", "--max-avg-outlier-coverage", type=float, help="Maximum average outlier coverage", default=0)
    parser.add_argument("-t", "--trim-proportion", type=float, help="Proportion of coverages to trim out", default=0.003)
    parser.add_argument("-i", "--ignore-coverage", action='store_true', help="Ignore the coverage filtering step")
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    file_prefix = config["classification_file_prefix"]
    max_outlier_coverage = config["max_avg_outlier_coverage"]
    proportion = config["trim_proportion"]
    
    create_coverage_outllier_table(file_prefix, min_freq, max_outlier_coverage, proportion)
    
    