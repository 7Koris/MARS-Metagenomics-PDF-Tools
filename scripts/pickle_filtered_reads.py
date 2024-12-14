import os
from os import path
import argparse
import read_data as rd
import pickle
import sys
import srs as srs    
    
def _get_coverage_outllier_list(file_prefix) -> list:
    mapping_units = rd.ReadDataMM(file_prefix, 0.0)
   
    mapping_units.load_coverage()
    mapping_units.filter_sig_bin_outliers(3, True)
    
    outliers = []
    for id in mapping_units.filtered_tax_ids.keys():
        if mapping_units.tax_id_is_outlier[id]:
            outliers.append(id)
            
    return outliers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Takes either MTSV or MetaMaps reads and loads into a pickle file for analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-mtsv", "--mtsv-file",  type=str, help="MTSV File")
    parser.add_argument("-mtsvl", "--mtsv-lookup-file",  type=str, help="MTSV Lookup File")
    parser.add_argument("-meta", "--meta-maps-file",  type=str, help="MetaMaps File")
    parser.add_argument("-metaref", "--meta-maps-reference-file", type=str, help="MetaMaps file to filter MTSV with")
    parser.add_argument("-mtsvref", "--mtsv-reference-file", type=str, help="MTSV file to filter MetaMaps reads with")
    parser.add_argument("-s", "--seed", type=int, help="Seed for random number generator", default=0)
    parser.add_argument("-a", "--alias", type=str, help="Alias for file", default=None)
    parser.add_argument("-C", "--clear", action="store_true", help="DANGER! Clears all files in reads directory") # DANGEROUS?
    parser.add_argument("-B", "--sig-bin", action="store_true", help="Enable sig-bin filtering (METAMAPS ONLY)", default=False)
    parser.add_argument("-r", "--rare", type=int, help="Rareify to given read count", default=None)
    parser.add_argument("-S", "--SRS", action="store_true", help="Enable SRS. Requires R to be set up in the current system's environment", default=False)
        
    args = parser.parse_args()
    config = vars(args)
    seed = config["seed"]
    
    mtsv_file = config["mtsv_file"]
    mtsv_lookup_file = config["mtsv_lookup_file"]
    meta_maps_file = config["meta_maps_file"]
    meta_maps_reference_file = config["meta_maps_reference_file"]
    mtsv_reference_file = config["mtsv_reference_file"]
    
    clear = config["clear"]
    sig_bin = config["sig_bin"]
    alias = config["alias"]
    rare = config["rare"]
    
    mtsv_present = (mtsv_file and mtsv_lookup_file)
    meta_present = (meta_maps_file)
    
    name = ""
    
    if mtsv_present and meta_present:
        sys.exit("Only one of MTSV or MetaMaps files can be present")
    
    if not mtsv_present and not meta_present:
        sys.exit("One of MTSV or MetaMaps files must be present")
    
    if not meta_maps_file and mtsv_file and not mtsv_lookup_file:
        sys.exit("MTSV lookup file must be present")
        
    read_data = rd.ReadData()
    
    if mtsv_file and mtsv_lookup_file:
        name = path.basename(mtsv_file)
        read_data.parse_mtsv_reads(mtsv_file, mtsv_lookup_file)
        read_data.resolve_lca()
        if meta_maps_reference_file:
            # Incidence filter
            read_data.parse_metamaps_reads_2_taxon(meta_maps_reference_file, True)
            read_data.prune_non_incidental_reads()
    elif meta_maps_file:
        name = path.basename(meta_maps_file)
        read_data.parse_metamaps_reads_2_taxon(meta_maps_file)
        if mtsv_reference_file:
            # Incidence filter
            read_data.parse_mtsv_reads_2_taxon(mtsv_reference_file, True)
            read_data.prune_non_incidental_reads()
    
    read_data.prune_by_level("genus") # Ensure all remaining reads can be resolved to a genus
    
    if alias is not None:
        name = alias
            
    # sig bin filter
    if meta_present and sig_bin:
        sig_bin_file = path.splitext(os.path.basename(meta_maps_file))[0]
        sig_bin_file = path.splitext(os.path.basename(sig_bin_file))[0]
        sig_bin_file_path = path.dirname(meta_maps_file)
        sig_bin_file = path.join(sig_bin_file_path, sig_bin_file)
        outlier_list = _get_coverage_outllier_list(sig_bin_file)
        outlier_dict = {}
        for outlier in outlier_list:
            outlier_dict[outlier] = True
        print("DISCARDING OUTLIERS")
        i = 0
        keys = list(read_data.reads.keys())
        for read in keys:
            currents_read = read_data.reads[read]
            if currents_read.assigned_taxon_id in outlier_dict:
                i += 1
                read_data.reads.pop(read)
        print(i, "outliers discarded")
    
    # rareify
    if rare:
        if not config["SRS"]:
            read_data.rarefy(rare, seed)
        else:
            read_data = srs.srs_on_read_data(read_data, rare, seed)
    
    # pickle reads
    currents_reads = read_data.reads
    processed_read_count = read_data.processed_read_count
    
    if clear:
        for file in os.listdir("../reads"):
            file_path = os.path.join("../reads", file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)
    
    if not os.path.exists("../reads"):
        os.makedirs("../reads")
     
    pickle.dump((currents_reads, processed_read_count), open("../reads/" + name + ".p", "wb"))
    print(len(read_data.reads), "reads pickled to reads/" + str(name) + ".p")
    