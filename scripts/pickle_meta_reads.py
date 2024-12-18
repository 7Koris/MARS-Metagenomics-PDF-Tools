import os
import argparse
import pickle
import sys
import read_data as rd
import utility as ut
from os import path
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Takes MetaMaps reads and loads into a pickle file for analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-meta", "--meta-maps-file",  type=str, help="MetaMaps Reads2Taxon File")
	parser.add_argument("-p", "--prune-level", type=str, help="Prune reads to given level", default=None)
	parser.add_argument("-s", "--seed", type=int, help="Seed for random number generator", default=0)
	parser.add_argument("-a", "--alias", type=str, help="Alias for file", default=None)
	parser.add_argument("-B", "--sig-bin", action="store_true", help="Enable sig-bin filtering (METAMAPS ONLY)", default=False)
	parser.add_argument("-r", "--rare", type=int, help="Rareify to given read count", default=None)
	parser.add_argument("-S", "--SRS", action="store_true", help="Enable SRS. Requires R to be set up in the current system's environment", default=False)
		
	args = parser.parse_args()
	config = vars(args)
	seed = config["seed"]
	meta_maps_file = config["meta_maps_file"]
	sig_bin = config["sig_bin"]
	alias = config["alias"]
	rare = config["rare"]
	prune_level = config["prune_level"]
	
	valid_ranks = {
		"superkingdom": 1,
		"family": 1,
		"phylum": 1,
		"class": 1,
		"order": 1,	
		"genus": 1,
		"kingdom": 1,
		"species": 1,
	}
	
	if valid_ranks[prune_level.lower()] != 1:
		sys.exit("Invalid rank")
		
	if (config["SRS"]):
		import srs as srs
	
	meta_present = (meta_maps_file)
	if not meta_present:
		sys.exit("MetaMaps file must be present")

	name = ""	
	read_data = rd.ReadData()
	
	if meta_maps_file:
		name = path.basename(meta_maps_file)
		read_data.parse_metamaps_reads_2_taxon(meta_maps_file)
	
	# TODO: Check that level is valid rank
	if prune_level is not None:
		read_data.prune_by_level(prune_level.lower()) # Ensure all remaining reads can be resolved to a genus
		
	if alias is not None:
		name = alias
			
	# sig bin filter
	if meta_present and sig_bin:
		sig_bin_file = path.splitext(os.path.basename(meta_maps_file))[0]
		sig_bin_file = path.splitext(os.path.basename(sig_bin_file))[0]
		sig_bin_file_path = path.dirname(meta_maps_file)
		sig_bin_file = path.join(sig_bin_file_path, sig_bin_file)
		outlier_list = ut.get_coverage_outllier_list(sig_bin_file)
		outlier_dict = {}
		for outlier in outlier_list:
			outlier_dict[outlier] = True
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
	
	if not os.path.exists("../reads"):
		os.makedirs("../reads")
	 
	pickle.dump((currents_reads, processed_read_count), open("../reads/" + name + ".p", "wb"))
	print(len(read_data.reads), "reads pickled to reads/" + str(name) + ".p")
	