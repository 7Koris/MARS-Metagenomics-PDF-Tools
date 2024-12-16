import os
import argparse
import pickle
import sys
import read_data as rd
from os import path

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Takes MTSV reads and loads into a pickle file for analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-mtsv", "--mtsv-file",  type=str, help="MTSV File")
	parser.add_argument("-mtsvl", "--mtsv-lookup-file",  type=str, help="MTSV Lookup File")
	parser.add_argument("-p", "--prune-level", type=str, help="Prune reads to given level", default=None)
	parser.add_argument("-s", "--seed", type=int, help="Seed for random number generator", default=0)
	parser.add_argument("-a", "--alias", type=str, help="Alias for file", default=None)
	parser.add_argument("-r", "--rare", type=int, help="Rareify to given read count", default=None)
	parser.add_argument("-S", "--SRS", action="store_true", help="Enable SRS. Requires R to be set up in the current system's environment", default=False)
		
	args = parser.parse_args()
	config = vars(args)
	seed = config["seed"]
	
	mtsv_file = config["mtsv_file"]
	mtsv_lookup_file = config["mtsv_lookup_file"]
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
	
	mtsv_present = (mtsv_file and mtsv_lookup_file)
	if not mtsv_present:
		sys.exit("One of MTSV or MetaMaps files must be present")	
	
	name = ""
	read_data = rd.ReadData()
	if mtsv_file and mtsv_lookup_file:
		name = path.basename(mtsv_file)
		read_data.parse_mtsv_reads(mtsv_file, mtsv_lookup_file)
		read_data.resolve_lca()
	
	read_data.prune_by_level("genus") # Ensure all remaining reads can be resolved to a genus
	
	if alias is not None:
		name = alias
	
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