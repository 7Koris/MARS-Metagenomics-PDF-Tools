import argparse
import os
import pickle
import sys
import pandas as pd
import read_data as rd
from plotnine import *
from os import listdir, path
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats as scistats

def generate_diverging_sk(output_name, all_read_data_names, all_read_data):
	if not os.path.exists("reads"):
		os.makedirs("reads")
	
	genus_reads = {}
	
	for idx, read_data in enumerate(all_read_data):
		current_name = all_read_data_names[idx]
		tax_count_dict = read_data.get_tax_count_dict()
		genus_level_counts = tax_count_dict["phylum"]
		genus_reads[current_name] = genus_level_counts
		
	for gr_name in genus_reads:
		for read in genus_reads[gr_name]:
			for gr_name2 in genus_reads:
				if gr_name == gr_name2:
					continue
				if read not in genus_reads[gr_name2]:
					genus_reads[gr_name2][read] = 0
	
	for gr_name in genus_reads:
		genus_reads[gr_name] = dict(sorted(genus_reads[gr_name].items()))
 
	for gr_name in genus_reads:
		reads = genus_reads[gr_name]
	
	vals = []
	for gr_name in genus_reads:
		reads = genus_reads[gr_name]
		vals.append(list(reads.values()))
		
	z_df = pd.DataFrame(genus_reads)
	z_df = z_df.transpose()
	z_df = z_df.apply(scistats.zscore)
	z_df = z_df.transpose()
	
	gg_plots = []
 
	for gr_name in genus_reads:
		df_subset = z_df[gr_name]
		df_subset = df_subset.sort_values(key=abs)
		df_subset = df_subset.reset_index()
		df_subset.columns = ["names", gr_name]
		df_subset = df_subset.dropna()

		gg = ggplot(df_subset, aes(x="reorder(names, " + gr_name + ")", y=gr_name)) + \
			geom_bar(stat='identity', width=.5) + \
			coord_flip() + \
			labs(title=f"Z-score of {gr_name}", x="Genus", y="Z-score")
		gg_plots.append(gg)

	with PdfPages(output_name + '.diverging.pdf') as pdf:
		for plot in gg_plots:
			pdf.savefig(plot.draw())
	   

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Diverging Phylum", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", "--files", nargs='+', help="Files to compare")
 
	args = parser.parse_args()
	config = vars(args)
	files = config["files"]
	
	# if no files given
	if files == None:
		files = []
		dir_files = listdir("../reads")
		for file in dir_files:
			if path.splitext(file)[1] != ".p":
				continue
			files.append(path.splitext(file)[0] + ".p")
	   
		if file is None:
			sys.exit("No files found in reads ~/reads directory.") 
	
	output_prefix = ""
	all_read_data = []
	all_read_data_names = []
	
	for file in files:
		current_name = path.splitext(file)[0]
		output_prefix += current_name + "_"
		read_data = rd.ReadData()
		read_data.load_data(pickle.load(open("../reads/" + file, "rb")))
		all_read_data.append(read_data)
		all_read_data_names.append(current_name)

	generate_diverging_sk( output_prefix, all_read_data_names, all_read_data)
	