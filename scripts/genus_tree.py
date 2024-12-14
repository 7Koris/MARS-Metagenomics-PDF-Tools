import argparse
import sys
import read_data as rd
import alphas as alphas
import pickle
from os import path, listdir
from plotnine import *
import plotly.express as px
from ete3 import NCBITaxa

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Alpha Diversity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", "--file", type=str, help="Pickle file name", default=None)
	args = parser.parse_args()
	config = vars(args)
	file = config["file"]
	
	if file is None:
		files = listdir("reads")
		for file in files:
			if path.splitext(file)[1] != ".p":
				continue
			file = path.splitext(file)[0] + ".p"
			break
		if file is None:
			sys.exit("No files found in ~/reads directory.") 
			
	data = pickle.load(open("../reads/" + file, "rb"))
	read_data = rd.ReadData()
	read_data.load_data(data)
 
	ncbi = NCBITaxa()
	id_dict = {}
	for read in read_data.reads:
		assignment = read_data.reads[read].get_assignment()
		if assignment not in id_dict:
			id_dict[assignment] = 1
		else:
			id_dict[assignment] += 1
	# SK
	# Phylum
	# Class
	# Order
	# Family
	# Genus
	rank_dict = {}
	for value in id_dict:
		lineage = ncbi.get_lineage(value)
		names = ncbi.get_taxid_translator(lineage)
		ranks = ncbi.get_rank(lineage)
  
		skrank = -1
		prank = -1
		crank = -1
		orank = -1
		frank = -1
		grank = -1


		for idx, rank in ranks.items():
			if rank == "superkingdom":
				skrank = idx
			elif rank == "phylum":
				prank = idx
			elif rank == "class":
				crank = idx
			elif rank == "order":
				orank = idx	
			elif rank == "family":
				frank = idx
			elif rank == "genus":
				grank = idx
		if (skrank == -1 or prank == -1 or crank == -1 or orank == -1 or frank == -1 or grank == -1):
			continue
		else:
			if names[skrank] not in rank_dict:
				rank_dict[names[skrank]] = {}
			if names[prank] not in rank_dict[names[skrank]]:
				rank_dict[names[skrank]][names[prank]] = {}
			if names[crank] not in rank_dict[names[skrank]][names[prank]]:
				rank_dict[names[skrank]][names[prank]][names[crank]] = {}
			if names[orank] not in rank_dict[names[skrank]][names[prank]][names[crank]]:
				rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]] = {}
			if names[frank] not in rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]]:
				rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]][names[frank]] = {}
			if names[grank] not in rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]][names[frank]]:
				rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]][names[frank]][names[grank]] = 0
			rank_dict[names[skrank]][names[prank]][names[crank]][names[orank]][names[frank]][names[grank]] += 1

	names = []
	parents = []
	values = []
							
	for superkingdom in rank_dict:
		names.append(superkingdom)
		parents.append("")
  
		d_sum = 0
		for phylum in rank_dict[superkingdom]:
			for class_ in rank_dict[superkingdom][phylum]:
				for order in rank_dict[superkingdom][phylum][class_]:
					for family in rank_dict[superkingdom][phylum][class_][order]:
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							d_sum += rank_dict[superkingdom][phylum][class_][order][family][genus]
		values.append(0)
  
	for superkingdom in rank_dict:
		for phylum in rank_dict[superkingdom]:
			names.append(phylum)
			parents.append(superkingdom)
	  
			d_sum = 0
			for class_ in rank_dict[superkingdom][phylum]:
				for order in rank_dict[superkingdom][phylum][class_]:
					for family in rank_dict[superkingdom][phylum][class_][order]:
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							d_sum += rank_dict[superkingdom][phylum][class_][order][family][genus]
			values.append(0)
  
	for superkingdom in rank_dict:
		for phylum in rank_dict[superkingdom]:
			for class_ in rank_dict[superkingdom][phylum]:
				names.append(class_)
				parents.append(phylum)
	  
				d_sum = 0
				for order in rank_dict[superkingdom][phylum][class_]:
					for family in rank_dict[superkingdom][phylum][class_][order]:
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							d_sum += rank_dict[superkingdom][phylum][class_][order][family][genus]
				values.append(0)
    
	for superkingdom in rank_dict:
		for phylum in rank_dict[superkingdom]:
			for class_ in rank_dict[superkingdom][phylum]:
				for order in rank_dict[superkingdom][phylum][class_]:
					names.append(order)
					parents.append(class_)
	  
					d_sum = 0
					for family in rank_dict[superkingdom][phylum][class_][order]:
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							d_sum += rank_dict[superkingdom][phylum][class_][order][family][genus]
					values.append(0)
     
	for superkingdom in rank_dict:
		for phylum in rank_dict[superkingdom]:
			for class_ in rank_dict[superkingdom][phylum]:
				for order in rank_dict[superkingdom][phylum][class_]:
					for family in rank_dict[superkingdom][phylum][class_][order]:
						names.append(family)
						parents.append(order)
	  
						d_sum = 0
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							d_sum += rank_dict[superkingdom][phylum][class_][order][family][genus]
						values.append(0)

	for superkingdom in rank_dict:
		for phylum in rank_dict[superkingdom]:
			for class_ in rank_dict[superkingdom][phylum]:
				for order in rank_dict[superkingdom][phylum][class_]:
					for family in rank_dict[superkingdom][phylum][class_][order]:
						for genus in rank_dict[superkingdom][phylum][class_][order][family]:
							names.append(genus)
							parents.append(family)
							values.append(rank_dict[superkingdom][phylum][class_][order][family][genus])
 
	fig = px.treemap(
		names=names,
		parents=parents,
		values=values,
	)
 
	fig.update_layout(margin = dict(t=50, l=25, r=25, b=25))
 
	fig.update_traces()
	fig.show()