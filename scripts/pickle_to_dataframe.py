import pandas as pd
import pickle as pk
import read_data as rd

def pickle_to_dataframe(pickle_file_path) -> pd.DataFrame:
	data = pk.load(open(pickle_file_path, "rb"))
	read_data = rd.ReadData()
	read_data.load_data(data)
	reads = [read.get_assignment() for read in read_data.reads.values()]
	hashes = [read.get_hash() for read in read_data.reads.values()]
	df = pd.DataFrame({'Hashes': hashes, 'Taxon ID': reads})
	print(df)
	
	return df

pickle_to_dataframe("C:\\Users\\Koris\\Documents\\MARS-Metagenomics-PDF-Tools\\reads\\svcap241_1k_sigbin.p")
