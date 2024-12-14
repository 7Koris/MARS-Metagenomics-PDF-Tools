import argparse
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from plotnine import *
from ete3 import NCBITaxa
from sklearn import metrics

percent_regex = r"((\(+\d+%+\)))"
percent_value_regex = r"(\d+%)"
digit_regex = r"(\d+)"

def extract_percent(df: pd.DataFrame, col: str) -> None:
	new_col = col + "-confidence-%"
	col_idx = df.columns.get_loc(col)
	df.insert(col_idx + 1, new_col, df[col].str.extract(percent_value_regex))
	df[new_col] = df[new_col].str.replace("%", "").astype(float)
	df[col] = df[col].str.replace(percent_regex, "", regex=True)
	df[col] = df[col].str.strip()
	
def parse_file(file: str) -> tuple[pd.DataFrame, int]:
	of = open(file, "r")
	
	line_count = 0
	for line in of:
		line_count += 1
		
	df = pd.read_fwf(file, infer_nrows = line_count - 1)
	if "Unnamed: 1" in df.columns:
		df = df.drop(columns=["Unnamed: 1"]) 
	
	extract_percent(df, "superkingdom")
	extract_percent(df, "phylum")
	extract_percent(df, "genus")
	
	df["code"] = df["id"].str.extract(digit_regex)
	
	return (df, line_count - 1)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Bertax Confusion Matrix", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-f", "--files", nargs='+', help="Bertax Analyses Files")
	parser.add_argument("-s", "--seed", type=str, help="Seed for random number generator", default=0)

	args = parser.parse_args()
	config = vars(args)
	bertax_files = config["files"]
	seed = config["seed"]

	ncbi = NCBITaxa()
	sdf = pd.DataFrame()
	adf = []

	true_sk = 0
	true_ph = 0
	true_g = 0
	total = 0

	actual_sk = []
	predicted_sk = []

	actual_ph = []
	predicted_ph = []

	actual_g = []
	predicted_g = []


	for analysis_file in bertax_files:
		(df, line_count) = parse_file(analysis_file)
		adf = df
		for (idx, row) in df.iterrows():
			total += 1
			sk_name = row["superkingdom"]
			ph_name = row["phylum"]
			g_name = row["genus"]
			true_species_code = row["code"]
			
			if (sk_name == "Viruses"):
				continue
			
			lineage = ncbi.get_lineage(true_species_code)
			ranks = ncbi.get_rank(lineage)
			rank2id = {rank: taxid for (taxid, rank) in ranks.items()}
			sk_code = rank2id["superkingdom"]
			ph_code = rank2id["phylum"]
			g_code = rank2id["genus"]
			names = ncbi.translate_to_names([sk_code, ph_code, g_code])
			
			true_sk_name = names[0]
			true_ph_name = names[1]
			true_g_name = names[2]
			
			if (true_sk_name.lower() != "none" ):
				actual_sk.append(true_sk_name)
				predicted_sk.append(sk_name)
			if (true_ph_name.lower() != "none" ):
				actual_ph.append(true_ph_name)
				predicted_ph.append(ph_name)
			if (true_g_name.lower() != "none" ):
				actual_g.append(true_g_name)
				predicted_g.append(g_name)
			
	with PdfPages('bertax_cm.pdf') as pdf:
		labels = list(set(actual_sk + predicted_sk))
		confusion_matrix = metrics.confusion_matrix(actual_sk, predicted_sk, labels=labels, normalize="true")
		cm_display: metrics.ConfusionMatrixDisplay = metrics.ConfusionMatrixDisplay(confusion_matrix=confusion_matrix, display_labels=labels)
		cm_display.plot()
		pdf.savefig(plt.draw())
		
		confusion_matrix = metrics.confusion_matrix(actual_ph, predicted_ph)
		cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix=confusion_matrix)
		cm_display.plot(include_values=False)
		cm_display.ax_.set_xticks([])
		cm_display.ax_.set_yticks([])
		pdf.savefig(plt.draw())
		
		confusion_matrix = metrics.confusion_matrix(actual_g, predicted_g)
		cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix=confusion_matrix)
		cm_display.plot(include_values=False)
		cm_display.ax_.set_xticks([])
		cm_display.ax_.set_yticks([])
		pdf.savefig(plt.draw())
		