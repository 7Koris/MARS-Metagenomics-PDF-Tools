
import argparse
import pandas as pd
import scripts.srs as srs
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from plotnine import *

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
	
    return (df, line_count - 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bertax Summary", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file", type=str, help="File name", default=None, required=True)
    parser.add_argument("-sk", "--skref", type=str, help="If data is from uniform superkingdom, compare", default="")
    
    args = parser.parse_args()
    config = vars(args)
    
    (df, count) = parse_file(config["file"])
    S = df["superkingdom"].value_counts()
    P = df["phylum"].value_counts()
    G = df["genus"].value_counts()
    
    S = S.reset_index()
    P = P.reset_index()
    G = G.reset_index()

    sdf = pd.DataFrame(S)
    pdf = pd.DataFrame(P)
    gdf = pd.DataFrame(G)
    
    for row in sdf.itertuples():
        sdf.loc[row.Index, "superkingdom"] = row.superkingdom + " (" + str(row.count) + ")"
    
    match config["skref"].lower():
        case "bacteria":
            acc = sdf.loc[sdf["superkingdom"] == "Bacteria"]["count"].values[0]
            acc = acc / float(count)
            print(f"Acc: {acc}")
        case "archaea":
            acc = sdf.loc[sdf["superkingdom"] == "Archaea"]["count"].values[0]
            acc = acc / float(count)
            print(f"Acc: {acc}")
        case "eukaryota":
            acc = sdf.loc[sdf["superkingdom"] == "Eukaryota"]["count"].values[0]
            acc = acc / float(count)
            print(f"Acc: {acc}")
        case _:
            pass
    
    
    
    gg1 = ggplot(sdf, aes(x="superkingdom", y="count")) + geom_bar(stat="identity") #+ theme(axis_text_x=element_text(rotation=90, hjust=1))
    
    with PdfPages('multiple_plots.pdf') as pdf:
        pdf.savefig(gg1.draw())
  
    
   