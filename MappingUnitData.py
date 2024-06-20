import sys
import re
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import argparse

class MappingUnitData:
    mapping_units: pd.DataFrame
    counts_per_unit: pd.Series
    freq_per_unit: pd.Series
    tax_id_2_mapping_units: dict
    mapping_unit_2_tax_id: dict
    filtered_tax_ids: dict 
    
    
    def __init__(self, read_file_prefix: str = None) -> None:
        if read_file_prefix is None:
            return
        lengths_file = read_file_prefix + ".EM.lengthAndIdentitiesPerMappingUnit"
        lengths_and_ids = pd.read_csv(lengths_file, delimiter="\t")
        lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"]
        self.mapping_units = lengths_and_ids
        self.counts_per_unit = pd.Series(lengths_and_ids["ID"].value_counts().sort_values(ascending=False)) # Count of each taxonomic ID's occurrences
        self.freq_per_unit = self.counts_per_unit / self.counts_per_unit.sum() # Frequency of each taxonomic ID's
        
        self.tax_id_2_mapping_units = {}
        self.mapping_unit_2_tax_id = {}
        self.filtered_tax_ids = {}
        for i in range(0, len(self.counts_per_unit)):
            current_id_label = self.counts_per_unit.index[i]
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label)
            
            if len(matches) != 1:
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
        
            taxon_id = matches[0]
            
            if (taxon_id not in self.tax_id_2_mapping_units.keys()):
                self.tax_id_2_mapping_units[taxon_id] = {}
            self.tax_id_2_mapping_units[taxon_id][current_id_label] = 1
            self.mapping_unit_2_tax_id[current_id_label] = taxon_id
            self.filtered_tax_ids[taxon_id] = 1
            
            
    def load_reads(self, file_prefix: str) -> None:
        self.__init__(self, file_prefix)


    def filter_by_frequency(self, min_freq: float) -> dict:
        temp_tax_id_2_mapping_units = {}
        temp_mapping_unit_2_tax_id = {}
        temp_filtered_tax_ids = {}
        for current_tax_id in self.filtered_tax_ids.keys():
            current_units = self.tax_id_2_mapping_units[current_tax_id]
            for unit in current_units:
                current_freq = self.freq_per_unit[unit]
                if (current_freq >= min_freq):
                    if (current_tax_id not in temp_tax_id_2_mapping_units.keys()):
                        temp_tax_id_2_mapping_units[current_tax_id] = {}
                    temp_tax_id_2_mapping_units[current_tax_id][unit] = 1
                    temp_mapping_unit_2_tax_id[unit] = current_tax_id
                    temp_filtered_tax_ids[current_tax_id] = 1
        self.tax_id_2_mapping_units = temp_tax_id_2_mapping_units
        self.mapping_unit_2_tax_id = temp_mapping_unit_2_tax_id
        self.filtered_tax_ids = temp_filtered_tax_ids
        

    def get_abundance_estimates(self) -> dict:
        read_count_dict = {}
        for current_tax_id in self.filtered_tax_ids.keys():
            current_units = self.tax_id_2_mapping_units[current_tax_id]
            for unit in current_units:
                if current_tax_id not in read_count_dict:
                    read_count_dict[current_tax_id] = self.counts_per_unit[unit]
                else:
                    read_count_dict[current_tax_id] += self.counts_per_unit[unit]
        sorted_read_count_dict = {k: v for k, v in sorted(read_count_dict.items(), key=lambda item: item[1])}
        return sorted_read_count_dict