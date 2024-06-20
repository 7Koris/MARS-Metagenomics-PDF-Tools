import sys
import re
from scipy import stats
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


# todo: Add read count integrity steps
class MappingUnitData:
    file_prefix: str
    mapping_units: pd.DataFrame
    coverage_data: pd.DataFrame
    counts_per_unit: pd.Series
    freq_per_unit: pd.Series
    tax_id_2_mapping_units: dict
    mapping_unit_2_tax_id: dict
    filtered_tax_ids: dict
    tax_id_2_filtered_contigs: dict
    tax_id_2_all_contigs: dict
    tax_id_2_name: dict
    contig_2_coverages: dict
    are_coverages_loaded: bool = False
    
    
    def __init__(self, read_file_prefix) -> None:
        self.file_prefix = read_file_prefix
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
        
    
    def load_coverage(self) -> None:            
        coverage_file = self.file_prefix + ".EM.contigCoverage"
        self.coverage_data = pd.read_csv(coverage_file, delimiter="\t")
        
        self.tax_id_2_filtered_contigs = {}
        self.tax_id_2_all_contigs = {}
        self.contig_2_coverages = {}
        self.tax_id_2_name = {}

        for idx, id, name, contig, start, stop, n_bases, coverage in self.coverage_data.itertuples():
            id = str(id)
            contig = str(contig)
            coverage = float(coverage)
            if id not in self.tax_id_2_name:
                self.tax_id_2_name[id] = name
            if id in self.filtered_tax_ids:
                # update tax_id_2_all_contigs
                if (id not in self.tax_id_2_all_contigs):
                    self.tax_id_2_all_contigs[id] = {}
                self.tax_id_2_all_contigs[id][contig] = 1
                
                # update tax_id_2_filtered_contigs
                if contig in self.tax_id_2_mapping_units[id].keys():
                    if (id not in self.tax_id_2_filtered_contigs.keys()):
                        self.tax_id_2_filtered_contigs[id] = {}
                    self.tax_id_2_filtered_contigs[id][contig] = 1

                # update contig_2_coverages
                if (contig not in self.contig_2_coverages.keys()):
                    self.contig_2_coverages[contig] = []
                self.contig_2_coverages[contig].append(coverage)
        self.are_coverages_loaded = True


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
        
        if self.are_coverages_loaded:
            self.load_coverage() # todo: rebuild hashmaps faster than re-reading file
    
    
    def filter_by_max_avg_coverage(self, max_outlier_coverage: float, proportion) -> dict:
        if self.are_coverages_loaded == False:
            print("Error: No coverage data loaded, please load coverage data before filtering by coverage")
            return
        
        o_count = 0
        
        temp_tax_id_2_mapping_units = {}
        temp_mapping_unit_2_tax_id = {}
        temp_filtered_tax_ids = {}
        for current_tax_id in self.filtered_tax_ids.keys():
            current_contigs = self.tax_id_2_all_contigs[current_tax_id]
            all_coverages = []
            for contig in current_contigs:
                current_coverages = self.contig_2_coverages[contig]
                all_coverages.extend(current_coverages)
            tm = stats.trim_mean(all_coverages, proportion)
            if (not tm <= max_outlier_coverage):
                if (current_tax_id not in temp_tax_id_2_mapping_units.keys()):
                    temp_tax_id_2_mapping_units[current_tax_id] = {}
                for contig in current_contigs:
                    temp_tax_id_2_mapping_units[current_tax_id][contig] = 1
                    temp_mapping_unit_2_tax_id[contig] = current_tax_id
                temp_filtered_tax_ids[current_tax_id] = 1
            else:
                o_count += 1
        self.tax_id_2_mapping_units = temp_tax_id_2_mapping_units
        self.mapping_unit_2_tax_id = temp_mapping_unit_2_tax_id
        self.filtered_tax_ids = temp_filtered_tax_ids
        print("Number of outliers: %s" %o_count)
        
        self.load_coverage() # todo: rebuild hashmaps faster than re-reading file

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
    
    
    def construct_plot_dict(self) -> dict:
        if self.are_coverages_loaded == False:
            print("Error: No coverage data loaded, please load coverage data before generating plot dict")
            
        plot_dict = {}
       
        return plot_dict
                
                