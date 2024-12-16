import sys
import re
import numpy as np
import pandas as pd
from taxonomy import TaxDict
from scipy import stats

class MetaMapsIdentityData:  
    """
    MetaMapsIdentityData is a class for loading and processing read data from MetaMaps output files. Supports filtering by frequency and coverage.
        Inelegant and slow, but functional.
    """
    
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
    tax_id_is_outlier: dict
    excluded_id_dict: dict
    min_freq: float = 0.0
    
    tax_dict: TaxDict
    reads: list
    
    def reset(self) -> None:
        self.file_prefix = None
        self.mapping_units = None
        self.coverage_data = None
        self.counts_per_unit = None
        self.freq_per_unit = None
        self.tax_id_2_mapping_units = None
        self.mapping_unit_2_tax_id = None
        self.filtered_tax_ids = None
        self.tax_id_2_filtered_contigs = None
        self.tax_id_2_all_contigs = None
        self.tax_id_2_name = None
        self.contig_2_coverages = None
        self.are_coverages_loaded = None
        self.tax_id_is_outlier = None
        self.excluded_id_dict = None
        self.min_freq = None
        self.tax_dict = TaxDict()
        self.reads = []
  
  
    def parse_metamaps_data(self, read_file_prefix: str, min_freq: float, excluded_taxon_ids: list[int] = []):
        lengths_file = read_file_prefix + ".EM.lengthAndIdentitiesPerMappingUnit"
        lengths_and_ids = pd.read_csv(lengths_file, delimiter="\t")
        lengths_and_ids = lengths_and_ids[lengths_and_ids["AnalysisLevel"] == "EqualCoverageUnit"]
        self.mapping_units = lengths_and_ids
        self.counts_per_unit = pd.Series(lengths_and_ids["ID"].value_counts().sort_values(ascending=False)) # Count of each taxonomic ID's occurrences
        self.freq_per_unit = self.counts_per_unit / self.counts_per_unit.sum() # Frequency of each taxonomic ID's
        
        self.excluded_id_dict = {}
        excluded_tax_count = 0
        for id in excluded_taxon_ids:
            self.excluded_id_dict[str(id)] = 1
            excluded_tax_count += 1
        
        self.tax_id_2_mapping_units = {}
        self.mapping_unit_2_tax_id = {}
        self.filtered_tax_ids = {}
        self.tax_id_2_filtered_contigs = {}
        for i in range(0, len(self.counts_per_unit)):
            current_id_label = self.counts_per_unit.index[i]
            matches = re.findall("kraken:taxid\\|(x?\\d+)\\|", current_id_label)
            
            if len(matches) != 1:
                print("Error: No recognizable taxonomic ID in " + current_id_label + ", pelase check formatting")
                sys.exit(1)
        
            taxon_id = matches[0]
            
            if taxon_id not in self.excluded_id_dict:
                if (taxon_id not in self.tax_id_2_mapping_units.keys()):
                    self.tax_id_2_mapping_units[taxon_id] = {}
                self.tax_id_2_mapping_units[taxon_id][current_id_label] = 1
                self.mapping_unit_2_tax_id[current_id_label] = taxon_id
                self.filtered_tax_ids[taxon_id] = 1
                
        self._filter_by_frequency(min_freq)
   
    
    def __init__(self, read_file_prefix, min_freq, excluded_taxon_ids: list[int] = []) -> None:
        self.reset()
        print("Initializing MappingUnitData")
        self.file_prefix = read_file_prefix
        self.min_freq = min_freq
        self.parse_metamaps_data(read_file_prefix, min_freq, excluded_taxon_ids)
        self._init_tax_dict()
        
        
    def _filter_by_frequency(self, min_freq: float) -> dict:
        temp_tax_id_2_mapping_units = {}
        temp_mapping_unit_2_tax_id = {}
        temp_filtered_tax_ids = {}
        for current_tax_id in self.filtered_tax_ids.keys():
            if current_tax_id in self.excluded_id_dict:
                continue
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
    
    
    def load_coverage(self) -> None:        
        print("Loading coverage data")    
        coverage_file = self.file_prefix + ".EM.contigCoverage"
        self.coverage_data = pd.read_csv(coverage_file, delimiter="\t")
        
        self.tax_id_2_filtered_contigs = {}
        self.tax_id_2_all_contigs = {}
        self.contig_2_coverages = {}
        self.tax_id_2_name = {}

        for idx, id, name, contig, start, stop, n_bases, coverage in self.coverage_data.itertuples():
            if id in self.excluded_id_dict:
                continue
            
            id = str(id)
            contig = str(contig)
            coverage = float(coverage)
            if id not in self.tax_id_2_name:
                self.tax_id_2_name[id] = name
            if id in self.filtered_tax_ids:
                # update mapping_unit_2_tax_id
                if contig in self.mapping_unit_2_tax_id.keys():
                    self.mapping_unit_2_tax_id[contig] = id
                
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
    
    
    def filter_sig_bin_outliers(self, non_zero_count_threshold: int=3, preserve_outliers: bool=False) -> int:
        print("Filtering by significant bins")
        if self.are_coverages_loaded == False:
            print("Error: No coverage data loaded, please load coverage data before filtering by coverage")
            return
        
        o_count = 0
        temp_filtered_tax_ids = dict(self.filtered_tax_ids)
        self.tax_id_is_outlier = {}
        for current_tax_id in temp_filtered_tax_ids.keys():
            if current_tax_id in self.excluded_id_dict:
                continue
            
            if current_tax_id not in self.tax_id_is_outlier:
                self.tax_id_is_outlier[current_tax_id] = False

            all_coverages = []
            current_contigs = self.tax_id_2_all_contigs[current_tax_id]
            for contig in current_contigs:
                current_coverages = self.contig_2_coverages[contig]
                all_coverages.extend(current_coverages)
            

            indices = []
            values = []
            i = 0
 
            for val in all_coverages:
                if val != 0:
                    indices.append(i)
                    values.append(val)
                i += 1
       
            method = 'fd'
            density = True
            hist, edges = np.histogram(indices, bins=600, density=density, weights=values)       
            nonzero_count = 0
            hist_devs = stats.zscore(hist)
            for i in range(0, len(hist)):
                if hist_devs[i] > .025:
                    nonzero_count += 1
                    
            if (nonzero_count <= non_zero_count_threshold):
                self.tax_id_is_outlier[current_tax_id] = True
                o_count += 1

        outliers = []
        for id in self.filtered_tax_ids.keys():
            if self.tax_id_is_outlier[id]:
                outliers.append(id)

        if not preserve_outliers:
            print("Rebuilding database with", o_count, "outliers removed")
            self.__init__(self.file_prefix, self.min_freq, outliers)
            self.load_coverage() # todo: rebuild hashmaps faster than re-reading file(s)
        
        return o_count
        
    
    def _init_tax_dict(self) -> None:
        self.tax_dict = TaxDict()
        for current_tax_id in self.filtered_tax_ids.keys():
            if current_tax_id in self.excluded_id_dict:
                continue
            current_units = self.tax_id_2_mapping_units[current_tax_id]
            for unit in current_units:
                self.tax_dict.insert_id(current_tax_id, self.counts_per_unit[unit])
                for i in range(0, self.counts_per_unit[unit]):
                    self.reads.append((current_tax_id))
        

    def get_abundance_estimates(self) -> dict:
        read_count_dict = {}
        for current_tax_id in self.filtered_tax_ids.keys():
            if current_tax_id in self.excluded_id_dict:
                continue
            current_units = self.tax_id_2_mapping_units[current_tax_id]
            for unit in current_units:
                if current_tax_id not in read_count_dict:
                    read_count_dict[current_tax_id] = self.counts_per_unit[unit]
                else:
                    read_count_dict[current_tax_id] += self.counts_per_unit[unit]
        sorted_read_count_dict = {k: v for k, v in sorted(read_count_dict.items(), key=lambda item: item[1])}
        return sorted_read_count_dict          