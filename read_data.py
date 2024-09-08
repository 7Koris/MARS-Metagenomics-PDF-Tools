import random
import sys
import re
import numpy as np
import pandas as pd
from taxonomy import TaxDict
from scipy import stats

class Read:
    """
    _summary_
    
    Class to represent a single read. Contains a hash to identify the read and an assigned taxon ID.
    """    
    
    assigned_taxon_id = None # Assigned taxon ID (Not necessarily species level)
    hash = None # Unique hash to identify read
    
    def get_assignment(self):
        return self.assigned_taxon_id
    
    def __init__(self, hash: str, id):
        self.hash = hash
        self.assigned_taxon_id = id
    
     
class ReadData:
    """
    _summary_

    ReadData is a class for loading and processing read data 
    from MetaMaps and MTSV output files. Supports filtering by 
    frequency and coverage. The usage of MTSV data requires a 
    lookup file to be present. The lookup file is used to match 
    read IDs to hashes between MTSV and MetaMaps output files.
    """    
    
    _id_blacklist_dict: dict = {} # Dict of ids to exclude/ignore of format: <numeric taxon id, Bool>, True if excluded. TODO: Implement
    _tax_dict: TaxDict # Dictionary of taxon IDs and read counts assigned for each ID
    _incidence_dict: dict = {} # Dictionary of taxon IDs based on being seen in reads <taxonID, bool>
    
    reads = dict() # Dictionary of reads given hash
    processed_read_count = 0 # Total number of reads processed including those discarded/pruned
    
    # TODO: Implement BLACKLIST
    def __init__(self, rebuild_database = False, id_blacklist_dict: dict = {}) -> None:
        self._id_blacklist_dict = id_blacklist_dict
        self._tax_dict = TaxDict(rebuild_database)
        self._incidence_dict = dict()
        self.reads = dict()

    
    def load_data(self, data: tuple) -> None:
        """
        _summary_

        Loads a tuple of data into the class. Used for loading data from a pickle file.

        Args:
            data (tuple): _description_: Tuple of data to load in the format (reads, processed_read_count)
        """        
        reads = data[0]
        processed_read_count = data[1]
        self.reads = reads
        self.processed_read_count = processed_read_count
        
    
    def get_tax_count_dict(self) -> TaxDict:
        """
        _summary_
        
        Get a dictionary of taxon IDs and read counts assigned for each ID for all available NCBI levels

        Returns:
            TaxDict: _description_: Dictionary of taxon IDs and read counts assigned for each ID
        """        
        
        self._tax_dict = TaxDict(False)
        for read in self.reads.values():
            if isinstance(read.assigned_taxon_id, list):
                print("Warning: Read %s assigned to multiple taxon IDs, cannot insert into tax_dict. Please resolve to LCA first using resolve_lca()" % read.hash)
                continue
            self._tax_dict.insert_id(read.assigned_taxon_id, 1)
        return self._tax_dict
     
        
    def insert_read(self, new_read: Read) -> None:
        """
        _summary_
        
        Insert an instantiated read object into the read database.

        Args:
            new_read (Read): _description_: Read object to insert
        """        
        self.processed_read_count += 1
        id = new_read.assigned_taxon_id
        if isinstance(id, list):
            pass
        elif id in self._id_blacklist_dict: # Skip blacklisted ids
            #print("Skipping blacklisted ID %s" % id)
            return 
        
        if new_read.hash in self.reads:
            print("Error: Read %s already parsed" % new_read.hash)
            sys.exit(1)
        self.reads[new_read.hash] = new_read
    
    
    def insert_hash_as_read(self, hash: str, assigned_id: int) -> None: 
        """
        _summary_
        
        Insert a hash as a read into the read database.

        Args:
            hash (str): _description_: Unique hash to identify read (as determined by MetaMaps/MTSV)
            assigned_id (int): _description_: Taxon ID assigned to read
        """        
        
        new_read = Read(hash, assigned_id)
        self.insert_read(new_read)
    
    
    def prune_non_incidental_reads(self) -> int:
        """
        _summary_
        
        Prune reads that are not incidental (not seen in incidence_dict). 
        This means that if the reads are MTSV, and the incidence dict was 
        built with MetaMaps, only reads that are seen in both will be kept.

        Returns:
            int: _description_: Number of reads pruned
        """        
        
        if len(self.reads) == 0:
            return

        read_keys = list(self.reads.keys())
        reads_deleted = 0
        start_read_count = len(read_keys)
        for read in read_keys:
            if read not in self._incidence_dict:
                self.reads.pop(read)
                reads_deleted += 1
        print("Pruned %d reads" % reads_deleted + " of %d" % start_read_count, "by incidence")
        return reads_deleted
    
    
    def prune_by_level(self, level: str) -> int:
        """
        _summary_
        
        Prune all reads that are missing classification at the level

        Args:
            level (str): _description_: NCBI rank level

        Returns:
            int: _description_: Number of reads pruned
        """        
        
        if len(self.reads) == 0:
            return
        
        read_keys = list(self.reads.keys())
        reads_deleted = 0
        start_read_count = len(read_keys)
        for read in read_keys:
            levels = self._tax_dict.ncbi.get_lineage(self.reads[read].assigned_taxon_id)
            level_dict = {}
            for idx, id in enumerate(levels):
                levels[idx] = self._tax_dict.id_2_rank(id)
                level_dict[levels[idx]] = 1
                
            if level not in level_dict:
                self.reads.pop(read)
                reads_deleted += 1
        print("Pruned %d reads" % reads_deleted + " of %d" % start_read_count, "at level", level)
        return reads_deleted
    
    
    def prune_reads_not_in_rank(self, rank_level: str, rank_name: str) -> int:
        """
        _summary_
        
        Remove all reads that do not have an NCBI rank of rank_name at rank_level

        Args:
            rank_level (str): _description_: NCBI rank level
            rank_name (str): _description_: NCBI rank name

        Returns:
            int: _description_: Number of reads pruned
        """        
        
        if len(self.reads) == 0:
            return
        
        read_keys = list(self.reads.keys())
        reads_deleted = 0
        start_read_count = len(read_keys)
        for read in read_keys:
            assigned_id = self.reads[read].assigned_taxon_id
            lineage = self._tax_dict.ncbi.get_lineage(assigned_id)
            ranks = self._tax_dict.ncbi.get_rank(lineage)
            rank_names = []
            for id in lineage:
                name = self._tax_dict.id_2_name(id)
                rank_names.append(name)
                
            if rank_name in rank_names:
                pass
            else:
                self.reads.pop(read)
                reads_deleted += 1
                continue
            
            for idx, id in enumerate(lineage):
                current_rank_name = self._tax_dict.id_2_rank(id)
                
                if current_rank_name == rank_level:
                    current_name = self._tax_dict.id_2_name(id)
                    if current_name == rank_name:
                        break
                    else:
                        self.reads.pop(read)
                        reads_deleted += 1
                        break
        print("Pruned %d reads" % reads_deleted + " of %d" % start_read_count, "not in rank", rank_name)
    
    
    
    def resolve_lca(self) -> None:
        """
        _summary_
        
        When there are multiple taxon IDs assigned to a read, resolve to the LCA of the taxon IDs.
        IT IS REQUIRED TO RUN THIS FUNCTION BEFORE USING THE READ DATA FOR ANALYSIS.
        
        """        
        
        for read in self.reads:
            current_read = self.reads[read]
            if isinstance(current_read.assigned_taxon_id, list):
                lca = self._tax_dict.get_lca_id(current_read.assigned_taxon_id)
                if lca == 0:
                    print("Error: No LCA found for read %s" % read)
                    sys.exit(1)
                current_read.assigned_taxon_id = lca
                
                
    def rarefy(self, min_reads: int, seed: int = 0) -> None:
        print("Rarefying reads down to %d" % min_reads, "with seed %d" % seed)
        keys = list(self.reads.keys())
        
        random.seed(seed)
        random.shuffle(keys)
        
        if len(keys) <= min_reads:
            return
        
        while len(keys) > min_reads:
            random_key = random.choice(keys)
            self.reads.pop(random_key)
            keys.remove(random_key)

        
    def parse_mtsv_reads(self, read_file_name: str, lookup_file_name: str, incidence_only: bool = False) -> None:
        lookup_dict = {} 
        file = lookup_file_name
        opened_file = open(file, "r")
        lines = opened_file.readlines()
        opened_file.close()
        
        for line in lines:
            (read_id, hash) = line.split()
            lookup_dict[read_id] = hash
        
        file = read_file_name
        opened_file = open(file, "r")
        lines = opened_file.readlines()
        opened_file.close()
        
        for line in lines:
            read_id = line.split(":")[0]
            hash = lookup_dict[read_id]
            assignments = re.findall(":(.*=.*|,)", line)
            assignments = re.split(",", assignments[0])            
            min_distance = None
            min_ids = []
            
            for assignment in assignments:
                assignment = re.split("=", assignment)
                
                if len(assignment) != 2:
                    sys.exit("Error: Invalid assignment format")
                    
                id = assignment[0]
                distance = assignment[1]
                
                if min_distance == None or distance < min_distance:
                    min_distance = distance
                    
            for assignment in assignments:
                assignment = re.split("=", assignment)
                id = assignment[0]
                distance = assignment[1]
                
                if distance == min_distance:
                    min_ids.append(id)
            
            if incidence_only:
                self._incidence_dict[hash] = 1
                continue

            if len(min_ids) == 1:
              new_read = Read(hash, min_ids[0])
              self.insert_read(new_read)
            elif len(min_ids) > 1:
                new_read = Read(hash, min_ids)
                self.insert_read(new_read)
            else:
                sys.exit("Error: No taxon ID found for read %s" % read_id, "Is there a formatting error?")


    def parse_metamaps_reads_2_taxon(self, reads_2_taxon_file_name: str, incidence_only: bool = False) -> None: 
        file = reads_2_taxon_file_name
        opened_file = open(file, "r")
        lines = opened_file.readlines()
        opened_file.close()
        for line in lines:
            (hash, id) = line.split("\t")
            
            id = id.strip()
            hash = hash.strip()
            
            if id == "0": # Skip unassigned reads
                self.processed_read_count += 1
                continue
            
            if incidence_only:
                self._incidence_dict[hash] = 1
                continue
            
            new_read = Read(hash, id)
            
            self.insert_read(new_read)
            


  
class ReadDataMM:  
    """
    ReadDataMM is a class for loading and processing read data from MetaMaps output files. Supports filtering by frequency and coverage.
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
    