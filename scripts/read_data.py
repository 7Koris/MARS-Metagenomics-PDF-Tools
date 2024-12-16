import random
import sys
import re
from taxonomy import TaxDict

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
    
    _tax_dict: TaxDict # Dictionary of taxon IDs and read counts assigned for each ID
    _incidence_dict: dict = {} # Dictionary of taxon IDs based on being seen in reads <taxonID, bool>
    
    reads = dict() # Dictionary of reads given hash
    processed_read_count = 0 # Total number of reads processed including those discarded/pruned
    
    def __init__(self, rebuild_database = False) -> None:
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
            


  

    