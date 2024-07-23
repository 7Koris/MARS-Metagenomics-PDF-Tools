from ete3 import NCBITaxa
from ete3 import Tree

class TaxDict(dict):
    ncbi = None
    hashes = {}
    
    def __init__(self, update_local: bool=False):
        self.ncbi = NCBITaxa()
        if update_local:
            self.ncbi.update_taxonomy_database()
        self.hashes = {}
            
        
    def insert_id(self, id: str, count: int, hash: str=None):
        lineage = self.ncbi.get_lineage(id)
        ranks = list(self.ncbi.get_rank(lineage).values())
        names = list(self.ncbi.get_taxid_translator(lineage).values())
        for idx, rank in enumerate(ranks):
            if rank == "no rank":
                continue
            if rank not in self:
                self[rank] = {}
            name = names[idx]
            if name not in self[rank]:
                self[rank][name] = 0
            self[rank][name] += count
            
            if rank not in self.hashes:
                self.hashes[rank] = {}
            self.hashes[rank][hash] = True
    
            
    def id_2_name(self, id: str) -> str:
        return self.ncbi.get_taxid_translator([id])[id]
    
    
    def id_2_rank(self, id: str) -> str:
        return self.ncbi.get_rank([id])[id]
    
    
    def id_2_level_id(self, id: str, level: str) -> str:
        levels = self.ncbi.get_lineage(id)
        ranks = self.ncbi.get_rank(levels)
        for idx, id in enumerate(levels):
            if ranks[id] == level:
                return id
        return None
        #exit("Level not found")
            
            
    def insert_lca(self, count, ids: list):
        lca_id = self.get_lca_id(ids)
        self.insert_id(lca_id, count)
        
    
    def get_lca_id(self, ids: list):
        lineages = []
        for id in ids:
            lineages.append(self.ncbi.get_lineage(id))
     
        last_tax_id = 1
        for tax_id in lineages[0]:
            id_in_siblings = True
            for lineage in lineages:
                if lineage != lineages[0]:
                    if tax_id not in lineage:
                        id_in_siblings = False
                        break
                    
            if not id_in_siblings:
                return last_tax_id
            else:
                last_tax_id = tax_id
        return last_tax_id
    

def get_lineage_counts(abundance_estimates: dict, update_local: bool=True) -> dict:
    # Abundance estimates has format <tax_id: str, count: int>
    tax_dict = {}
    
    ncbi = NCBITaxa()
    if update_local:
        ncbi.update_taxonomy_database()
        
    for id in abundance_estimates.keys():
        lineage = ncbi.get_lineage(id)
        ranks = list(ncbi.get_rank(lineage).values())
        names = list(ncbi.get_taxid_translator(lineage).values())
        for idx, rank in enumerate(ranks):
            if rank == "no rank":
                continue
            if rank not in tax_dict:
                tax_dict[rank] = {}
            name = names[idx]
            if name not in tax_dict[rank]:
                tax_dict[rank][name] = 0
            tax_dict[rank][name] += abundance_estimates[id]
    print(tax_dict.keys())
            
    return tax_dict
        
    