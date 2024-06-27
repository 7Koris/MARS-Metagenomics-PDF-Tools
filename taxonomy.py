from ete3 import NCBITaxa

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
            
    return tax_dict
        


    
    