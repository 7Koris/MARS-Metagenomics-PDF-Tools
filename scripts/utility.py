import math
import math
import numpy as np
import mm_identity_data as mmid

def get_coverage_outllier_list(file_prefix) -> list:
    mapping_units = mmid.MetaMapsIdentityData(file_prefix, 0.0)
   
    mapping_units.load_coverage()
    mapping_units.filter_sig_bin_outliers(3, True)
    
    outliers = []
    for id in mapping_units.filtered_tax_ids.keys():
        if mapping_units.tax_id_is_outlier[id]:
            outliers.append(id)
            
    return outliers

"""
Collection of alpha diversity metric functions for microbial communities.
"""
class Alphas:
    def chao1(read_counts: list[float]):
    # https://www.cd-genomics.com/microbioseq/the-use-and-types-of-alpha-diversity-metrics-in-microbial-ngs.html#:~:text=1.,The%20index%20of%20Community%20richness&text=The%20Chao%20estimator%20is%20employed,species%20diversity%20within%20the%20sample.

        species_count = len(read_counts)
        singles = 0
        doubles = 0
        for count in read_counts:
            if count == 1:
                singles += 1
            elif count == 2:
                doubles += 1
            
        return species_count + (singles * (singles - 1)) / (2 * (doubles + 1))


    def shannons_alpha(proportions: list[float]):
    # Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
        n = sum(proportions)
        s = len(proportions)
        
        p = []
        d = 0.0
        
        for i in proportions:
            p.append(i / n)
            d += i * (i - 1)
        d = d / (n * (n - 1.0))
        
        h = []
        for i in p:
            h.append(i * math.log(i))

        return (-1 *sum(h))


    def berger_parkers_alpha(proportions: list[float]):
    # Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
        n = sum(proportions)
        s = len(proportions)
        
        p = []
        d = 0.0
        
        for i in proportions:
            p.append(i / n)
            d += i * (i - 1)
        d = d / (n * (n - 1.0))

        return max(p)


    def simpsons_alpha(proportions: list[float]):
    # Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
        n = sum(proportions)
        s = len(proportions)
        
        p = []
        d = 0.0
        
        for i in proportions:
            p.append(i / n)
            d += i * (i - 1)
        d = d / (n * (n - 1.0))
        
        return 1-d


    def inverse_simpsons_alpha(proportions: list[float]):
    # Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
        n = sum(proportions)
        s = len(proportions)
        
        p = []
        d = 0.0
        
        for i in proportions:
            p.append(i / n)
            d += i * (i - 1)
        d = d / (n * (n - 1.0))
        
        return 1/d

"""
Collection of beta diversity metric functions for microbial communities.
"""
class Betas:
    def bray_curtis_distance(mu1_abundances, mu2_abundances, all_ids: dict[str, int]):
        cij = 0
        sij = 0
        
        for id in all_ids.keys():
            if id in mu1_abundances and id in mu2_abundances:
                cij += min(mu1_abundances[id], mu2_abundances[id])
                sij += mu1_abundances[id] + mu2_abundances[id]
            elif id in mu1_abundances:
                sij += mu1_abundances[id]
            elif id in mu2_abundances:
                sij += mu2_abundances[id]
        if cij == 0:
            return 1 # no overlap
        if sij == 0:
            return -1 # no reads
        
        return 1 - ( (2 * cij) / sij)


    def jaccard_distance(mu1_abundances, mu2_abundances, all_ids: dict[str, int]):
        num_common = 0
        for id in all_ids.keys():
            if id in mu1_abundances and id in mu2_abundances:
                if mu1_abundances[id] > 0 and mu2_abundances[id] > 0:
                    num_common += 1
            else:
                pass

        return 1 - (num_common / len(all_ids))


    def city_block_distance(mu1_abundances, mu2_abundances, all_ids: dict[str, int]):
        distance = 0
        for id in all_ids.keys():
            if id in mu1_abundances and id in mu2_abundances:
                distance += abs(mu1_abundances[id] - mu2_abundances[id])
            elif id in mu1_abundances:
                distance += mu1_abundances[id]
            elif id in mu2_abundances:
                distance += mu2_abundances[id]
        return distance
        
        