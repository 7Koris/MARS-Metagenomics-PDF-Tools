import math
import composition_stats as cs
import numpy as np

"""
Collection of beta diversity metric functions for microbial communities.
"""

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


def aitchison_distance(mu1_abundances, mu2_abundances, all_ids: dict[str, int]):
    clr_mu1 = cs.clr(np.array(list(mu1_abundances.values())))
    clr_mu1_dict = {}
    for i, id in enumerate(mu1_abundances.keys()):
        clr_mu1_dict[id] = clr_mu1[i]
  
    clr_mu2 = cs.clr(np.array(list(mu2_abundances.values())))
    clr_mu2_dict = {}
    for i, id in enumerate(mu2_abundances.keys()):
        clr_mu2_dict[id] = clr_mu2[i]
    
    sum = 0
    for id in all_ids.keys():
        if id in clr_mu1_dict and id in clr_mu2_dict:
            sum += (clr_mu1_dict[id] - clr_mu2_dict[id]) ** 2
        elif id in clr_mu1_dict:
            sum += clr_mu1_dict[id] ** 2
        elif id in clr_mu2_dict:
            sum -= clr_mu2_dict[id] ** 2
    
    return math.sqrt(sum)
    
    