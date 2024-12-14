import math

"""
Collection of alpha diversity metric functions for microbial communities.
"""

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