import math


# https://www.cd-genomics.com/microbioseq/the-use-and-types-of-alpha-diversity-metrics-in-microbial-ngs.html#:~:text=1.,The%20index%20of%20Community%20richness&text=The%20Chao%20estimator%20is%20employed,species%20diversity%20within%20the%20sample.
def chao1(read_counts: list[float]):
    species_count = len(read_counts)
    singles = 0
    doubles = 0
    for count in read_counts:
        if count == 1:
            singles += 1
        elif count == 2:
            doubles += 1
        
    return species_count + (singles * (singles - 1)) / (2 * (doubles + 1))

# def ace(read_counts: list[float]):
#     species_count = len(read_counts)
#     singles = 0
#     doubles = 0
#     for count in read_counts:
#         if count == 1:
#             singles += 1
#         elif count == 2:
#             doubles += 1

#     return species_count + (singles * (singles - 1)) / (2 * (doubles + 1))


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def shannons_alpha(proportions: list[float]):
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
    #print("Shannon's diversity: %s" %(-1 *sum(h)))
    return (-1 *sum(h))


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def berger_parkers_alpha(proportions: list[float]):
    n = sum(proportions)
    s = len(proportions)
    
    p = []
    d = 0.0
    
    for i in proportions:
        p.append(i / n)
        d += i * (i - 1)
    d = d / (n * (n - 1.0))
    # print("Berger-parker's diversity: %s" %max(p))
    return max(p)


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def simpsons_alpha(proportions: list[float]):
    n = sum(proportions)
    s = len(proportions)
    
    p = []
    d = 0.0
    
    for i in proportions:
        p.append(i / n)
        d += i * (i - 1)
    d = d / (n * (n - 1.0))
    # print("Simpson's index of diversity: %s" %(1-d))
    return 1-d


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def inverse_simpsons_alpha(proportions: list[float]):
    n = sum(proportions)
    s = len(proportions)
    
    p = []
    d = 0.0
    
    for i in proportions:
        p.append(i / n)
        d += i * (i - 1)
    d = d / (n * (n - 1.0))
    # print("Simpson's Reciprocal Index: %s" %(1/d))
    return 1/d


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
# def fishers_alpha(abundance_estimates: dict[str, int]):	
#     n = sum(abundance_estimates)
#     s = len(abundance_estimates)
    
#     p = []
#     d = 0.0
    
#     for i in abundance_estimates:
#         p.append(i / n)
#         d += i * (i - 1)
#     d = d / (n * (n - 1.0))
#     global np
#     import numpy as np
#     from scipy.optimize import fsolve
 
#     global N_f
#     N_f = sum(abundance_estimates.values())
#     global S_f
#     S_f = len(abundance_estimates)
    
#     def eqn_output(a):
#         return a * np.log(1+N_f/a) - S_f
    
#     fish = fsolve(eqn_output, 1)
    
#     print("Fisher's index: %s" %fish[0])
#     return fish