import math
import argparse

import MappingUnitData as mu


# Abundance estimate dict should have the format <tax id, count>
def get_alpha_diversity(abundance_estimates: dict):
    n = sum(abundance_estimates.values())
    s = len(abundance_estimates)
    
    p = []
    d = 0
    
    for i in abundance_estimates.values():
        p.append(i / n)
        d += i * (i - 1)
    d = d / (n * (n - 1))
    
    shannons_alpha(p)
    berger_parkers_alpha(p)
    simpsons_alpha(d)
    inverse_simpsons_alpha(d)
    fishers_alpha(abundance_estimates)
    
    
# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def shannons_alpha(p):
    h = []
    for i in p:
        h.append(i * math.log(i))
    print("Shannon's diversity: %s" %(-1 *sum(h)))
    return (-1 *sum(h))


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def berger_parkers_alpha(p):
    print("Berger-parker's diversity: %s" %max(p))
    return max(p)


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def simpsons_alpha(D):
    print("Simpson's index of diversity: %s" %(1-D))
    return 1-D


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def inverse_simpsons_alpha(D):
    print("Simpson's Reciprocal Index: %s" %(1/D))
    return 1/D


# Copied from https://github.com/jenniferlu717/KrakenTools/blob/master/DiversityTools/alpha_diversity.py
def fishers_alpha(abundance_estimates: dict):	
    global np
    import numpy as np
    from scipy.optimize import fsolve
 
    global N_f
    N_f = sum(abundance_estimates.values())
    global S_f
    S_f = len(abundance_estimates)
    
    def eqn_output(a):
        return a * np.log(1+N_f/a) - S_f
    
    fish = fsolve(eqn_output, 1)
    
    print("Fisher's index: %s" %fish[0])
    return fish


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Abundance Estimates", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("classification_file_prefix", help="Prefix of the classification file")
    parser.add_argument("-f", "--min-frequency", type=float, help="Minimum frequency of taxon label to plot", default=0.0)
    args = parser.parse_args()
    config = vars(args)
    min_freq = config["min_frequency"]
    file_prefix = config["classification_file_prefix"]
    
    mapping_units = mu.MappingUnitData(file_prefix)
    mapping_units.filter_by_frequency(min_freq)
    abundance_estimates = mapping_units.get_abundance_estimates()
    
    get_alpha_diversity(abundance_estimates)