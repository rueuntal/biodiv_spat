from __future__ import division
import numpy as np
import scipy.stats

def global_range_size(mu, sigma, S, max_size):
    """Generate a global range size distribution with S species from a
    
    lognormal distribution, with an upper limit. 
    
    """
    sizes = scipy.stats.lognorm.rvs(sigma, scale = np.exp(mu), size = S)
    sizes = sizes[sizes <= max_size]
    while len(sizes) < S:
        additional_ranges = scipy.stats.lognorm.rvs(sigma, scale = np.exp(mu), size = S - len(sizes))
        additional_ranges = additional_ranges[additional_ranges <= max_size]
        sizes = np.append(sizes, additional_ranges)
    return sizes

def sim_range_size_landscape(width, height, mu, sigma, S, Nsim):
    """Simulate species range sizes on a landscape and explore their correlations with overall diversity.
    
    Input:
    width, height - dimension of the grid.
    mu, sigma - parameters for the Poisson lognormal distribution from which the global range size distribution is simulated.
    S - total number of species
    
    Output: 
    r2_quartile -  a list with 4 R^2 values showing correlation between richness in each range size quartile and overall 
        richness, from the smallest species to the largest species.
    
    """