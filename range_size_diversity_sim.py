from __future__ import division
import numpy as np
import scipy.stats
import copy
from scipy.stats.stats import pearsonr

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

def convert_2D_to_1D(i, j, width):
    """Convert [i, j] in a grid with known width and height to 1-D index."""
    return j * width + i

def convert_1D_to_2D(index, width):
    """Convert a 1-D index back to [i, j] in grid."""
    j = int(np.floor(index / width))
    i = int(index - j * width)
    return i, j

def new_edge_cells(i, j, width, height):
    """Given the position of the focal cell [i, j], return the location of its 4 surrouding cells, 
    
    taking into account grid boundaries.
    
    """
    edge_list = []
    if i >= 1: edge_list.append([i - 1, j])
    if i <= width - 2: edge_list.append([i + 1, j])
    if j >= 1: edge_list.append([i, j - 1])
    if j <= height - 2: edge_list.append([i, j + 1])
    return edge_list

def ind_range_continuous(width, height, Nsize, env_p = None):
    """Subfunction to generate a continuous species range.
    
    Inputs:
    width, height - dimension of the landscape.
    Nsize - size of the focal range, in number of grid cells.
    env_p - 1-D list of probability of occupancy. If None, all cells are equally likely.
    
    """
    np.random.seed()
    range_set = set([])
    edge_set = set([])
    loc = np.random.choice(range(width * height), size = 1, p = env_p)
    range_set.update(loc)
    while len(range_set) < Nsize:
        i, j = convert_1D_to_2D(loc, width)
        new_edges = new_edge_cells(i, j, width, height)
        new_edges_1D = [convert_2D_to_1D(cell[0], cell[1], width) for cell in new_edges]
        new_edges_1D = [cell for cell in new_edges_1D if cell not in range_set]
        edge_set.update(new_edges_1D)
        # Obtain p for edge cells
        edge_list = list(edge_set)
        if env_p: 
            p_edge = [env_p[index] for index in edge_list]
            p_edge = np.array(p_edge)  / np.sum(p_edge)
        else: p_edge = None
        loc = np.random.choice(edge_list, size = 1, p = p_edge)
        edge_set.remove(int(loc))
        range_set.update(loc)
    
    return np.array(list(range_set))
        
def ind_range_generator(width, height, size, continuous = False, env = 0, env_landscape = None, r = None):
    """Generate the spatial range of one species on the landscale.
    
    Inputs:
    width, height - dimension of the landscape.
    size - range size of the focal species.
    continuous - whether the range of each species is continuous or scattered. 
    env - whether the range correlates with an environmental variable.
    env_landscape - the landscape of the environmental variable. Only needed when env = True.
    r - correlation coefficient between (transformed) probability of species presence and
        environmental variable. Only needed when env = True.
    
    Output:
    spatial_range - an array of shape (height, width) where True is for presence and False is for absence.
    """
    np.random.seed()
    spatial_range = np.empty([height, width], dtype = bool)
    spatial_range.fill(False)
    size_grid = int(np.ceil(size))
    if env: # Convert the environment landscape to probability landscale
        env_flat = np.array([env_landscape[j][i] for j in range(height) for i in range(width)])
        noise_term = scipy.stats.norm.rvs(size = len(env_flat))
        noise_term = (noise_term - np.mean(noise_term)) / np.std(noise_term, ddof = 1)
        a = r / np.sqrt(1 - r **2)
        env_flat_standardize = (env_flat - np.mean(env_flat)) / np.std(env_flat, ddof = 1)
        p_transform = a * env_flat_standardize + noise_term
        p_raw = 1 / (1 + np.exp(-p_transform)) # Use logistic function to convert transformed env to p
        p = p_raw / np.sum(p_raw)
    else: p = None
    if not continuous:
        loc_sp = np.random.choice(range(width * height), size = size_grid, replace = False, p = p)
    else:
        loc_sp = ind_range_continuous(width, height, size_grid, env_p = p)
    for loc in loc_sp:
        i, j = convert_1D_to_2D(loc, width)
        spatial_range[j][i] = True
    return spatial_range

def sim_range_size_landscape(width, height, mu, sigma, S, continuous = False, env = 0, r = 0):
    """Simulate species range sizes on a landscape and explore their correlations with overall diversity.
    
    Input:
    width, height - dimension of the landscape (grid cells).
    mu, sigma - parameters for the Poisson lognormal distribution from which the global range size distribution is simulated.
    S - total number of species
    continuous - whether the range of each species is continuous or scattered. If True, the range is in one polygon. If False,
        each cell is randomly assigned.
    env - number of environmental factors correlated with spatial distribution of ranges. 0 (default) means that the ranges
        do not correlate with any environmental factors and are randomly distributed. 1 means that all ranges are driven
        by the same factor. 2 means that 1/4 of the species having the largest ranges are driven by one factor, while the other
        species are driven by another. 
    r - correlation coefficient between (transformed) probability of species presence and
        environmental variable. Only needed when env = True. Called in subfunction ind_range_generator() when there is 
        environmental gradient.
    
    Output: 
    r_quartile -  a list with 4 correation coefficient r showing correlation between richness in each range size quartile and overall 
        richness, from the smallest species to the largest species.
    r_low - a list of length S, with r between total diversity and cumulative diversity from the smallest-ranged species to
        the largest-ranged species.
    r_high - a list of length S, with r between total diversity and cumulative diversity from the largest-ranged species to
        the smallest-ranged species.
    
    """
    sp_range_list = np.sort(global_range_size(mu, sigma, S, height * width))
    sp_range_array_list = []
    richness_landscape = np.empty([height, width], dtype = int)
    richness_landscape.fill(0)
    richness_from_low = copy.deepcopy(richness_landscape)
    richness_from_high = copy.deepcopy(richness_landscape)
    richness_quar1 = copy.deepcopy(richness_landscape)
    richness_quar2 = copy.deepcopy(richness_landscape)
    richness_quar3 = copy.deepcopy(richness_landscape)
    richness_quar4 = copy.deepcopy(richness_landscape)
    r_low, r_high, r_quartile = [], [], []
    
    if env == 0: env_list = [0, 0]
    # elif env == 1: generate landscape, env_list = [landscape, landscape]
    # elif env == 2: env_list = [generate landscape 1, generate landscape 2]
    
    for i, size in enumerate(sp_range_list):
        if i < S * 0.75: sp_range_landscape = ind_range_generator(width, height, size, continuous = continuous, \
                                                                  env = env, env_landscape = env_list[0], r = r)
        else: sp_range_landscape = ind_range_generator(width, height, size, continuous = continuous, \
                                                                env = env, env_landscape = env_list[1], r = r)
        sp_range_array_list.append(sp_range_landscape)
        richness_landscape += sp_range_landscape
    
    # Analysis
    for i in range(S):
        richness_from_low += sp_range_array_list[i]
        r_low.append(pearsonr(np.ravel(richness_from_low), np.ravel(richness_landscape))[0])
        richness_from_high += sp_range_array_list[S - i - 1]
        r_high.append(pearsonr(np.ravel(richness_from_high), np.ravel(richness_landscape))[0])
        if i < S * 0.25: richness_quar1 += sp_range_array_list[i]
        elif i < S * 0.5: richness_quar2 += sp_range_array_list[i]
        elif i < S * 0.75: richness_quar3 += sp_range_array_list[i]
        else: richness_quar4 += sp_range_array_list[i]
    
    for quar in [richness_quar1, richness_quar2, richness_quar3, richness_quar4]:
        r_quartile.append(pearsonr(np.ravel(quar), np.ravel(richness_landscape))[0])
    
    return r_quartile, r_low, r_high

def sim_range_size_landscape_Niter(width, height, mu, sigma, S, Niter, out_dir, out_name, 
                                   continuous = False, env = 0, r = 0):
    """Run sim_range_size_landscape_Niter multiple times and save the three output lists into txt files 
    
    with names out_dir + out_name + '_quartile.txt'/'_low.txt'/'_high.txt'. 
    
    """
    out_file_extension = ['_quartile.txt', '_low.txt', '_high.txt']
    for i in range(Niter):
        three_rs = sim_range_size_landscape(width, height, mu, sigma, S, 
                                                             continuous = continuous, env = env, r = r)
        for j in range(len(out_file_extension)):
            out_file = out_dir + '/' + out_name + out_file_extension[j]
            out_file_write = open(out_file, 'a')
            r_list = three_rs[j]
            print>>out_file_write, '\t'.join([str(round(x, 5)) for x in r_list])
            out_file_write.close()
            