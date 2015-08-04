"""Examine correlation of diversity across taxa and across continents using IUCN data, 
to see if what Jetz and Rahbek found always hold."""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import toy_iucn as ti
import psycopg2
from osgeo import ogr
import cPickle
import numpy as np
import glob
import shapely
from scipy.stats.stats import pearsonr
import os.path
import copy
import range_size_diversity_sim as rsim

def corr_richness_taxon_continent(taxon, continent, sp_filter = 'all'):
    """Obtain the correlation between overall richness and partial richness (by adding species one by one)
    
    for one taxon on a given continent.
    Inputs:
    taxon - string, name of the taxon.
    continent - string, name of the continent. Can take one of seven values: Asia, North America,
        Europe, Africa, South  America, Oceania, Australia.
    sp_filter - species used in the analysis. If 'all', all species are included. If 'lower', only species with the lower 3/4 ranges
        on the given continent are included. 
        
    Outputs:
    Two txt files saved to disc. 
    ind_sp_corr.txt or ind_sp_corr_lower.txt: r versus species accumulation. 
        Columns: taxon, continent, global/continent, high/low, r from species 1 to species S.
    quart_corr.txt or quart_corr_lower.txt: r between overall richness and richness of each quartile of range sizes. 
        Columns: taxon, continent, global/continent, four or three r values.
    """
    proj_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\'
    pixel_size = 100000
    
    if taxon is not 'birds': 
        tables = psycopg2.connect("dbname = 'IUCN_range_maps' port = 5432 user = 'postgres' password = 'lavender' host = 'localhost'")
        postgis_cur = tables.cursor()
    
    sp_list_array = ti.import_pickle_file(proj_dir + 'IUCN_sp_lists\\' + taxon + '_100000.pkl')
    conts_dir = proj_dir + 'continents\\continent.shp'
    cont_geom  = ti.reproj_geom(ti.import_shapefile(conts_dir, Attr = 'CONTINENT', AttrFilter = continent)[0])
    cont_array = ti.create_array_for_raster(ti.proj_extent('behrmann'), geom = cont_geom, pixel_size = pixel_size) 
    sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if cont_array[i][j] == 1]
    sp_list_flat = [grid for grid in sp_list_flat if len(grid) > 0] # Remove empty grid cells with no species from the given taxon.
    sp_set = list(set([sp for grid in sp_list_flat for sp in grid]))
    
    # Obtain range size on continent and save to disk
    dic_dir = proj_dir + 'emp_range_corr\\' + taxon + '_' + continent + '_range.pkl'
    if os.path.isfile(dic_dir):  # So that range sizes on continents don't have to be computed over and over again
        taxon_cont_range_dic = ti.import_pickle_file(dic_dir)
        sp_cont_range_list = [taxon_cont_range_dic[sp] for sp in sp_set]
    else:
        taxon_cont_range_dic = {}    
        cont_shape = shapely.wkt.loads(cont_geom)
        sp_cont_range_list = []
        for sp in sp_set: 
            if taxon is not 'birds':
                sp_range = shapely.wkt.loads(ti.sp_reproj(postgis_cur, taxon, sp))
            else: 
                genus, spname = sp.split(' ')
                sp_shp = proj_dir + 'IUCN_range_maps\\BIRDS\\' + genus + '_' + spname + '*.shp'
                sp_dir = glob.glob(sp_shp)[0]
                sp_geom_list = ti.import_shapefile(sp_dir)
                sp_geom_shapes = [shapely.wkt.loads(x) for x in sp_geom_list]
                try:
                    sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
                except: 
                    sp_geom_shapes = [x.buffer(0) for x in sp_geom_shapes]
                    sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))            
                sp_range = shapely.wkt.loads(ti.reproj_geom(sp_geom_wkt))
            
            # Try shapely instead of ogr
            try: sp_cont_range = (sp_range.intersection(cont_shape)).area
            except:
                sp_range = sp_range.buffer(0)
                sp_cont_range = (sp_range.intersection(cont_shape)).area
            sp_cont_range_list.append(sp_cont_range)
            taxon_cont_range_dic[sp] = sp_cont_range
       
        out_file = open(dic_dir, 'wb')
        cPickle.dump(taxon_cont_range_dic, out_file, protocol = 2)
        out_file.close()
       
   # Rank species based on their range sizes on the continent 
    sp_order_cont = [sp for (area, sp) in sorted(zip(sp_cont_range_list, sp_set))]
    if sp_filter is not 'all': sp_order_cont = sp_order_cont[:int(np.floor(len(sp_order_cont) * 0.75))]
    taxon_cont_richness = [len([x for x in grid if x in sp_order_cont]) for grid in sp_list_flat] 
    orders = [sp_order_cont, sp_order_cont[::-1]] # 07/29/15: Only look at continental range distributions
    
    sp_accumu_ind = np.zeros((2, len(sp_list_flat)))
    r_ind = np.zeros((2, len(sp_order_cont)))
    
    if sp_filter is 'all': 
        sp_accumu_quart_cont = np.zeros((4, len(sp_list_flat)))
        r_quart = np.zeros((1, 4))
    else: 
        sp_accumu_quart_cont = np.zeros((3, len(sp_list_flat)))
        r_quart = np.zeros((1, 3))        
    
    for j in range(len(sp_order_cont)): # Loop through species
        for i in range(2): # Loop through two ways to accumulate species - from the lower end, or from the higher end
            sp = orders[i][j]
            sp_dist = np.array([1 if sp in grid else 0 for grid in sp_list_flat])
            sp_accumu_ind[i] += sp_dist
            r_ind[i][j] = pearsonr(sp_accumu_ind[i], taxon_cont_richness)[0]
            if i == 0: # If the order is low to high, continent
                sp_accumu_quart_cont[np.floor(j / len(sp_set) * 4)] += sp_dist
    
    for i in range(len(sp_accumu_quart_cont)):
        r_quart[0][i] = pearsonr(sp_accumu_quart_cont[i], taxon_cont_richness)[0]
        
    # Save output 
    ind_header = [['continent', 'low'], ['continent', 'high']]
    if sp_filter is 'all': 
        out_ind = open(proj_dir + 'emp_range_corr\\ind_sp_corr.txt', 'a')
        out_quart = open(proj_dir + 'emp_range_corr\\quart_corr.txt', 'a')
    else:
        out_ind = open(proj_dir + 'emp_range_corr\\ind_sp_corr_lower.txt', 'a')
        out_quart = open(proj_dir + 'emp_range_corr\\quart_corr_lower.txt', 'a')
        
    for i, r_ind_row in enumerate(r_ind):
        out_row = [taxon, continent] + ind_header[i] + list(r_ind_row)
        print>>out_ind, '\t'.join(map(str, out_row))
    for j, r_quart_row in enumerate(r_quart):
        out_row_quart = [taxon, continent] + ['continent'] + list(r_quart_row)
        print>>out_quart, '\t'.join(map(str, out_row_quart))
    out_ind.close()
    out_quart.close()

def sim_taxon_continent_landscape(taxon, continent, continuous = False):
    """This function is modified from sim_range_size_landscale() from module range_size_diversity_sim.
    
    Output: 
    r_quartile_r4, r_quantile_r3 -  lists with 4 or 3 correation coefficient r showing correlation between richness in each range size quartile and overall 
        richness for all species or the lower 3/4 species.
    r_low_r4, r_low_r3 - lists of length S or 3/4 * S, with r between total diversity and cumulative diversity from the smallest-ranged species to
        the largest-ranged species.
    r_high_r4, r_low_r3 - lists of length S or 3/4 * S, with r between total diversity and cumulative diversity from the largest-ranged species to
        the smallest-ranged species.
    
    """
    proj_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\'
    pixel_size = 100000
    
    sp_list_array = ti.import_pickle_file(proj_dir + 'IUCN_sp_lists\\' + taxon + '_100000.pkl')
    conts_dir = proj_dir + 'continents\\continent.shp'
    cont_geom  = ti.reproj_geom(ti.import_shapefile(conts_dir, Attr = 'CONTINENT', AttrFilter = continent)[0])
    cont_array = ti.create_array_for_raster(ti.proj_extent('behrmann'), geom = cont_geom, pixel_size = pixel_size) 
    sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if cont_array[i][j] == 1]
    sp_list_flat = [grid for grid in sp_list_flat if len(grid) > 0] # Remove empty grid cells with no species from the given taxon.
    # Read in list of species range sizes on continent
    dic_dir = proj_dir + 'emp_range_corr\\' + taxon + '_' + continent + '_range.pkl'
    taxon_cont_range_dic = ti.import_pickle_file(dic_dir)
    sp_range_list = sorted(taxon_cont_range_dic.values())
    
    sp_range_array_list = []
    sim_len = int(np.ceil(np.sqrt(len(sp_list_flat)))) # For simplicity, simulate a square grid 
    richness_landscape = np.empty([sim_len, sim_len], dtype = int)
    richness_landscape.fill(0)
    richness_from_low = copy.deepcopy(richness_landscape)
    richness_from_high_4 = copy.deepcopy(richness_landscape)
    richness_from_high_3 = copy.deepcopy(richness_landscape)
    richness_quar1 = copy.deepcopy(richness_landscape)
    richness_quar2 = copy.deepcopy(richness_landscape)
    richness_quar3 = copy.deepcopy(richness_landscape)
    richness_quar4 = copy.deepcopy(richness_landscape)
    r_low_4, r_low_3, r_high_4, r_high_3, r_quartile_4, r_quartile_3 = [], [], [], [], [], []
    
    for i, size in enumerate(sp_range_list):
        sp_range_landscape = rsim.ind_range_generator(sim_len, sim_len, size / pixel_size ** 2, continuous = continuous, 
                                                      env = 0, env_landscape = 0, r = 1)
        sp_range_array_list.append(sp_range_landscape)
        richness_landscape += sp_range_landscape
        if i == int(np.floor(len(sp_range_list) * 3/4)): richness_landscape_3 = copy.deepcopy(richness_landscape)
        
    # Analysis
    S = len(sp_range_list)
    for i in range(S):
        richness_from_low += sp_range_array_list[i]
        r_low_4.append(pearsonr(np.ravel(richness_from_low), np.ravel(richness_landscape))[0])
        richness_from_high_4 += sp_range_array_list[S - i - 1]
        r_high_4.append(pearsonr(np.ravel(richness_from_high_4), np.ravel(richness_landscape))[0])
        if i <= int(np.floor(len(sp_range_list) * 3/4)):
            richness_from_high_3 += sp_range_array_list[int(np.floor(len(sp_range_list) * 3/4)) - i]
            r_high_3.append(pearsonr(np.ravel(richness_from_high_3), np.ravel(richness_landscape_3))[0])
            r_low_3.append(pearsonr(np.ravel(richness_from_low), np.ravel(richness_landscape_3))[0])
        if i < S * 0.25: richness_quar1 += sp_range_array_list[i]
        elif i < S * 0.5: richness_quar2 += sp_range_array_list[i]
        elif i < S * 0.75: richness_quar3 += sp_range_array_list[i]
        else: richness_quar4 += sp_range_array_list[i]
    
    for quar in [richness_quar1, richness_quar2, richness_quar3, richness_quar4]:
        r_quartile_4.append(pearsonr(np.ravel(quar), np.ravel(richness_landscape))[0])
        r_quartile_3.append(pearsonr(np.ravel(quar), np.ravel(richness_landscape_3))[0])
    del r_quartile_3[-1]
    
    return r_quartile_4, r_quartile_3, r_low_4, r_low_3, r_high_4, r_high_3

def sim_range_size_landscape_Niter(taxon, continent, Niter, out_dir, out_name, continuous = False):
    """Run sim_taxon_continent_landscape multiple times and save the three output lists into txt files 
    
    with names out_dir + out_name + '_quartile(_lower).txt'/'_low(_lower).txt'/'_high(_lower).txt'. 
    
    """
    out_file_extension = ['_quartile.txt', '_quartile_lower.txt', '_low.txt', '_low_lower.txt', '_high.txt', '_high_lower.txt']
    for i in range(Niter):
        six_rs = sim_taxon_continent_landscape(taxon, continent, continuous = continuous)
        for j in range(len(out_file_extension)):
            out_file = out_dir + '\\' + out_name + out_file_extension[j]
            out_file_write = open(out_file, 'a')
            r_list = six_rs[j]
            print>>out_file_write, '\t'.join([str(round(x, 5)) for x in r_list])
            out_file_write.close()
            
def import_ind_rows(taxon, cont, rank, file_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\emp_range_corr\\ind_sp_corr.txt'):
    """Return the r^2 values from ind_sp_corr.txt with given taxon, continent, and ranking method."""
    with open(file_dir) as f:
        content = f.read().splitlines()
    content = [x.split('\t') for x in content]
    out = []
    for x in content:
        if x[0:3] == [taxon, cont, rank]:
            out.append([float(y) for y in x[4:]])
    return out

def import_quart_file(taxon, cont, rank, file_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\emp_range_corr\\quart_corr.txt'):
    """Return the r^2 values from quart_corr.txt with given taxon, continent, and ranking method."""
    with open(file_dir) as f:
        content = f.read().splitlines()
    content = [x.split('\t') for x in content]
    for x in content:
        if x[0:3] == [taxon, cont, rank]:
            return [float(y) for y in x[3:]]
    
def plot_quartile(list_of_r, ax):
    """Plot the r^2 of the 4 quartiles."""
    plt.ylim(-0.5, 1)
    ax.plot(range(1, 5), list_of_r, 'bo-')
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.locator_params(nbins=5)
    ax.set_xlabel('Quantile', labelpad = 4, size = 8)
    ax.set_ylabel('r value', labelpad = 4, size = 8)
    return ax

def plot_quartile_comp(list_of_r_4, list_of_r_3, ax):
    """Compare the r-value obtained for quantiles when all species are included, 
    
    or when the species with the largest ranges are excluded.
    
    """
    plt.ylim(-0.5, 1.5)
    ax.plot(range(1, 5), list_of_r_4, 'bo-')
    ax.plot(range(1, 4), list_of_r_3, 'ro-')
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.locator_params(nbins=5)
    ax.set_xlabel('Quantile', labelpad = 4, size = 8)
    ax.set_ylabel('r value', labelpad = 4, size = 8)
    return ax

def plot_ind_accum(list_of_lists_of_r, ax):
    """Plot the increase of r^2 with the increase of number of species, from both ends."""
    plt.ylim(-0.5, 1)
    ax.plot(range(len(list_of_lists_of_r[0])), list_of_lists_of_r[0], color = 'red', linewidth = 2)
    ax.plot(range(len(list_of_lists_of_r[1])), list_of_lists_of_r[1], color = 'blue', linewidth = 2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.locator_params(nbins=5)
    ax.set_xlabel('Number of species', labelpad = 4, size = 8)
    ax.set_ylabel('r value', labelpad = 4, size = 8)
    return ax
 
def plot_ind_accum_comp(list_of_lists_of_r_4, list_of_lists_of_r_3, ax):
    """Compare the increase of r^2 with the increase of number of species, from both ends, 
    
    when the most broad-ranged species are included versus excluded.
    
    """
    plt.ylim(-0.5, 1.5)
    ax.plot(range(len(list_of_lists_of_r_4[0])), list_of_lists_of_r_4[0], color = 'blue', linewidth = 2)
    ax.plot(range(len(list_of_lists_of_r_4[1])), list_of_lists_of_r_4[1], color = 'blue', linewidth = 2, linestyle = '--')
    ax.plot(range(len(list_of_lists_of_r_3[0])), list_of_lists_of_r_3[0], color = 'red', linewidth = 2)
    ax.plot(range(len(list_of_lists_of_r_3[1])), list_of_lists_of_r_3[1], color = 'red', linewidth = 2, linestyle = '--')    
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    ax.locator_params(nbins=5)
    ax.set_xlabel('Number of species', labelpad = 4, size = 8)
    ax.set_ylabel('r value', labelpad = 4, size = 8)
    return ax
    
if __name__ == '__main__':        
    taxon_list = ['amphibians', 'reptiles', 'birds', 'terrestrial_mammals']
    # For simplicity, use default continent in shp file
    continent_list = ['Asia', 'North America', 'Europe', 'Africa', 'South America']
    for taxon in taxon_list:
        for continent in continent_list:
            corr_richness_taxon_continent(taxon, continent)
            corr_richness_taxon_continent(taxon, continent, sp_filter = 'lower')
   
    ## Run simulations to set the expected bounds for empirical results
    #out_dir_sim = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\emp_range_corr\\sim'
    #for taxon in taxon_list:
        #for continent in continent_list:
            #out_name_scatter = taxon + '_' + continent + '_scattered'
            #sim_range_size_landscape_Niter(taxon, continent, 1000, out_dir_sim, out_name_scatter)
            #out_name_continuous = taxon + '_' + continent + '_continuous'
            #sim_range_size_landscape_Niter(taxon, continent, 1000, out_dir_sim, out_name_continuous, continuous = True)
            
    ## Plot comparison
    #out_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\emp_r\\'
    #ranking = 'continent'
    #for taxon in taxon_list:
        #fig = plt.figure(figsize = (4, 10))
        #for i, continent in enumerate(continent_list):
            #quartile_r2_list_r4 = import_quart_file(taxon, continent, ranking)
            #quartile_r2_list_r3 = import_quart_file(taxon, continent, ranking, 
                                                    #file_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\emp_range_corr\\quart_corr_lower.txt')
            #ind_r2_list_r4 = import_ind_rows(taxon, continent, ranking)
            #ind_r2_list_r3 = import_ind_rows(taxon, continent, ranking, 
                                             #file_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\emp_range_corr\\ind_sp_corr_lower.txt')
            #ax1 = plt.subplot(7, 2, 2 * i + 1)
            #sub = plot_quartile_comp(quartile_r2_list_r4, quartile_r2_list_r3, ax1)
            #plt.title(continent, size = 16)
            #ax2 = plt.subplot(7, 2, 2 * i + 2)
            #sub = plot_ind_accum_comp(ind_r2_list_r4, ind_r2_list_r3, ax2)
            #if i == 0:
                #r4_label = mlines.Line2D([], [], color = 'blue', linestyle = '-', linewidth = 1, marker = 'o', 
                                     #label = 'All species')
                #r3_label = mlines.Line2D([], [], color = 'red', linestyle = '-', linewidth = 1, marker = 'o', 
                                                     #label = '3/4 species')                                
                #ax1.legend(loc = 'lower right', handles = [r4_label, r3_label], prop={'size':6})
                #high_label = mlines.Line2D([], [], color = 'black', linestyle = '-', linewidth = 0.5, label = 'From high')
                #low_label = mlines.Line2D([], [], color = 'black', linestyle = '--', linewidth = 0.5, label = 'From low')
                #ax2.legend(loc = 'lower right', handles = [high_label, low_label], prop={'size':6}, handlelength=3)
                
        #plt.subplots_adjust(top = 0.95, bottom = -0.3,  wspace = 0.4, hspace = 0.65)       
        #out_name = out_dir + taxon + '_' + ranking + '_comp.pdf'
        #plt.savefig(out_name, format = 'pdf', dpi = 400)
        #plt.close(fig)            
            
    # Plot results
    #ranking_list = ['continent'] 
    #for ranking in ranking_list:
        #for taxon in taxon_list:
            #fig = plt.figure(figsize = (4, 14))
            #for i, continent in enumerate(continent_list):
                #quartile_r2_list = import_quart_file(taxon, continent, ranking)
                #ind_r2_list = import_ind_rows(taxon, continent, ranking)
                #ax1 = plt.subplot(7, 2, 2 * i + 1)
                #sub = plot_quartile(quartile_r2_list, ax1)
                #plt.title(continent, size = 16)
                #ax2 = plt.subplot(7, 2, 2 * i + 2)
                #sub = plot_ind_accum(ind_r2_list, ax2)
            #plt.subplots_adjust(wspace = 0.4, hspace = 0.55)
            #out_name = out_dir + taxon + '_' + ranking + '.pdf'
            #plt.savefig(out_name, format = 'pdf', dpi = 400)
            #plt.close(fig)
            