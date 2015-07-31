"""Examine correlation of diversity across taxa and across continents using IUCN data, 
to see if what Jetz and Rahbek found always hold."""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import toy_iucn as ti
import psycopg2
from osgeo import ogr
import cPickle
import numpy as np
import glob
import shapely
from scipy.stats.stats import pearsonr

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
        taxon_cont_range_dic[sp] = sp_range.area
   
    out_file = open(proj_dir + 'emp_range_corr\\' + taxon + '_' + continent + '_range.pkl', 'wb')
    cPickle.dump(taxon_cont_range_dic, out_file, protocol = 2)
    out_file.close()
   
   # Rank species based on their range sizes on the continent 
    sp_order_cont = [sp for (area, sp) in sorted(zip(sp_cont_range_list, sp_set))]
    if sp_filter is not 'all': sp_order_cont = sp_order_cont[:int(np.floor(len(sp_order_cont) * 0.75))]
    taxon_cont_richness = [len([x for x in grid if x in sp_order_cont]) for grid in sp_list_flat] 
    orders = [sp_order_cont, sp_order_cont[::-1]] # 07/29/15: Only look at continental range distributions
    
    sp_accumu_ind = np.zeros((2, len(sp_list_flat)))
    r_ind = np.zeros((2, len(sp_set)))
    
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
    
if __name__ == '__main__':        
    taxon_list = ['amphibians', 'reptiles', 'birds', 'terrestrial_mammals']
    # For simplicity, use default continent in shp file
    continent_list = ['Asia', 'North America', 'Europe', 'Africa', 'South America']
    for taxon in taxon_list:
        for continent in continent_list:
            corr_richness_taxon_continent(taxon, continent)
            corr_richness_taxon_continent(taxon, continent, sp_filter = 'lower')
            
    # Plot results
    #ranking_list = ['continent']
    #out_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\emp_r\\'
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
            