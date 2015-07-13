"""Examine correlation of diversity across taxa and across continents using IUCN data, 
to see if what Jetz and Rahbek found always hold."""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox')

import toy_iucn as ti
import psycopg2
from osgeo import ogr
import numpy as np
import glob
import shapely
from scipy.stats.stats import pearsonr

def corr_richness_taxon_continet(taxon, continent):
    """Obtain the correlation between overall richness and partial richness (by adding species one by one)
    
    for one taxon on a given continent.
    Inputs:
    taxon - string, name of the taxon.
    continent - string, name of the continent. Can take one of seven values: Asia, North America,
        Europe, Africa, South  America, Oceania, Australia.
        
    Outputs:
    
    """
    proj_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\'
    pixel_size = 100000
    
    if taxon is not 'birds': 
        tables = psycopg2.connect("dbname = 'IUCN_range_maps' port = 5432 user = 'postgres' password = 'lavender' host = 'localhost'")
        postgis_cur = tables.cursor()
    
    sp_list_array = ti.import_pickle_file(proj_dir + 'IUCN_sp_lists\\' + taxon + '_100000.pkl')
    sp_global_range = ti.import_pickle_file(proj_dir + 'IUCN_sp_lists\\' + taxon + '_range_size.pkl')
    conts_dir = proj_dir + 'continents\\continent.shp'
    cont_geom  = ti.reproj_geom(ti.import_shapefile(conts_dir, Attr = 'CONTINENT', AttrFilter = continent)[0])
    cont_array = ti.create_array_for_raster(ti.proj_extent('behrmann'), geom = cont_geom, pixel_size = pixel_size) 
    sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if cont_array[i][j] == 1]
    sp_set = list(set([sp for grid in sp_list_flat for sp in grid]))
    sp_global_range_list = [sp_global_range[sp] for sp in sp_set]
    
    # Obtain range size on continent
    cont_shape = ogr.CreateGeometryFromWkt(cont_geom)
    sp_cont_range_list = []
    for sp in sp_set: 
        if taxon is not 'birds':
            sp_range = ogr.CreateGeometryFromWkt(ti.sp_reproj(postgis_cur, taxon, sp))
        else: 
            genus, spname = sp.split(' ')
            sp_shp = proj_dir + 'IUCN_range_maps\\BIRDS\\' + genus + '_' + spname + '*.shp'
            sp_dir = glob.glob(sp_shp)[0]
            sp_geom_list = ti.import_shapefile(in_dir)
            sp_geom_shapes = [shapely.wkt.loads(x) for x in sp_geom_list]
            try:
                sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
            except: 
                sp_geom_shapes = [x.buffer(0) for x in sp_geom_shapes]
                sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))            
            sp_range = ogr.CreateGeometryFromWkt(ti.reproj_geom(sp_geom_wkt))
        try: sp_cont_range = cont_shape.Intersection(sp_range).GetArea()
        except:
            sp_range = sp_range.Buffer(0)
            sp_cont_range = cont_shape.Intersection(sp_range).GetArea()
        sp_cont_range_list.append(sp_cont_range)
    
    taxon_cont_richness = [len(grid) for grid in sp_list_flat]
    sp_order_global = [sp for (area, sp) in sorted(zip(sp_global_range_list, sp_set))]
    sp_order_cont = [sp for (area, sp) in sorted(zip(sp_cont_range_list, sp_set))]
    orders = [sp_order_global, sp_order_global[::-1], sp_order_cont, sp_order_cont[::-1]]
    
    sp_accumu_ind = np.zeros((4, len(sp_list_flat)))
    sp_accumu_quart_cont = np.zeros((4, len(sp_list_flat)))
    sp_accumu_quart_global = np.zeros((4, len(sp_list_flat)))
    
    r2_ind = np.zeros((4, len(sp_set)))
    r2_quart = np.zeros((2, 4))
    
    for j in range(len(sp_set)): # Loop through species
        for i in range(4): # Loop through four different ways to accumulate species
            sp = orders[i][j]
            sp_dist = np.array([1 if sp in grid else 0 for grid in sp_list_flat])
            sp_accumu_ind[i] += sp_dist
            r2_ind[i][j] = pearsonr(sp_accumu_ind[i], taxon_cont_richness)[0]
            if i == 0: # If the order is low to high, global
                sp_accumu_quart_global[np.floor(j / len(sp_set) * 4)] += sp_dist
            elif i == 2: # If the order is low to high, continent
                sp_accumu_quart_cont[np.floor(j / len(sp_set) * 4)] += sp_dist
    
    for i in range(4):
        r2_quart[0][i] = pearsonr(sp_accumu_quart_global[i], taxon_cont_richness)[0]
        r2_quart[1][i] = pearsonr(sp_accumu_quart_cont[i], taxon_cont_richness)[0]
        
    # Save output 
    ind_header = [['global', 'low'], ['global', 'high'], ['continent', 'low'], ['continent', 'high']]
    quart_header = [['global'], ['continent']]
    out_ind = open(proj_dir + 'emp_range_corr\\ind_sp_corr.txt', 'a')
    out_quart = open(proj_dir + 'emp_range_corr\\quart_corr.txt', 'a')
    for i, r2_ind_row in enumerate(r2_ind):
        out_row = [taxon, continent] + ind_header[i] + list(r2_ind_row)
        print>>out_ind, '\t'.join(map(str, out_row))
    for j, r2_quart_row in enumerate(r2_quart):
        out_row_quart = [taxon, continent] + quart_header[j] + list(r2_quart_row)
        print>>out_quart, '\t'.join(map(str, out_row_quart))
    out_ind.close()
    out_quart.close()

if __name__ == '__main__':        
    taxon_list = ['amphibians', 'reptiles', 'birds', 'terrestrial_mammals']
    # For simplicity, use default continent in shp file
    continent_list = ['Asia', 'North America', 'Europe', 'Africa', 'South America', 
                      'Oceania', 'Australia']
    for taxon in taxon_list:
        for continent in continent_list:
            corr_richness_taxon_continet(taxon, continent)