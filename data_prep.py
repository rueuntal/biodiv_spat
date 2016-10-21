"""General data preparation"""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\biodiv_spat')
import spatial_functions as spat
import cPickle as pickle
import shapely.wkt, shapely.ops
import numpy as np
import csv

mammal_rangemap_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\TERRESTRIAL_MAMMALS\\TERRESTRIAL_MAMMALS.shp'
amphibian_rangemap_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\AMPHIBIANS\\AMPHIBIANS.shp'
out_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\final_results\\'
continent_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\continents\\continent.shp'
taxa = ['amphibians', 'mammals']

# 1. Reproject species range maps and convert their presence into an array
#for i, rangemap_dir in enumerate([amphibian_rangemap_dir, mammal_rangemap_dir]):
    #range_geom_dic = spat.import_shapefile_field(rangemap_dir, "'presence'", "(1)")
    ## Save the dictionary with species ranges in wkt as multipolygons
    #pickle.dump(range_geom_dic, open(out_folder + 'sp_dist\\' + taxa[i] + '_range_wkt.pkl', 'wb'))
    #out_dir = out_folder + 'sp_dist\\' + taxa[i] + '_dist.pkl'
    #spat.create_array_sp_list(range_geom_dic, out_dir)

# 2. Species range size (global & by continent)
continents_geom = spat.import_shapefile_field(continent_dir, field = 'CONTINENT')
asia, europe = continents_geom['Asia'], continents_geom['Europe']
eurasia = shapely.ops.cascaded_union([shapely.wkt.loads(asia), shapely.wkt.loads(europe)])
eurasia_wkt = shapely.wkt.dumps(eurasia)
continents_geom['Eurasia'] = eurasia_wkt
continent_list = ['Africa', 'North America', 'South America', 'Eurasia']
for taxon in taxa:
    range_geom_dic = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_range_wkt.pkl')
    global_range = spat.create_sp_range_dic(range_geom_dic)
    for sp in global_range.keys():
        global_range[sp] = [global_range[sp]]
    for continent in continent_list:
        cont_geom = continents_geom[continent]
        cont_range = spat.create_sp_range_dic(range_geom_dic, cont_geom)
        for sp in global_range.keys():
            global_range[sp].append(cont_range[sp])
    
    out_dir = out_folder + 'sp_dist\\' + taxon + '_range_size.csv'
    out_write = open(out_dir, 'wb')
    out = csv.writer(out_write)
    out.writerow(['sp', 'global'] + continent_list)
    for sp in sorted(global_range_keys()):
        row = [sp] + global_range[sp]
        out.writerow(row)
    out_write.close()
