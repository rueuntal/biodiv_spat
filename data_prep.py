"""General data preparation"""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\biodiv_spat')
import spatial_functions as spat
import cPickle as pickle
import shapely.wkt, shapely.ops
import numpy as np
import csv
import gdal
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
STAP = SignatureTranslatedAnonymousPackage

mammal_rangemap_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\TERRESTRIAL_MAMMALS\\TERRESTRIAL_MAMMALS.shp'
amphibian_rangemap_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\AMPHIBIANS\\AMPHIBIANS.shp'
bien_rangemaps_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\BIEN_TREES\\'
out_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\final_results\\'
env_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\pred_vars_raw\\'
continent_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\continents\\continent.shp'
#taxa = ['amphibians', 'mammals', 'plants']
# Define R functions
with open('C:\\Users\\Xiao\\Documents\\GitHub\\biodiv_spat\\rfunctions.R', 'r') as f:
    string = f.read()
rfunctions = STAP(string, 'rfunctions')     

## 1. Reproject species range maps and convert their presence into an array
#for i, rangemap_dir in enumerate([amphibian_rangemap_dir, mammal_rangemap_dir, bien_rangemaps_folder]):
    #if i is not 2:
        #range_geom_dic = spat.import_shapefile_field(rangemap_dir, "'presence'", "(1)")
    #else:
        #range_geom_dic = spat.import_shapefile_folder(rangemap_dir)
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
taxa = ['plants']
for taxon in taxa:
    range_geom_dic = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_range_wkt.pkl')
    global_range = spat.create_sp_range_dic(range_geom_dic)
    for sp in global_range.keys():
        global_range[sp] = [global_range[sp]]
    for continent in continent_list:
        if (taxon is 'plants') and (continent in ['Africa', 'Eurasia']): # plant data only available in the new world
            for sp in global_range.keys():
                global_range[sp].append(0)
        else:
            cont_geom = continents_geom[continent]
            cont_range = spat.create_sp_range_dic(range_geom_dic, cont_geom)
            for sp in global_range.keys():
                global_range[sp].append(cont_range[sp])
    
    out_dir = out_folder + 'sp_dist\\' + taxon + '_range_size.csv'
    out_write = open(out_dir, 'wb')
    out = csv.writer(out_write)
    out.writerow(['sp', 'global'] + continent_list)
    for sp in sorted(global_range.keys()):
        row = [sp] + global_range[sp]
        out.writerow(row)
    out_write.close()

## 3. Biodiversity rasters
#for taxon in taxa: 
    #taxon_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    #spat.richness_to_raster(taxon_array, out_folder + 'sp_dist\\' + taxon + '_biodiv.tif')
## 3. Reproject environmental layers
## 3a. AET 
#in_file_AET = env_folder + 'AET_YR\\aet_yr\\w001001.adf'
#out_file_AET = out_folder + 'env\\AET.tif'
#spat.reproj_raster_to_match(in_file_AET, out_file_AET, out_folder + 'sp_dist\\mammals_biodiv.tif')
## 3b. temporal variation in mean annual T (bio4)
#in_file_tempT = env_folder + 'bio_10m_bil\\bio4.bil'
#out_file_tempT = out_folder + 'env\\tempT.tif'
#spat.reproj_raster_to_match(in_file_tempT, out_file_tempT, out_folder + 'sp_dist\\mammals_biodiv.tif')
## 3c. spatial variation in mean annual T (standard deviation in 30s bio1, 20*20 = 10m)
#in_file_30sT = env_folder + 'bio_30s_meanT\\bio_1.bil'
#out_file_Tsd_WGS = out_folder + 'env\\spatT_WGS.tif'
#out_file_Tsd = out_folder + 'env\\spatT.tif'
#rfunctions.rastersd(in_file_30sT, out_file_Tsd)
#spat.reproj_raster_to_match(out_file_Tsd_WGS, out_file_Tsd, out_folder + 'sp_dist\\mammals_biodiv.tif')


