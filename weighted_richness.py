"""Script to examine how the importance of environmental variables shifts with weights on species range sizes."""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox')

import spatial_functions as spat
import psycopg2
import numpy as np

if __name__ == '__main__':
    # 1. Map species ranges onto 100km*100km grid cells, and obtain dictionaries of range sizes
    # Range maps of terrestrial mammals, birds, and amphibians are obtained from IUCN
    # Projection: Behrmann equal area
    taxa = ['terrestrial_mammals', 'amphibians']
    con_ranges = psycopg2.connect("dbname = 'IUCN_range_maps' port = 5432 user = 'postgres' password = 'lavender' host = 'localhost'")
    postgis_cur = con_ranges.cursor()
    pixel_size = 100000
    out_folder_sp_dist = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\sp_dist\\'
    
    for taxon in taxa:
        spat.create_array_sp_list(postgis_cur, taxon, out_folder_sp_dist + taxon + '_' + str(pixel_size) + '.pkl', pixel_size = pixel_size)
        spat.create_sp_range_dic(postgis_cur, taxon, out_folder_sp_dist + taxon + '_range_size.pkl')
        
    birds_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\IUCN_range_maps\\BIRDS\\'
    spat.create_array_sp_list_birds(birds_folder, out_folder_sp_dist + 'birds_' + str(pixel_size) + '.pkl', pixel_size = pixel_size)
    spat.create_sp_range_dic_bird(birds_folder, out_folder_sp_dist + 'birds_range_size.pkl')    

    # 2. Compute the weighted richness
    out_folder_weightedS = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\weighted_S\\'
    for taxon in taxa:
        array_taxon = spat.import_pickle_file(out_folder_sp_dist + taxon + '_' + str(pixel_size) + '.pkl')
        range_dic_taxon = spat.import_pickle_file(out_folder_sp_dist + taxon + '_range_size.pkl')
        for q in np.arange(0, 1.1, 0.1):
            spat.weighted_richness_to_raster(array_taxon, range_dic_taxon, -q, out_folder_weightedS, taxon)
     
    # Birds is a bit complicated b/c of migration
    
    # 3. Obtain rasters of environmental variables in same projection & resolution
    # Variables included: mean annual T, PET, AET, NDVI, altitudinal range, seasonality
    # mean annual T & seasonality are obtained from Bioclim (BIO1 and BIO4)
    in_folder_env = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\pred_vars_raw\\'
    out_folder_env = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\pred_vars\\'
    in_folder_bioclim = in_folder_env + 'bio_10m_bil\\'
    