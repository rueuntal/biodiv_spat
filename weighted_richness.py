"""Script to examine how the importance of environmental variables shifts with weights on species range sizes."""
from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox')

import spatial_functions as spat
import psycopg2
from osgeo import gdal
import numpy as np
import subprocess
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# R script with function for multiple regression
STAP = SignatureTranslatedAnonymousPackage
with open('C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox\\regression.R', 'r') as f:
    string = f.read()
    
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
    
    # A separate map of birds with only resident species, & dictionary with species ranges for resident species
    spat.create_array_sp_list_birds(birds_folder, out_folder_sp_dist + 'birds_resident_' + str(pixel_size) + '.pkl',
                                    Attr = "SEASONAL", Attr_filter = "1", pixel_size = pixel_size)
    spat.create_sp_range_dic_bird(birds_folder, out_folder_sp_dist + 'birds_resident_range_size.pkl', Attr = "SEASONAL", 
                                  Attr_filter = "1")
    
    # 2. Compute the weighted richness
    out_folder_weightedS = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\weighted_S\\'
    taxa.extend(['birds_resident'])  # Now only looking at resident birds; TO BE MODIFIED LATER
    seabirds_list = np.genfromtxt('C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\seabirds.csv', dtype = None, \
                                   delimiter = ',', skip_header = 1, names = ['common name', 'binomial'])    
    for taxon in taxa:
        array_taxon = spat.import_pickle_file(out_folder_sp_dist + taxon + '_' + str(pixel_size) + '.pkl')
        range_dic_taxon = spat.import_pickle_file(out_folder_sp_dist + taxon + '_range_size.pkl')
        for q in np.arange(-10, 11, 1)/10: # This is to prevent weird approximation problem at zero
            if 'birds' in taxon:
                spat.weighted_richness_to_raster(array_taxon, range_dic_taxon, q, out_folder_weightedS, taxon, remove_sp_list = seabirds_list['binomial'])
            else:
                spat.weighted_richness_to_raster(array_taxon, range_dic_taxon, q, out_folder_weightedS, taxon)
     
    # Birds is a bit complicated b/c of migration
    
    # 3. Obtain rasters of environmental variables in same projection & resolution
    # Variables included: mean annual T, PET, AET, NDVI, altitudinal range, seasonality
    in_folder_env = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\pred_vars_raw\\'
    out_folder_env = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\pred_vars\\'
    match_file = out_folder_weightedS + taxa[0] + '_weighted_richness_' + str(pixel_size) + '_' + str(q) + '.tif'
    
    # mean annual T & seasonality are obtained from WordClim (BIO1 and BIO4)
    in_folder_bioclim = in_folder_env + 'bio_10m_bil\\'
    spat.reproj_raster_to_match(in_folder_bioclim + 'bio1.bil', out_folder_env + 'mean_annual_T.tif', match_file)
    spat.reproj_raster_to_match(in_folder_bioclim + 'bio4.bil', out_folder_env + 'seasonality.tif', match_file)
    
    # Altitudinal range is also obtained from WorldClim but requires additional processing in R
    in_folder_alt = in_folder_env + 'alt_30s_bil\\'
    spat.reproj_raster_pixel_size(in_folder_alt + 'alt.bil', out_folder_env + 'alt_behrmann_1000.tif', pixel_size= 1000) 
    # Get min and max in each pixel using R
    r_cmd = "C:\\Program Files\\R\\R-3.1.0\\bin\\x64\\Rscript.exe"
    r_script = "C:\\Users\\Xiao\\Documents\\GitHub\\gis_sandbox\\get_alt_min_max.R"
    subprocess.call([r_cmd, r_script])
    # Obtain altitudinal range from min and max
    alt_max_dir = out_folder_env + 'alt_behrmann_max.tif'
    alt_min_dir = out_folder_env + 'alt_behrmann_min.tif'
    alt_range_dir = out_folder_env + 'alt_range.tif'
    spat.get_range_raster(alt_max_dir, alt_min_dir, alt_range_dir)
    
    # PET is obtained from CGIAR-CSI
    in_folder_pet = in_folder_env + 'PET_he_annual\\pet_he_yr\\'
    spat.reproj_raster_to_match(in_folder_pet + 'w001001.adf', out_folder_env + 'PET.tif', match_file)
    
    # AET is obtained from MOD16
    file_AET = in_folder_env + 'mean_annual_AET\\MOD16A3_ET_2000_to_2013_mean.tif'
    # Replace the null values (>65500) with zeros
    AET_raster_orig = gdal.Open(file_AET)
    AET_band_orig  = AET_raster_orig.GetRasterBand(1)
    AET_array = AET_band_orig.ReadAsArray()
    AET_array[AET_array > 65500] = 0
    # Write modifed array back to file
    AET_geotrans = AET_raster_orig.GetGeoTransform()
    xmin, ymax, pixel_size = AET_geotrans[0], AET_geotrans[3], AET_geotrans[1]
    spat.convert_array_to_raster(AET_array, (xmin, ymax), out_folder_env + 'AET_raw.tif', pixel_size, 
                                 out_proj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    # Convert to Behrmann
    spat.reproj_raster_to_match(out_folder_env + 'AET_raw.tif', out_folder_env + 'AET.tif', match_file)
    
    # NDVI is obtained from FAO: http://www.fao.org/geonetwork/srv/en/metadata.show?id=37058
    file_NDVI = in_folder_env + 'Global annual sum NDVI Fig S5a\\Fig S5a.img'
    spat.reproj_raster_to_match(file_NDVI, out_folder_env + 'NDVI.tif', match_file)
    
    # 4. Simple regression of weighted richness versus predictors
    # Predictors are divided into three groups: temperature (mean annual T & PET), productivity (AET & NDVI),
    # and Janzen's hypothesis (seasonsality, altitudinal range, and interaction)
    # This section requires the module "rpy2" and calls the R script regression.R for multiple regression.
    out_analysis_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\output\\'
    multilin = STAP(string, 'multilin')    
    for taxon in taxa:
        for q in np.arange(-10, 11, 1)/10:
            S_dir = out_folder_weightedS + taxon + '_weighted_richness_' + str(pixel_size) + '_' + str(q) + '.tif'
            # The output is a list with 14 values: r-squared for the full and the three reduced models; 
            # unique contributions of the three reduced models; and p-value of the seven variables in the full model.
            out_taxon_q = [taxon, q] + [x for x in multilin.multilin(S_dir)]
            # Save output to a file 
            out_file = open(out_analysis_folder + 'multilin.txt', 'a')
            print>>out_file, '\t'.join(map(str, out_taxon_q))
            out_file.close()
    
    # Plot q versus r-squared for the four models, as well as unique contributions from each set of variables, for each taxon
    out_plot_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\plots\\'
    results_multilin_dir = out_analysis_folder + 'multilin.txt'
    for taxon in taxa:
        out_plot_dir = out_plot_folder + taxon + '_r2_multilin.png'
        spat.plot_r2_multilin(results_multilin_dir, taxon, out_plot_dir)
        