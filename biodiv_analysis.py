from __future__ import division
import sys
sys.path.append('C:\\Users\\Xiao\\Documents\\GitHub\\biodiv_spat')
import spatial_functions as spat
import numpy as np
import shapely.ops, shapely.wkt

out_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\final_results\\'
continent_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\continents\\continent.shp'
taxa = ['amphibians', 'mammals']
continent_list = ['Africa', 'North America', 'South America', 'Eurasia']

# 1. Correlation b/w partial S & overall S for each of the four quartiles, globally and on each of the continents
continents_geom = spat.import_shapefile_field(continent_dir, field = 'CONTINENT')
asia, europe = continents_geom['Asia'], continents_geom['Europe']
eurasia = shapely.ops.cascaded_union([shapely.wkt.loads(asia), shapely.wkt.loads(europe)])
eurasia_wkt = shapely.wkt.dumps(eurasia)
continents_geom['Eurasia'] = eurasia_wkt

#for taxon in taxa: 
    #taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                #delimiter = ',', names = True)
    #taxon_range.dtype.names = ['sp', 'global'] + continent_list
    #taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    #out_dir = out_folder + 'partial_S\\' + taxon + '.csv'
    #spat.corr_richness_taxon_continent(taxon_range, taxon_list_array, 'global', out_dir)
    #for continent in continent_list:
        #cont_geom = continents_geom[continent]
        #spat.corr_richness_taxon_continent(taxon_range, taxon_list_array, continent, out_dir, cont_geom)

# 2. Correlation b/w partial S & overall S for each of the four quartiles, globally and on each of the continents
#        when species ranges are randomly distributed (either scattered or contiguous)
# 2a. Species ranges are scattered 
#for taxon in taxa: 
    #taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                #delimiter = ',', names = True)
    #taxon_range.dtype.names = ['sp', 'global'] + continent_list
    #taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    #out_dir = out_folder + 'partial_S\\' + taxon + '_sim_scattered.csv'
    #spat.corr_null_range(taxon_range, taxon_list_array, 'global', out_dir, Niter = 100)
    #for continent in continent_list:
        #cont_geom = continents_geom[continent]
        #spat.corr_null_range(taxon_range, taxon_list_array, continent, out_dir, cont_geom = cont_geom, Niter = 100)
## 2b. Species ranges are contiguous
#for taxon in taxa: 
    #taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                #delimiter = ',', names = True)
    #taxon_range.dtype.names = ['sp', 'global'] + continent_list
    #taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    #out_dir = out_folder + 'partial_S\\' + taxon + '_sim_contiguous.csv'
    #spat.corr_null_range(taxon_range, taxon_list_array, 'global', out_dir, sim_type = 'contiguous', Niter = 50)
    #for continent in continent_list:
        #cont_geom = continents_geom[continent]
        #spat.corr_null_range(taxon_range, taxon_list_array, continent, out_dir, cont_geom = cont_geom, sim_type = 'contiguous', Niter = 50)

# 3. Correlation betwen S(q) and S, both globally and on each continent
#for taxon in taxa: 
    #taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                #delimiter = ',', names = True)
    #taxon_range.dtype.names = ['sp', 'global'] + continent_list
    #taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    #out_dir = out_folder + 'sq\\' + taxon + '_corr.csv'
    #out_file = open(out_dir, 'ab')
    #for q in np.arange(-1, 1.1, 0.1):
        #corr_global = spat.corr_sq_s_continent(taxon_range, taxon_list_array, 'global', q)
        #print>>out_file, ','.join(map(str, ['global'] + [q, corr_global['pearson'], corr_global['spearman']]))
        #for continent in continent_list: 
            #cont_geom = continents_geom[continent]
            #corr_cont = spat.corr_sq_s_continent(taxon_range, taxon_list_array, continent, q, cont_geom)
            #print>>out_file, ','.join(map(str, [continent, "%.1f" % q, corr_cont['pearson'], corr_cont['spearman']]))
    #out_file.close()
    
# 4. Modling S(q) with environmental variables
AET_array = spat.import_raster_as_array(out_folder + 'env\\AET.tif')
spatT_array = spat.import_raster_as_array(out_folder + 'env\\spatT.tif')
tempT_array = spat.import_raster_as_array(out_folder + 'env\\tempT.tif')
for taxon in taxa: 
    out_dir = out_folder + 'sq\\' + taxon + '_sq_lm.csv'
    out_file = open(out_dir, 'ab')
    print>>out_file, ','.join(['continent', 'q', 'coef_AET', 'coef_spatT', 'coef_tempT', 'coef_int', \
                               'r2_full', 'r2_AET', 'r2_spatT', 'r2_tempT', 'r2_int'])
    taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                delimiter = ',', names = True)
    taxon_range.dtype.names = ['sp', 'global'] + continent_list
    taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    for q in np.arange(-1, 1.1, 0.1):
        reg_ans = spat.model_sq(taxon_range, taxon_list_array, 'global', AET_array, spatT_array, tempT_array, q)
        reg_ans = [x for x in reg_ans]
        print>>out_file, ','.join(map(str, ['global', "%.1f" % q] + reg_ans))
        for continent in continent_list:
            cont_geom = continents_geom[continent]
            reg_ans = spat.model_sq(taxon_range, taxon_list_array, continent, AET_array, spatT_array, tempT_array, q, cont_geom)
            reg_ans = [x for x in reg_ans]
            print>>out_file, ','.join(map(str, [continent, "%.1f" % q] + reg_ans))
    out_file.close()

