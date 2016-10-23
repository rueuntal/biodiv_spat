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

for taxon in taxa: 
    taxon_range = np.genfromtxt(out_folder + 'sp_dist\\' + taxon + '_range_size.csv', dtype = None, \
                                delimiter = ',', names = True)
    taxon_range.dtype.names = ['sp', 'global'] + continent_list
    taxon_list_array = spat.import_pickle_file(out_folder + 'sp_dist\\' + taxon + '_dist.pkl')
    out_dir = out_folder + 'partial_S\\' + taxon + '.csv'
    spat.corr_richness_taxon_continent(taxon_range, taxon_list_array, 'global', out_dir)
    for continent in continent_list:
        cont_geom = continents_geom[continent]
        spat.corr_richness_taxon_continent(taxon_range, taxon_list_array, continent, out_dir, cont_geom)
    