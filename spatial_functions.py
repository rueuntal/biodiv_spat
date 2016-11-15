from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
import cPickle
from osgeo import ogr, osr, gdal
import shapely.wkt, shapely.ops
import multiprocessing
from matplotlib.mlab import PCA
import re
from math import sin, cos, sqrt, atan2
import csv
from matplotlib.pyplot import cm
from scipy.stats.stats import pearsonr, spearmanr
from collections import Counter
import random
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
import rpy2.robjects as rob
STAP = SignatureTranslatedAnonymousPackage
with open('C:\\Users\\Xiao\\Documents\\GitHub\\biodiv_spat\\rfunctions.R', 'r') as f:
    string = f.read()
rfunctions = STAP(string, 'rfunctions')     

def import_pickle_file(in_dir):
    """Read in pickled file."""
    in_file = open(in_dir, 'rb')
    in_obj = cPickle.load(in_file)
    in_file.close()
    return in_obj

def import_shapefile(in_dir, Attr = None, AttrFilter = None): 
    """Read in a .shp file and save the geometry of each feature in a list of WKT."""
    driver = ogr.GetDriverByName('ESRI Shapefile')
    datasource = driver.Open(in_dir, 0)
    layer = datasource.GetLayer()
    if Attr:
        filter_string = Attr + "='" + AttrFilter + "'"
        layer.SetAttributeFilter(filter_string)
    geom_list = []
    for feature in layer:
        geom_list.append(feature.GetGeometryRef().ExportToWkt())
    return geom_list

def import_shapefile_field(in_dir, Attr = None, AttrFilter = None, field = 'binomial'):
    """Read in a .shp file and save the geometry(ies) for unique values of field as a dictionary."""
    driver = ogr.GetDriverByName('ESRI Shapefile')
    datasource = driver.Open(in_dir, 0)
    layer = datasource.GetLayer()
    if Attr:
        filter_string = Attr + " IN " + AttrFilter
        layer.SetAttributeFilter(filter_string)
    geom_dic = {}
    for feature in layer:
        val = feature.GetField(field)
        geom_wkt = feature.GetGeometryRef().ExportToWkt()
        if val not in geom_dic: geom_dic[val] = [geom_wkt]
        else: geom_dic[val].append(geom_wkt)
    for val in geom_dic.keys():
        geom_shapes =  [shapely.wkt.loads(x) for x in geom_dic[val]]
        try:
            geom_wkt_union = shapely.wkt.dumps(shapely.ops.cascaded_union(geom_shapes))
        except: 
            geom_shapes = [x.buffer(0) for x in geom_shapes] # Adding zero buffer somehow helps
            geom_wkt_union = shapely.wkt.dumps(shapely.ops.cascaded_union(geom_shapes))
        geom_dic[val] = geom_wkt_union
    return geom_dic

def import_shapefile_folder(in_folder, Attr = None, AttrFilter = None):
    """Similar to import_shapefile_field, but for the case where sperate .shp files for different species 
    
    are in one folder (e.g., BIEN). 
    
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    file_list = os.listdir(in_folder)
    sp_list = sorted(list(set([x.split('.')[0] for x in file_list])))
    geom_dic = {}
    for sp in sp_list:
        sp_dir = in_folder + sp + '.shp'
        try:
            geom_wkt = import_shapefile(sp_dir, Attr, AttrFilter)[0]
            geom_dic[sp] = geom_wkt      
        except: 
            pass
    return geom_dic
    
def import_raster_as_array(raster_dir, nodata = None):
    raster = gdal.Open(raster_dir)
    band = raster.GetRasterBand(1).ReadAsArray()
    if nodata is not None:
        nodata_orig = raster.GetRasterBand(1).GetNoDataValue()
        band[band == nodata_orig] = nodata
    return band

def proj_name_to_proj4(proj_name):
    """Given name of a projection, return its specification in PROJ4."""
    name_and_proj = {'latlong': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                     'behrmann': '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
                     }
    return name_and_proj[proj_name]

def proj_extent(proj_name):
    """Give the global extent (xmin, xmax, ymin, ymax) of a projection for raster. Only Behrmann is available now."""
    in_proj = 'latlong'
    transform = reproj(in_proj, proj_name)
    x, y, z = transform.TransformPoint(180, 90)
    return [-x, x, -y, y]

def prob_of_presence(focal_sp, dic_of_sp_range, S, Niter = 10000):
    """Using weighted random sampling without replacement to compute the probability of 
    
    presence of a given species in a grid cell with local richness S. 
    Note that no analytical solution exists for this problem. 
    Inputs: 
    focal_sp: species of interest, has to be in dic_of_sp_range
    dic_of_sp_range: a dictionary of range sizes for species of interest. Range sizes are in unit of number of grid cells.
    Niter: number of repeated samples to compute probability.
    
    Output: 
    Probability of presence, a single float number.
    
    """
    sp_list, range_list = [], []
    for sp in dic_of_sp_range:
        sp_list.append(sp)
        range_list.append(dic_of_sp_range[sp])
    p = np.array(range_list) / np.sum(range_list)
    counter = 0
    for i in range(Niter):
        sample = np.random.choice(sp_list, size = S, replace = False, p = p)
        if focal_sp in sample: counter += 1
    return counter / Niter

def write_raster_to_file(out_dir, wide, high, geotrans, proj_wkt, nodata = 0, dtype = gdal.GDT_Float32):    
    """Create an empty raster at a specified path (in memory if path is None), define geotransform, projection, nodata etc. 
    
    and return the object (filled with no data).
    
    """
    if out_dir: outRaster = gdal.GetDriverByName('GTiff').Create(out_dir, wide, high, 1, dtype)
    else: outRaster = gdal.GetDriverByName('MEM').Create('', wide, high, 1, dtype)
    if geotrans: outRaster.SetGeoTransform(geotrans)
    if proj_wkt: outRaster.SetProjection(proj_wkt)
    outRaster.GetRasterBand(1).SetNoDataValue(nodata)
    outRaster.GetRasterBand(1).Fill(nodata)
    return outRaster    
 
def get_distance_latlon(lat1, lon1, lat2, lon2):
    """Compute the distance (in meters) between two points given their latitudes and longitudes (in degrees)."""
    R = 6373 * 10**3
    lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c
    return distance
    
def reproj(in_proj = 'latlong', out_proj = 'behrmann'):
    """General function for reprojection which can be used for both vectors and rasters."""
    in_proj4 = proj_name_to_proj4(in_proj)
    out_proj4 = proj_name_to_proj4(out_proj)
    in_proj = osr.SpatialReference()
    in_proj.ImportFromProj4(in_proj4)
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(out_proj4)
    transform = osr.CoordinateTransformation(in_proj, out_proj)
    return transform

def reproj_geom(in_geom, in_proj = 'latlong', out_proj = 'behrmann'):
    """Function to reproject geometries (defined by WKT). 
    
    Default input is unprojected WGS84 (EPSG 4326), default output is Behrmann equal area.
    
    """
    in_geom_ogr = ogr.CreateGeometryFromWkt(in_geom)
    transform = reproj(in_proj, out_proj)
    in_geom_ogr.Transform(transform)
    return in_geom_ogr.ExportToWkt()

def weighted_sample_range_size(range_size_list, sample_size):
    """Generate a random sample of a certain size (length) from a list of range sizes,
    
    weighted by the range sizes.
    
    """
    np.random.seed()
    sum_size = np.sum(range_size_list)
    selected_range = np.random.choice(range_size_list, size = sample_size, replace = False, \
                                      p = np.array(range_size_list) / sum_size)
    return selected_range

def create_array_for_raster(extent, geom = None, no_value = 0, pixel_size = 100000):
    """Createa an array that can be manipulated and/or later converted to raster.
    
    Input:
    extend: a list in the form of [xmin, xmax, ymin, ymax]. Can be self-defined, or take world extend from proj_extend() for a projection.
    geom: Default is empty, leading to an array filled with non-values. Or could be a multipolygon (WKT) in the desired projection,
        leading to 1's on the landscape.
    no_value: value to fill the grid with no data. Default is 0.
    pixel_size: size of each pixel. Default is 100km. 
    
    """
    xmin, xmax, ymin, ymax = extent
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    
    geotrans = (xmin, pixel_size, 0, ymax, 0, -pixel_size)
    target_ds = write_raster_to_file(None, x_res, y_res, geotrans, None, nodata = no_value, dtype = gdal.GDT_Int32)
    if geom: 
        # Convert geom to layer
        virtual_shp = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('/vsimem/virtual.shp')
        virtual_layer = virtual_shp.CreateLayer('sp_range', geom_type = ogr.wkbPolygon)
        feature = ogr.Feature(virtual_layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkt(geom))
        virtual_layer.CreateFeature(feature)
        # This currently throws a warning about missing spatial ref
        gdal.RasterizeLayer(target_ds, [1], virtual_layer, None, None, [1],  ['ALL_TOUCHED=TRUE'])
    band = target_ds.GetRasterBand(1)
    return band.ReadAsArray()

def convert_array_to_raster(array, rasterOrigin, out_file, pixel_size, no_value = 0, out_proj = 'behrmann'):
    """Convert an array to raster."""
    cols = array.shape[1]
    rows = array.shape[0]
    xmin, ymax = rasterOrigin
    geotrans = (xmin, pixel_size, 0, ymax, 0, -pixel_size)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(proj_name_to_proj4(out_proj))
    proj_wkt = outRasterSRS.ExportToWkt()
    
    outRaster = write_raster_to_file(out_file, cols, rows, geotrans, proj_wkt, nodata = no_value)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRaster.FlushCache()
    outRaster = None
    return None
    
def richness_to_raster(sp_array, out_dir, proj = 'behrmann', pixel_size = 100000, remove_sp_list = []):  
    """Convert diversity to raster"""
    xmin, xmax, ymin, ymax = proj_extent(proj)
    table_richness = np.zeros(shape = sp_array.shape)
    for j in range(len(sp_array)):
        for i in range(len(sp_array[0])):
            sp_grid = sp_array[j][i]
            richness = len([sp for sp in sp_grid if sp not in remove_sp_list])
            table_richness[j][i] = richness
    convert_array_to_raster(table_richness, [xmin, ymax], out_dir, pixel_size, out_proj = proj)
    return None

#def richness_to_raster_by_family(postgis_cur, table_name, out_dir, pixel_size = 100000):
    #"""Convert an IUCN shapefile with range maps of a taxon into rasters of richness, one for each family, 
    
    #assuming that the shapefile has already been loaded into a PostGIS database.
    #Inputs:
    #postgis_cur - cursor pointing to the PostGIS database
    #table_name - name of the table in the PostGIS databse
    #out_dir - output directory for the resulting raster
    #pixel_size - grain size of raster, default is 100km * 100km
    
    #"""
    #xmin, xmax, ymin, ymax = proj_extent('behrmann')
    #family_list_exec = 'SELECT DISTINCT family_nam FROM ' + table_name
    #postgis_cur.execute(family_list_exec)
    #family_list = [x[0] for x in postgis_cur.fetchall()]
    
    #for family in family_list:
        #out_file = out_dir + '/' + table_name + '_' + family + '_' + str(int(pixel_size)) + '.tif'
        ## Check to make sure file doesn't exist already - useful to resume from unexpected breaks
        #if not os.path.exists(out_file):
            #sp_list_exec = "SELECT DISTINCT binomial FROM " + table_name + " WHERE family_nam='" + family + "'"
            #postgis_cur.execute(sp_list_exec)
            #sp_list = [x[0] for x in postgis_cur.fetchall()]
            #family_landscape = create_array_for_raster(proj_extent('behrmann'), pixel_size = pixel_size)
            ## Create empty array 
            #for sp in sp_list:
                #wkt_reproj = sp_reproj(postgis_cur, table_name, sp)
                ## Convert species range to raster array
                #sp_landscape = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)        
                #family_landscape += sp_landscape
            #convert_array_to_raster(family_landscape, [xmin, ymax], out_file, pixel_size)
    #return None

def create_sp_range_dic(geom_dic_sp, continent = None, proj = 'behrmann'):
    """Obtain the range size of species (in m^2) either globally (when "continent" = None) or 
    
    on one continent ("continent" has to be wkt of the polygon in WGS 84). The output is a dictionary. 
    
    """
    sp_range_dic = {}
    if continent is not None: reproj_continent = reproj_geom(continent, out_proj = proj)
    for sp in geom_dic_sp.keys():
        wkt_sp = geom_dic_sp[sp]
        reproj_sp = reproj_geom(wkt_sp, out_proj = proj)
        sp_range = shapely.wkt.loads(reproj_sp)
        if continent is not None:
            continent_range = shapely.wkt.loads(reproj_continent)
            try:
                sp_range_continent = sp_range.intersection(continent_range)
            except:
                continent_range = continent_range.buffer(0)
                sp_range = sp_range.buffer(0)
                sp_range_continent = sp_range.intersection(continent_range)
        else: sp_range_continent = sp_range
        sp_range_dic[sp]  = sp_range_continent.area
    return sp_range_dic

def create_array_sp_list(geom_dic_sp, out_file_name, pixel_size = 100000, out_proj = 'behrmann'):
    """Create an array with species list in each grid covering the globe."""
    xmin, xmax, ymin, ymax = proj_extent(out_proj)
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    array_list = np.array([[[] for i in range(x_res)] for j in range(y_res)])
    
    for sp in geom_dic_sp.keys():
        sp_geom_wkt = geom_dic_sp[sp]             
        wkt_reproj = reproj_geom(sp_geom_wkt, out_proj = out_proj)        
        # Convert species range to raster array
        sp_array = create_array_for_raster(proj_extent(out_proj), geom = wkt_reproj, pixel_size = pixel_size)        
        array_list = np.array([[list(array_list[j][i]) + [sp] if sp_array[j][i] > 0 else list(array_list[j][i]) for i in range(x_res)] for j in range(y_res)])
    
    # Save to file
    out_file = open(out_file_name, 'wb')
    cPickle.dump(array_list, out_file, protocol = 2)
    out_file.close()
    return None

def metric_dist(dist, metric):
    """Sub-function of compare_range_size_dists calcuting the designated metric for a distribution. 
    
    Inputs:
    dist: a list of range sizes, either emprirical or randomly generated.
    metric: metric  to compute. See docstring for compare_range_size_dists() for details. 
    
    """
    log_dist = np.log(dist)
    if metric == 'mean': return np.mean(log_dist)
    elif metric == 'sd': return np.std(log_dist, ddof = 1)
    elif metric == 'skew': return stats.skew(log_dist)

def obtain_metrics_single(input_list):
    """Subfunction for multiprocessing, called in compare_range_size_dists."""
    range_size_list, sample_size, metrics = input_list
    rand_dist_grid = weighted_sample_range_size(range_size_list, sample_size)
    metrics_single =  [metric_dist(rand_dist_grid, metric) for metric in metrics]
    return metrics_single

def range_size_dists_raw(sp_list_array, range_size_dic, out_dir, out_name, metrics = ['mean', 'sd', 'skew'], marine_list = [], threshold = 5):    
    """Obtain the mean, sd, and skewness of the range size distribution in each pixel and save to file."""
    array_out_list = [np.empty([len(sp_list_array), len(sp_list_array[0])], dtype = float) for k in range(len(metrics))]
        
    for j in range(len(sp_list_array)):
        for i in range(len(sp_list_array[0])):
            sp_grid = sp_list_array[j][i]
            # Remove marine species from local species list
            sp_grid = [sp for sp in sp_grid if sp not in marine_list]
            richness_grid = len(sp_grid)
            if richness_grid < threshold: 
                for k in range(len(metrics)):
                    array_out_list[k][j][i] = -1 # Fill unanalyzed grids with -1
            else:
                emp_range_dist_grid = [range_size_dic[sp] for sp in sp_grid]        
                for k in range(len(metrics)):
                    array_out_list[k][j][i] = metric_dist(emp_range_dist_grid, metrics[k])
    
    # Save to file
    for k, metric in enumerate(metrics):
        out_file_name = out_dir + '/' + out_name + '_' + metric + '.pkl'
        out_file = open(out_file_name, 'wb')
        cPickle.dump(array_out_list[k], out_file, protocol = 2)
        out_file.close()
    return None
    
def compare_range_size_dists(sp_list_array, range_size_dic, out_dir, out_name, Nsample, metrics = ['mean', 'sd', 'skew'], marine_list = [], threshold = 5):
    """Compare with empirical range size distribution in each grid with the expected distribution, 
    
    and save the output into an array.
    
    Inputs:
    sp_list_array: array with species list in each grid, created by create_array_sp_list().
    range_size_dic: dictionary of range size for each species, created by create_sp_range_dic().
    out_dir: output directory
    out_name: name of the output file, e.g., "terrestrial_mammals", "reptiles", etc. 
    Nsample: number of random samples to draw from the full distribution.
    metric: list of metrics used to compare the empirical and randomly generated size distributions, which can contain the following values:
    "mean", "sd" (standard deviation), "skew" (skewness). All values are calcuated on log scale.
    marine_list: if there are any marine species to be removed from the analysis
    threshold: minimal species richness of a grid, below which it is not analyzed.
    
    Output:
    File with an array of the same dimension as sp_list_array, where the value of each grid is the quantile of the empirical size distribution among the 
    randomly generated size distributions, measured with the designated metric. 
    The output file is written to the designated output directory, with the same out_name + "_" + metric + ".pkl".
    
    """
    array_out_list = [np.empty([len(sp_list_array), len(sp_list_array[0])], dtype = float) for k in range(len(metrics))]
    # Remove marine species from the dictionary of range sizes  
    for marine_sp in marine_list: del range_size_dic[marine_sp]
    range_size_list = range_size_dic.values()  
        
    for j in range(len(sp_list_array)):
        for i in range(len(sp_list_array[0])):
            sp_grid = sp_list_array[j][i]
            # Remove marine species from local species list
            sp_grid = [sp for sp in sp_grid if sp not in marine_list]
            richness_grid = len(sp_grid)
            if richness_grid < threshold: 
                for k in range(len(metrics)):
                    array_out_list[k][j][i] = -1 # Fill unanalyzed grids with -1
            else:
                emp_range_dist_grid = [range_size_dic[sp] for sp in sp_grid]
                emp_metrics = [metric_dist(emp_range_dist_grid, metric) for metric in metrics]
                
                pool = multiprocessing.Pool()
                rand_metrics_list = pool.map(obtain_metrics_single, [[range_size_list, richness_grid, metrics] for m in range(Nsample)])
                pool.close()
                pool.join()
                
                # Compute quantile
                for k in range(len(metrics)):
                    rand_metric_list = [rand_single_iter[k] for rand_single_iter in rand_metrics_list]
                    quan = len([x for x in rand_metric_list if x <= emp_metrics[k]]) / Nsample
                    array_out_list[k][j][i] = quan
    
    # Save to file
    for k, metric in enumerate(metrics):
        out_file_name = out_dir + '/' + out_name + '_' + metric + '.pkl'
        out_file = open(out_file_name, 'wb')
        cPickle.dump(array_out_list[k], out_file, protocol = 2)
        out_file.close()
    return None

def reproj_raster_generic(in_file, out_dir, wide, high, geotrans, in_proj, out_proj, nodata, alg = gdal.GRA_Bilinear):
    """Generic function to reproject raster."""
    outRaster = write_raster_to_file(out_dir, wide, high, geotrans, out_proj, nodata = nodata)
    res = gdal.ReprojectImage(in_file, outRaster, in_proj, out_proj, alg)
    outRaster.FlushCache()
    return outRaster

def reproj_raster_pixel_size(in_dir, out_dir, pixel_size = 100000, alg = gdal.GRA_Bilinear, \
                             in_proj_name = 'latlong', out_proj_name = 'behrmann'):
    """Reproject a raster with projection and pixel size spcified."""
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(proj_name_to_proj4(out_proj_name))
    out_proj_wkt = out_proj.ExportToWkt()
    in_file = gdal.Open(in_dir)
    try: in_proj_wkt = in_file.GetProjection()
    except: 
        in_proj = osr.SpatialReference()
        in_proj.ImportFromProj4(proj_name_to_proj4(in_proj_name))
        in_proj_wkt = in_proj.ExportToWkt()
    in_nodata = in_file.GetRasterBand(1).GetNoDataValue()
    out_extent = proj_extent(out_proj_name)
    xmin, xmax, ymin, ymax = out_extent
    wide = int((xmax - xmin) / pixel_size)
    high = int((ymax - ymin) / pixel_size)
    out_geotrans = (xmin, pixel_size, 0, ymax, 0, -pixel_size)
    reproj_raster_generic(in_file, out_dir, wide, high, out_geotrans, in_proj_wkt, out_proj_wkt, in_nodata)
    return None
    
def reproj_raster_to_match(in_dir, out_dir, match_dir, alg = gdal.GRA_Bilinear, no_data = None, \
                  in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', \
                  out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Similar to reproj_raster(), but the reprojected/resampled raster matches the grids of another raster."""
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(out_proj4)    
    in_file = gdal.Open(in_dir)
    try: in_proj_wkt = in_file.GetProjection()
    except: 
        in_proj = osr.SpatialReference()
        in_proj.ImportFromProj4(in_proj4)
        in_proj_wkt = in_proj.ExportToWkt()
    in_nodata = in_file.GetRasterBand(1).GetNoDataValue()
    if no_data is not None: # Self-define nodata value
        in_wide = in_file.RasterXSize
        in_high = in_file.RasterYSize
        in_geotrans = in_file.GetGeoTransform()
        in_array = in_file.GetRasterBand(1).ReadAsArray()
        in_array[in_array == in_nodata] = no_data
        new_file = write_raster_to_file(None, in_wide, in_high, in_geotrans, in_proj_wkt, nodata = no_data)
        new_file.GetRasterBand(1).WriteArray(in_array)
        in_file = new_file
        in_nodata = no_data
  
    match_file = gdal.Open(match_dir)
    match_proj_wkt = match_file.GetProjection()
    match_geotrans = match_file.GetGeoTransform()
    wide = match_file.RasterXSize
    high = match_file.RasterYSize
    
    out_raster = reproj_raster_generic(in_file, out_dir, wide, high, match_geotrans, in_proj_wkt, match_proj_wkt, in_nodata)
    return out_raster

def get_range_raster(in_dir_high, in_dir_low, out_dir):
    """Obtain the difference between two rasters and write to a new raster. 
    
    All files are in the same projection and resolution.
    
    """
    high_file = gdal.Open(in_dir_high)
    high_array = high_file.GetRasterBand(1).ReadAsArray()
    low_file = gdal.Open(in_dir_low)
    low_array = low_file.GetRasterBand(1).ReadAsArray()
    out_array = high_array - low_array
    
    proj_wkt = high_file.GetProjection()   
    geotrans = high_file.GetGeoTransform()
    wide = high_file.RasterXSize
    high = high_file.RasterYSize
    nodata_high = high_file.GetRasterBand(1).GetNoDataValue()
    nodata_low = low_file.GetRasterBand(1).GetNoDataValue()
    for i in range(len(out_array)): # Fill in nodata values
        for j in range(len(out_array[0])):
            if high_array[i][j] == nodata_high or low_array[i][j] == nodata_low:
                out_array[i][j] = nodata_high

    outRaster = write_raster_to_file(out_dir, wide, high, geotrans, proj_wkt, nodata = nodata_high)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(out_array)
    outRaster.FlushCache()
    outRaster = None
    return None
    
def raster_reproj_flat(in_dir, match_dir, log = False):
    """Subfunction. Reproject a raster to match another file, then flatten the array to 1-D to be used in PCA."""
    in_file = reproj_raster_to_match(in_dir, None, match_dir)
    in_array = in_file.GetRasterBand(1).ReadAsArray()
    in_nodata = in_file.GetRasterBand(1).GetNoDataValue()
    in_flat = np.ravel(in_array)
    in_flat = in_flat.astype('float')
    in_flat[in_flat == in_nodata] = np.nan
    if log:
        if min(in_flat) < 0: 
            print "Error: cannot log-transform negative values."
            return in_flat
        else:
            in_flat[in_flat == 0] = 1 # replace zeros with 1's
            in_flat = np.log(in_flat)       
    return in_flat
    
def PCA_with_NA(array, num_axes):
    """Perform PCA with an array (remove rows with nan), and keep the specified number of axes."""
    array_no_nan = array[~np.isnan(array).any(axis = 1)]
    array_PCA = PCA(array_no_nan)
    array_list = []
    for i in range(array.shape[0]):
        if np.isnan(array[i]).any(): array_list.append([np.nan] * num_axes)
        else: 
            row_proj = array_PCA.project(array[i])
            array_list.append(list(row_proj[:num_axes]))       
    return array_list

def convert_point_latlon(x, y, proj = 'behrmann'):
    """Convert a point to lat/lon."""
    point = ogr.CreateGeometryFromWkt("POINT (" + str(x) + " " + str(y) + ")")
    in_proj4 = proj_name_to_proj4(proj)
    out_proj4 = proj_name_to_proj4('latlong')
    point_reproj = reproj_geom(point.ExportToWkt(), out_proj4 = out_proj4, in_proj4 = in_proj4)    
    lon, lat = re.findall(r'[+-]?\d+.\d+', point_reproj)
    return float(lon), float(lat)

def pixel_center_to_latlon(i, j, pixel_size = 100000, proj = 'behrmann'):
    """Obtain the center of the pixel (ith row jth column from top left corner) and convert it to lat/lon."""
    extent = proj_extent(proj)
    xmin =  extent[0]
    ymax = extent[3] # The top left corner is xmin-ymax
    xcoord = xmin + (j + 0.5) * pixel_size
    ycoord = ymax - (i + 0.5) * pixel_size
    lon, lat = convert_point_latlon(xcoord, ycoord, proj = proj)
    return lon, lat

def comp_pixel_to_pixel(pca_out_flat, i, j, m, n, ncol, radius, pixel_size):
    """Subfunction to obtain the Euclidean distance in environmental PCA between two pixels.
    
    inputs:
    pca_out_flat - pca output summarizing the environmental variables in each pixel.
    i, j, m, n - x and y coordinates of the two pixels in the raster.
    ncol - number of columns in raster
    radius - if distance between two pixels is larger than radius, return nan.
    
    """
    if np.isnan(pca_out_flat[m * ncol + n]).any(): return np.nan
    elif i == m and j == n: return np.nan
    else:
        if n >= ncol: n -= ncol # wrap to the other side if column is out of bound
        elif n < 0: n += ncol
        lon1, lat1 =  pixel_center_to_latlon(i, j, pixel_size = pixel_size)
        lon2, lat2 = pixel_center_to_latlon(m, n, pixel_size = pixel_size)
        dist = get_distance_latlon(lat1, lon1, lat2, lon2)
        if dist > radius: return np.nan
        else: return np.linalg.norm(np.array(pca_out_flat[i*ncol + j]) - np.array(pca_out_flat[m*ncol + n]))
         
def get_unique_raster(bio_dir, num_axes, radius, match_dir, out_dir, \
                      proj = 'behrmann', buffer = 1.3):
    """Obtain the measure of uniqueness of each cell as compared to cells nearby, sensu Morueta et al. 
    
    Inputs:
    bio_dir: folder with bioclim layers
    num_axes: number of axes to keep in PCA
    radius: radius of the circle centered on each cell, within which cells are compared with the focus cell
    match_file: raster file as a template for calculations and the output
    out_dir: directory of the output file
    buffer: 
    
    """
    # Get PCA
    full_bio = []
    log_trans_list = [3, 12, 13, 14, 16, 17, 18, 19] # Layers to be log-transformed 
    for i in range(1, 20):
        bio_i_dir = bio_dir + '\\bio' + str(i) + '.bil'
        if i in log_trans_list: 
            bio_i_flat = raster_reproj_flat(bio_i_dir, match_dir, log = True)
        else: bio_i_flat = raster_reproj_flat(bio_i_dir, match_dir)
        full_bio.append(bio_i_flat)
     
    full_bio_array = np.asarray(full_bio)
    full_bio_array = full_bio_array.T
    pca_out_flat = PCA_with_NA(full_bio_array, num_axes)

    match_file = gdal.Open(match_dir)
    match_geotrans = match_file.GetGeoTransform()
    xmin = match_geotrans[0]
    ymax = match_geotrans[3]
    pixel_size = match_geotrans[1]
    no_value = match_file.GetRasterBand(1).GetNoDataValue()
    
    # Comparison
    out_array = create_array_for_raster([xmin, -xmin, -ymax, ymax], no_value = no_value, pixel_size = pixel_size)
    out_array = out_array.astype(float)
    nrow, ncol = out_array.shape
    num_cells = int(np.ceil(radius * buffer / pixel_size))
    for i in range(nrow):
        for j in range(ncol):
            if not np.isnan(pca_out_flat[i*ncol + j]).any(): 
                ij_neighbour_distance = np.array([comp_pixel_to_pixel(pca_out_flat, i, j, m, n, ncol, radius, pixel_size) \
                                         for m in range(max(0, i - num_cells), min(nrow, i + num_cells + 1)) \
                                         for n in range(j - num_cells, j + num_cells + 1)])
                ij_neighbour_distance = ij_neighbour_distance[~np.isnan(ij_neighbour_distance)]
                out_array[i][j] = np.mean(ij_neighbour_distance)
    
    # Save to raster
    convert_array_to_raster(out_array, [xmin, ymax], out_dir, pixel_size, no_value = no_value)
    return None

def weighted_richness_to_raster(array_sp_list, range_dic, q, out_dir, out_name, pixel_size = 100000, remove_sp_list = []):
    """Obtains the raw or weighted richness, given an array with species list in each cell, 
    
    a dictionary with species range sizes, and weight (q). 
    
    The weighing function takes the form S(q) = S_tot**q * sum((range_i / range_tot) ** q)
    So S(q = 0) is the raw richness, while S(q = -1) is where each species is inversely weighted by their range, 
    i.e., all species get the same total count across the landscape.
    The scaling factor, S_tot**q, ensures that when all species are present in a cell, and when they all have the 
    same range size, S(q) = S_tot regardless of q. 
    
    The output is saved to disc as a raster in Behrmann projection.
    """
    weighted_S_array = np.zeros(shape = array_sp_list.shape)
    S_list = []
    for col in array_sp_list:
        for row in col:
            S_list.extend(row)
    S_list_unique = [sp for sp in set(S_list) if sp not in remove_sp_list]
    S_tot = len(S_list_unique)
    range_tot = np.sum([range_dic[x] for x in S_list_unique])
    for i, col in enumerate(array_sp_list):
        for j, row in enumerate(col):
            sp_row = [sp for sp in row if sp not in remove_sp_list]
            weighted_S_array[i][j] = S_tot ** (-q) * np.sum([(range_tot/range_dic[x]) ** q  for x in sp_row])
    
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    out_file = out_dir + '/' + out_name + '_weighted_richness_' + str(pixel_size) + '_' + str(q) + '.tif'
    convert_array_to_raster(weighted_S_array, [xmin, ymax], out_file, pixel_size)
    return None

def plot_r2_weighted_S(weighted_S_folder, taxon, file_ext, ax = None):
    """Obtain r2 between S and weighted S at different values of q, 
    
    and plot them against q. 
    Inputs:
    weighted_S_folder: folder where the raster files of weighted S are located
    taxon: taxon of interest
    file_ext: "middle part" of file name, current in the form of "_weighted_richness_100000_"
        The full file name would be taxon + file_ext + q + '.tif'
    ax: plotting frame. If None, a new one will be created.
    Returns ax with the plot.
    
    """
    r2_list = []
    orig_S_dir = weighted_S_folder + taxon + file_ext + str(0.0) + '.tif'
    orig_S = np.ndarray.flatten(import_raster_as_array(orig_S_dir))
    orig_S_no_zero = orig_S[orig_S > 0]
    for q in np.arange(-10, 11, 1)/10:
        weighted_S_dir = weighted_S_folder + taxon + file_ext + str(q) + '.tif'
        weighted_S = np.ndarray.flatten(import_raster_as_array(weighted_S_dir))
        weighted_S_no_zero = weighted_S[orig_S > 0]
        r = stats.pearsonr(orig_S_no_zero, weighted_S_no_zero)[0] 
        r2_list.append(r)
    # Plot
    plt.plot(np.arange(-10, 11, 1)/10, r2_list, c = '#000000', linewidth = 2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('Correlation r between weighted S and raw S', fontsize = 8)
    return ax 
        
def plot_r2_multilin(results_dir, taxon, out_dir, legend = True):
    """Plot the weighing parameter q against R^2 values (full model, two reduced models, 
    
    and the unique contributions.
    The output is a 1*3 plot, with R^2 of the three models in the left subplot, unique contributions 
    of the two reduced models in the middle subplot, and unique contributions of the three parameters
    in Janzen's hypothesis in the right subplot.
    Inputs:
    results_dir: txt file with the output from multiple regression, generated in weighted_richness.py
    taxon: string, name of the taxon
    out_dir: directory for the output figure
    legend: whether legend is to be included in the figure
    
    """
    results = np.genfromtxt(results_dir, dtype = None, delimiter = '\t')
    results_taxon = results[results['f0'] == taxon]
    q, tot_r2, prod_r2, janzen_r2, prod_r2_part, janzen_r2_part, alt_r2_part, seas_r2_part, int_r2_part = \
        zip(*results_taxon[['f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9']])
    r2 = [tot_r2, prod_r2, janzen_r2]
    unique_r2 = [prod_r2_part, janzen_r2_part]
    unique_r2_Janzen = [alt_r2_part, seas_r2_part, int_r2_part]
    col = ['#000000','#1874CD', '#8B7765']
    col_Janzen = ['#CD853F', '#EE799F', '#912CEE']
    
    fig = plt.figure(figsize = (11, 4))
    ax1 = plt.subplot(1, 3, 1) 
    line_list = [None] * 3
    for i, r2_val in enumerate(r2):
        line_list[i],  = plt.plot(q, r2_val, c = col[i], linewidth = 2)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)
    if legend:
        plt.legend(line_list, ['Full', 'Productivity', 'Janzen'], loc = 2, prop = {'size': 10})    
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('R-squared', fontsize = 8)
    plt.ylim((0, 1))
    
    ax2 = plt.subplot(1, 3, 2)
    for j, unique_r2_val in enumerate(unique_r2):
        plt.plot(q, unique_r2_val, c = col[j + 1], linewidth = 2)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('Unique R-squared', fontsize = 8)
    plt.ylim((0, 0.7))
    
    ax3 = plt.subplot(1, 3, 3)
    line_list_Janzen = [None] * 3
    for k, unique_r2_val in enumerate(unique_r2_Janzen):
        line_list_Janzen[k], = plt.plot(q, unique_r2_val, c = col_Janzen[k], linewidth = 2)
    ax3.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('Unique R-squared', fontsize = 8)
    plt.ylim((0, 0.7))
    if legend:
        plt.legend(line_list_Janzen, ['Altitude', 'Seasonality', 'Interaction'], loc = 2, prop = {'size': 10})    
    
    plt.subplots_adjust(wspace = 0.29)
    plt.savefig(out_dir, dpi = 600)

def plot_quantile_constraint(out_dir, q_list, r2_list, pred_vars, pred_var_prop_list):
    """This function creates a 2-panel plot for the output from quantile regression on constraining factors.
    
    The left panel shows how r^2 between predicted vs observed S(q) changes with q.
    The right panel shows how the proportion of predictors change among grid cells with q.
    Inputs:
    out_dir - directory where the results will be saved.
    q_list, r2_list - lists of q's and corresponding r2's
    pred_vars - a list of predictor variables
    pred_var_prop_list - a list of lists, each sublist contains the proportions (length q) for a predictor variable
    
    """
    fig = plt.figure(figsize = (8, 4))
    ax1 = plt.subplot(1, 2, 1) 
    plt.plot(q_list, r2_list, c = 'black', linewidth = 2)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('R-squared between observed and predicted S(q)', fontsize = 8)
    
    ax2 = plt.subplot(1, 2, 2)
    line_list_constraint = [None] * len(pred_vars)
    color = cm.rainbow(np.linspace(0, 1, len(pred_vars)))
    for i, pred_var_prop in enumerate(pred_var_prop_list):
        line_list_constraint[i],  = plt.plot(q_list, pred_var_prop, c = color[i], linewidth = 2)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('Proportion of grid cells', fontsize = 8)
    plt.ylim((0, 0.7))
    plt.legend(line_list_constraint, pred_vars, loc = 2, prop = {'size': 10})
    
    plt.subplots_adjust(wspace = 0.29)
    plt.savefig(out_dir, dpi = 600)    
    
def obtain_mean_annual_ndvi(in_folder, start_year, end_year, out_dir):
    """Mean annual NDVI across the grid using monthly data from multiple years.
    
    This is a highly specific function that only works with NDVI data in the same form
    as those from GIMMS. 
    """
    start_year, end_year = int(start_year), int(end_year)
    raster_list = []
    nodata_list = [-99, -88, -77]
    # First sum over months within each year, then take average across years
    for year in range(start_year, end_year + 1):
        year_list = []
        for month in range(1, 13):
            if len(str(month)) == 1: month_str = '0' + str(month)
            else: month_str = str(month)
            file_name = 'gimms_ndvi_qd_' + str(year) + month_str + '00.asc'
            file_array = import_raster_as_array(in_folder + file_name)
            for item in nodata_list: file_array[file_array == item] = float('nan')
            year_list.append(file_array)
        year_sum = np.nansum(np.array(year_list), axis = 0)
        raster_list.append(year_sum)
    mean_raster = np.nanmean(np.array(raster_list), axis = 0)
    mean_raster[np.isnan(mean_raster)] = -99
    
    last_file = gdal.Open(in_folder + file_name)
    proj = osr.SpatialReference()
    proj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    proj_wkt = proj.ExportToWkt()
    geotrans = last_file.GetGeoTransform()
    outRaster = write_raster_to_file(out_dir, len(mean_raster[0]), len(mean_raster), geotrans, proj_wkt, nodata = -99)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(mean_raster)
    outRaster.FlushCache()
    outRaster = None
    return None

def obtain_monthly_avg_ndvi(in_folder, start_year, end_year, out_folder):
    """NDVI for each month across multiple years.
    The output has the file name "NDVI_month_WGS84.tif', with "month" running from "01" to "12".
    This is a highly specific function that only works with NDVI data in the same form
    as those from GIMMS. 
    
    """
    start_year, end_year = int(start_year), int(end_year)
    nodata_list = [-99, -88, -77]
    for month in range(1, 13):
        if len(str(month)) == 1: month_str = '0' + str(month)
        else: month_str = str(month)
        raster_list = []
        for year in range(start_year, end_year + 1):
            file_name = 'gimms_ndvi_qd_' + str(year) + month_str + '00.asc'
            file_array = import_raster_as_array(in_folder + file_name)
            file_raster = gdal.Open(in_folder + file_name)
            file_nodata = file_raster.GetRasterBand(1).GetNoDataValue()
            for item in nodata_list: file_array[file_array == item] = float('nan')
            raster_list.append(file_array)
        month_mean = np.nanmean(np.array(raster_list), axis = 0)
        month_mean[np.isnan(month_mean)] = file_nodata
        
        out_dir_month = out_folder + 'NDVI_' + month_str + '_WGS84.tif'
        proj = osr.SpatialReference()
        proj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        proj_wkt = proj.ExportToWkt()
        geotrans = file_raster.GetGeoTransform()        
        outRaster = write_raster_to_file(out_dir_month, len(month_mean[0]), len(month_mean), geotrans, proj_wkt, nodata = file_nodata)
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(month_mean)
        outRaster.FlushCache()
        outRaster = None
    return None

def reproj_monthly_file(in_folder, in_file_1, in_file_2, out_folder, out_file_name, match_file, with_zero = False):
    """Reproject the monthly average environmental rasters to the designated projection and resolution.
    
    The input file has the path in_folder + in_file_1 + month + in_file_2
    The output file has the path out_folder + out_file_name + month + '.tif'
    
    """
    for month in range(1, 13):
        if len(str(month)) == 1: month_str = '0' + str(month)
        else: month_str = str(month)
        if with_zero is False:
            in_dir = in_folder + in_file_1 + str(month) + in_file_2
        else: in_dir = in_folder + in_file_1 + month_str + in_file_2
        out_dir = out_folder + out_file_name + month_str + '.tif'
        reproj_raster_to_match(in_dir, out_dir, match_file)
    return None

def obtain_annual_monthly_max(in_folder, in_file_1, in_file_2, out_folder, out_name, with_zero = False):
    """Given 12 monthly layers of one environmental variable, 
    
    create a raster with the max value across months for each cell, as well as a raster
    that stores the month identifier (1-12) for each cell. 
    The input file is assumed to have the path in_folder + in_file_1 + month (with or without zero) + in_file_2
    The output environemntal raster has the path out_folder + out_name + '_max.tif'
    The output month identifier raster has the path out_folder + out_name + '_max_month.tif'
    
    """
    monthly_list = []
    for month in range(1, 13):
        if len(str(month)) == 1: month_str = '0' + str(month)
        else: month_str = str(month)
        if with_zero is False:
            in_dir = in_folder + in_file_1 + str(month) + in_file_2
        else: in_dir = in_folder + in_file_1 + month_str + in_file_2
        month_file = gdal.Open(in_dir)
        month_array = import_raster_as_array(in_dir)
        monthly_list.append(month_array)
 
    monthly_max = np.nanmax(np.array(monthly_list), axis = 0)
    max_month = np.nanargmax(np.array(monthly_list), axis = 0) + 1
    
    # Replace cells with no data values with 0 for month
    nodata = month_file.GetRasterBand(1).GetNoDataValue()
    max_month[monthly_max == nodata] = 0
     
    # Write to file
    proj_wkt = month_file.GetProjection()
    geotrans = month_file.GetGeoTransform()
    out_dir_max = out_folder + out_name + '_max.tif'
    outRaster_max = write_raster_to_file(out_dir_max, len(monthly_max[0]), len(monthly_max), geotrans, proj_wkt, nodata = nodata)
    outband_max = outRaster_max.GetRasterBand(1)
    outband_max.WriteArray(monthly_max)
    outRaster_max.FlushCache()
    outRaster_max = None
    
    out_dir_month = out_folder + out_name + '_max_month.tif'
    outRaster_month = write_raster_to_file(out_dir_month, len(max_month[0]), len(max_month), geotrans, proj_wkt, nodata = 0)
    outband_month = outRaster_month.GetRasterBand(1)
    outband_month.WriteArray(max_month)
    outRaster_month.FlushCache()
    outRaster_month = None
    return None

def obtain_max_min_var(in_folder, in_file_1, in_file_2, out_folder, out_name, max_month_dir, max = True, with_zero = False):
    """Given 12 monthly rasters of one environmental variable, obtain a new raster corresponding to the months given by the raster
    
    at max_month_dir. 
    Inputs:
    in_folder, in_file_1, in_file_2: inputs to define input file paths, which are in the format in_folder + in_file_1 + month + in_file_2.
    out_folder: folder to which the output will be saved.
    max_month_dir: path of the raster given the desirable months ("max" months). Should be in same resolution and projection as input files.
    out_folder, out_name: defines the output path, which is in the format out_folder + out_name + '_' + max or min + '.tif'
    max: whether the output is the max (corresponding to the given months) or the min (corresponding to 6 months from the max month)
    with_zero: whether "month" in the input files contains zeros or not
    
    """
    monthly_var_list = []
    for month in range(1, 13):
        if len(str(month)) == 1: month_str = '0' + str(month)
        else: month_str = str(month)
        if with_zero is False:
            in_dir = in_folder + in_file_1 + str(month) + in_file_2
        else: in_dir = in_folder + in_file_1 + month_str + in_file_2
        monthly_var_list.append(import_raster_as_array(in_dir))
    
    monthly_var_list = np.array(monthly_var_list)
    var_file = gdal.Open(in_dir)
    nodata = var_file.GetRasterBand(1).GetNoDataValue()
    month_file = gdal.Open(max_month_dir)
    nodata_month = month_file.GetRasterBand(1).GetNoDataValue()
    
    month_array = import_raster_as_array(max_month_dir)
    if max is not True:  # Minimum instead of maximum month
        month_array[month_array == nodata_month] = float('nan')
        month_array = (month_array + 6)%12 
        month_array[month_array == 0] = 12
        month_array[np.isnan(month_array)] = nodata_month
    # Self-reminder below: k is "column", j is "row"
    out = [[monthly_var_list[month_array[j][k] - 1][j][k] if month_array[j][k] is not nodata_month \
            else nodata for k in range(len(month_array[0]))] for j in range(len(month_array))]
    out = np.array(out)
    # Write to file
    proj_wkt = var_file.GetProjection()
    geotrans = var_file.GetGeoTransform()
    if max: out_dir = out_folder + out_name + '_max.tif'
    else: out_dir = out_folder + out_name + '_min.tif'
    outRaster = write_raster_to_file(out_dir, len(out[0]), len(out), geotrans, proj_wkt, nodata = nodata)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(out)
    outRaster.FlushCache()
    outRaster = None
    return None

def corr_richness_quartiles(sp_list_flat, sp_sort):
    """Sub-function to compute the correlation between overall S & partial S of each quartile, given:
    
    sp_list_flat - 1-d (flattened) array with species list for each grid cell
    sp_sort - species list sorted by range size (from large to small)
    
    Output is a list with four r-values.
    
    """
    S_flat = [len(grid) for grid in sp_list_flat]
    Nquart = int(np.round(len(sp_sort) / 4))
    out = []
    for i in range(4):
        sp_quart = sp_sort[(Nquart * i): (min(Nquart * (i + 1), len(sp_sort)))]
        S_quart = [len([x for x in grid if x in sp_quart]) for grid in sp_list_flat]
        out.append(pearsonr(S_flat, S_quart)[0])
    return out
    
def corr_richness_taxon_continent(sp_range_array, sp_list_array, continent, out_dir, cont_geom = None):
    """Obtain the correlation between overall richness and partial richness for each of the four quartiles
    
    defind by range size on the continent.
    Inputs:
    sp_range_array: a structured array with the first column ('sp') as sp name, the second column ('global') 
        as the global range size of the sp, and subsequent columns as range of the sp on different continents
    continent - string, name of the continent. Can take one of the five values: 'global', 'North America', 'South America', 
        'Africa', 'Eurasia'.
    out_dir - output directory in 
    cont_geom - multipolygon of the continent in WKT, WGS 84. If coninent is 'global', it does not have to be defined.
    
    Saves the result to out_dir with five values in a row: continent, pearson correlation between overall S and partial S 
        for each of the four quartiles. 
        
     """
    if continent is not 'global':
        cont_geom_reproj = reproj_geom(cont_geom)
        cont_array = create_array_for_raster(proj_extent('behrmann'), geom = cont_geom_reproj)
        sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if cont_array[i][j] == 1]
    else: sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array))]
    sp_list_flat = [grid for grid in sp_list_flat if len(grid) > 0] # Remove empty grid cells with no species from the given taxon.
    
    sp_cont = sp_range_array[['sp', continent]][sp_range_array[continent] > 0]
    sp_sort = sp_cont['sp'][sp_cont[continent].argsort()[::-1]]
    quart_r = corr_richness_quartiles(sp_list_flat, sp_sort)
    
    out_file = open(out_dir, 'ab')
    print>>out_file, ','.join(map(str, [continent] + quart_r))
    out_file.close()

def gen_sp_range_contiguous(n_sp, sp_cont_cells, wrapped = True, ymax = None):
    """Sub-function to generate range for a species using the spreading dye algorithm with four neighbours."""
    xnew, ynew = random.sample(sp_cont_cells, 1)[0]
    sp_range = set([(xnew, ynew)])
    neighbours = set([])
    while len(sp_range) < n_sp:
        new_neighbour = [(xnew + 1, ynew), (xnew - 1, ynew)]
        if ynew == ymax: 
            new_neighbour += [(xnew, 0), (xnew, ynew  -1)]
        elif ynew == 0: 
            new_neighbour += [(xnew, ynew + 1), (xnew, ymax)]
        else: 
            new_neighbour += [(xnew, ynew + 1), (xnew, ynew - 1)]
        new_neighbour = [xy for xy in new_neighbour if (xy in sp_cont_cells and xy not in sp_range)]
        neighbours.update(new_neighbour)
        if len(neighbours) > 0: 
            xnew, ynew = random.sample(neighbours, 1)[0]
            neighbours.remove((xnew, ynew))
        else:  # when the neighbours are exhausted, jump to a new location (e.g., island)
            cont_remain = [cell for cell in sp_cont_cells if cell not in sp_range]
            xnew, ynew = random.sample(cont_remain, 1)[0]
        sp_range.update([(xnew, ynew)])
    return sp_range

def corr_null_range(sp_range_array, sp_list_array, continent, out_dir, sim_type = 'scattered', cont_geom = None, Niter = 200):
    """Randomly distribute species on the landscape, keeping the range size of each species, and compute correlation
    
    between overall S and partial S on these simulated landscapes.
    Species ranges are either scattered, or contiguous (spreading dye algorithm with four neighbours).
    
    """
    if continent is not 'global':
        cont_geom_reproj = reproj_geom(cont_geom)
        cont_array = create_array_for_raster(proj_extent('behrmann'), geom = cont_geom_reproj)
        sp_cont_cells = [(i, j) for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if (cont_array[i][j] == 1 and len(sp_list_array[i][j]) > 0)]
        sp_list_flat = [sp_list_array[i][j] for (i, j) in sp_cont_cells]
    else: 
        sp_cont_cells = [(i, j) for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) if len(sp_list_array[i][j]) > 0]
        sp_list_flat = [sp_list_array[i][j] for (i, j) in sp_cont_cells]
        
    sp_cont = sp_range_array[['sp', continent]][sp_range_array[continent] > 0] # Species with range on continent
    sp_sort = sp_cont['sp'][sp_cont[continent].argsort()[::-1]]
    sp_list_onelist = [sp for grid in sp_list_flat for sp in grid]
    sp_Ncell = Counter(sp_list_onelist)
    
    for i in range(Niter):
        sp_list_sim = [[] for grid in sp_list_flat]
        for sp in sp_Ncell.keys():
            n_sp = sp_Ncell[sp]
            if sim_type is 'scattered':
                sp_range_sim = np.random.choice(range(len(sp_list_flat)), n_sp, replace = False)
                sp_list_sim = [sp_list_sim[j] + [sp] if j in sp_range_sim else sp_list_sim[j] for j in range(len(sp_list_sim))]
            else:
                sp_range_sim = gen_sp_range_contiguous(n_sp, sp_cont_cells, ymax = len(sp_list_array[0]) - 1)
                sp_list_sim = [sp_list_sim[j] + [sp] if sp_cont_cells[j] in sp_range_sim else sp_list_sim[j] for j in range(len(sp_list_sim))]
        out_sim = corr_richness_quartiles(sp_list_sim, sp_sort)
        out_file = open(out_dir, 'ab')
        print>>out_file, ','.join(map(str, [continent] + out_sim))
        out_file.close()
        
def corr_sq_s_continent(sp_range_array, sp_list_array, continent, q, cont_geom = None):
    """Obtain the correlation between richness (S) and S(q) on continent or global scale. 
    
    If computed for continent, the continent-level range sizes are used. 
    S(q) is computed as S(q) = S_tot**(-q) * sum((range_i / range_tot) ** (-q))
    Cells with zero richness are removed.
    The output is a dictionary with two correlation values, 'spearman' and 'pearson'. 
    
     """
    if continent is not 'global':
        cont_geom_reproj = reproj_geom(cont_geom)
        cont_array = create_array_for_raster(proj_extent('behrmann'), geom = cont_geom_reproj)
        sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                    if (cont_array[i][j] == 1 and len(sp_list_array[i][j]) > 0)]
    else: sp_list_flat = [sp_list_array[i][j] for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array)) \
                          if len(sp_list_array[i][j]) > 0]
    sp_unique = np.unique([sp for grid in sp_list_flat for sp in grid])
    
    sp_cont = sp_range_array[['sp', continent]][sp_range_array[continent] > 0]
    sp_range_dic = {}
    for row in sp_cont: sp_range_dic[row[0]] = row[1]
    for sp in sp_unique: # Some species with limited ranges are missed by the continent shapefile
        if sp not in sp_range_dic.keys():
            sp_range_dic[sp] = sp_range_array['global'][sp_range_array['sp'] == sp][0]
    range_tot = sum(sp_range_dic.values())
    
    S_list = [len(grid) for grid in sp_list_flat]
    Sq_list = []
    for grid in sp_list_flat: # Note here the constant is S_tot ** (-q), not S ** (-q), which wouldn't be a constant
        Sq = len(sp_cont) ** (-q) * np.sum([(range_tot / sp_range_dic[sp]) ** q for sp in grid])
        Sq_list.append(Sq)
    corr_dic = {}
    corr_dic['spearman'] = spearmanr(Sq_list, S_list)[0]
    corr_dic['pearson'] = pearsonr(Sq_list, S_list)[0]
    return corr_dic

def model_sq(sp_range_array, sp_list_array, continent, AET_array, spatT_array, tempT_array, 
             q, cont_geom = None):
    """Model S(q), globally and on each continent, using environmental variables. 
    
    Current model: S(q) ~ AET + spatT * tempT
    S(q) is log-transformed, the indep vars are square-root-transformed. All are then z-normalized.
    
    Return: 
    
     """
    # Note that for other variables the condition >=0 to screen nodata values may have to change
    if continent is not 'global':
        cont_geom_reproj = reproj_geom(cont_geom)
        cont_array = create_array_for_raster(proj_extent('behrmann'), geom = cont_geom_reproj)
        valid_index = [1 if (cont_array[i][j] == 1 and len(sp_list_array[i][j]) > 0 and AET_array[i][j] >= 0 and \
                           spatT_array[i][j] >= 0 and tempT_array[i][j] >= 0) else 0 \
                           for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array))]
    else:
        valid_index = [1 if (len(sp_list_array[i][j]) > 0 and AET_array[i][j] >= 0 and \
                           spatT_array[i][j] >= 0 and tempT_array[i][j] >= 0) else 0 \
                           for j in range(len(sp_list_array[0])) for i in range(len(sp_list_array))]
    sp_list_flat = [x for (index, x) in zip(valid_index, list(sp_list_array.flatten('F'))) if index > 0]
    AET_flat = rob.FloatVector([x for (index, x) in zip(valid_index, list(AET_array.flatten('F'))) if index > 0])
    spatT_flat = rob.FloatVector([x for (index, x) in zip(valid_index, list(spatT_array.flatten('F'))) if index > 0])
    tempT_flat = rob.FloatVector([x for (index, x) in zip(valid_index, list(tempT_array.flatten('F'))) if index > 0])
    sp_unique = np.unique([sp for grid in sp_list_flat for sp in grid])
    
    sp_cont = sp_range_array[['sp', continent]][sp_range_array[continent] > 0]
    sp_range_dic = {}
    for row in sp_cont: sp_range_dic[row[0]] = row[1]
    for sp in sp_unique: # Some species with limited ranges are missed by the continent shapefile
        if sp not in sp_range_dic.keys():
            sp_range_dic[sp] = sp_range_array['global'][sp_range_array['sp'] == sp][0]
    range_tot = sum(sp_range_dic.values())
    
    Sq_list = []
    for grid in sp_list_flat: # Note here the constant is S_tot ** (-q), not S ** (-q), which wouldn't be a constant
        Sq = len(sp_cont) ** (-q) * np.sum([(range_tot / sp_range_dic[sp]) ** q for sp in grid])
        Sq_list.append(Sq)
    
    Sq_list = rob.FloatVector(Sq_list)
    out = rfunctions.sq_regression(Sq_list, AET_flat, spatT_flat, tempT_flat)
    out = [x for x in out]
    return out
