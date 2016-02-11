from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
import cPickle
from osgeo import ogr, osr, gdal
import psycopg2
import shapely.wkt, shapely.ops
import multiprocessing
from matplotlib.mlab import PCA
import re
from math import sin, cos, sqrt, atan2
import csv
import statsmodels.api as sm

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

def import_raster_as_array(raster_dir):
    raster = gdal.Open(raster_dir)
    band = raster.GetRasterBand(1).ReadAsArray()
    return band

def range_dic_to_csv(in_dir, out_dir, remove_list = []):
    """Read in a dictionary of range sizes for a taxon and convert it to csv for processing in R."""
    range_dic = import_pickle_file(in_dir)
    out_write = open(out_dir, 'wb')
    out = csv.writer(out_write)
    
    for sp in range_dic.keys():
        if sp not in remove_list:
            results = np.zeros((1, ), dtype = ('S25, f15'))
            results['f0'] = sp
            results['f1'] = range_dic[sp]
            out.writerows(results)
    out_write.close()

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
    
def reproj(in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', \
           out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """General function for reprojection which can be used for both vectors and rasters."""
    in_proj = osr.SpatialReference()
    in_proj.ImportFromProj4(in_proj4)
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(out_proj4)
    transform = osr.CoordinateTransformation(in_proj, out_proj)
    return transform

def reproj_geom(in_geom, in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', \
           out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Function to reproject geometries (defined by WKT). 
    
    Default input is unprojected WGS84 (EPSG 4326), default output is Behrmann equal area.
    
    """
    in_geom_ogr = ogr.CreateGeometryFromWkt(in_geom)
    transform = reproj(in_proj4 = in_proj4, out_proj4 = out_proj4)
    in_geom_ogr.Transform(transform)
    return in_geom_ogr.ExportToWkt()

def proj_name_to_proj4(proj_name):
    """Given name of a projection, return its specification in PROJ4."""
    name_and_proj = {'latlong': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                     'behrmann': '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
                     }
    return name_and_proj[proj_name]

def proj_extent(proj_name):
    """Give the global extent (xmin, xmax, ymin, ymax) of a projection for raster. Only Behrmann is available now."""
    in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' # Lat and long
    out_proj4 = proj_name_to_proj4(proj_name)
    transform = reproj(in_proj4 = in_proj4, out_proj4 = out_proj4)
    
    x, y, z = transform.TransformPoint(180, 90)
    return [-x, x, -y, y]

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

def convert_array_to_raster(array, rasterOrigin, out_file, pixel_size, no_value = 0, \
                            out_proj = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Convert an array to raster."""
    cols = array.shape[1]
    rows = array.shape[0]
    xmin, ymax = rasterOrigin
    geotrans = (xmin, pixel_size, 0, ymax, 0, -pixel_size)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(out_proj)
    proj_wkt = outRasterSRS.ExportToWkt()
    
    outRaster = write_raster_to_file(out_file, cols, rows, geotrans, proj_wkt, nodata = no_value)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRaster.FlushCache()
    outRaster = None
    return None
    
def sp_reproj(postgis_cur, table_name, sp):
    """Subfunction to obtain the range map of a species and reproject it to Behrmann Equal Area."""
    sp_geom_exec = "SELECT ST_AsText(geom) FROM " + table_name + " WHERE binomial='" + sp + "'"
    postgis_cur.execute(sp_geom_exec)
    sp_geom_list = [x[0] for x in postgis_cur.fetchall()] # Here sp_geom_list is a list of geoms in WKT, each a multipolygon
    sp_geom_shapes = [shapely.wkt.loads(x) for x in sp_geom_list]
    try:
        sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
    except: # Fails for Gulo gulo, because of holes?
        sp_geom_shapes = [x.buffer(0) for x in sp_geom_shapes] # Adding zero buffer somehow helps
        sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
    # Reproject geom into Behrmann                
    wkt_reproj = reproj_geom(sp_geom_wkt)
    return wkt_reproj

def sp_reproj_birds(folder, file_name, Attr = None, Attr_filter = None):
    """sp_reproj() for birds. This function returns two strings - species binomial, and its reprojected range."""
    in_dir = folder + '/' + file_name
    file_split = file_name.split('_')
    sp_name = file_split[0] + ' ' + file_split[1]
    if isinstance(Attr_filter, basestring): # If Attr_filter is a string, i.e., a single value
        sp_geom_list = import_shapefile(in_dir, Attr = Attr, AttrFilter = Attr_filter)
    else: 
        sp_geom_list = []
        for ind_filter in Attr_filter: 
            sp_geom_list.extend(import_shapefile(in_dir, Attr = Attr, AttrFilter = ind_filter))
    if len(sp_geom_list) > 0: # If Attribute exists
        sp_geom_shapes = [shapely.wkt.loads(x) for x in sp_geom_list]
        try:
            sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
        except: 
            sp_geom_shapes = [x.buffer(0) for x in sp_geom_shapes]
            sp_geom_wkt = shapely.wkt.dumps(shapely.ops.cascaded_union(sp_geom_shapes))
        wkt_reproj = reproj_geom(sp_geom_wkt) 
        return sp_name, wkt_reproj
    else: return sp_name, None

def richness_to_raster(postgis_cur, table_name, out_dir, pixel_size = 100000, remove_sp_list = []):
    """Convert an IUCN shapefile with range maps of a taxon into a raster file of richness, with a given list of species removed."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    sp_list_exec = 'SELECT DISTINCT binomial FROM ' + table_name
    postgis_cur.execute(sp_list_exec)
    sp_list = [x[0] for x in postgis_cur.fetchall() if x[0] not in remove_sp_list]
    
    table_landscape = create_array_for_raster(proj_extent('behrmann'), pixel_size = pixel_size)
    for sp in sp_list:
            wkt_reproj = sp_reproj(postgis_cur, table_name, sp)
            # Convert species range to raster array
            sp_landscape = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)        
            table_landscape += sp_landscape
            
    out_file = out_dir + '/' + table_name + '_richness_' + str(int(pixel_size)) + '.tif'
    convert_array_to_raster(table_landscape, [xmin, ymax], out_file, pixel_size)
    return None
 
def richness_to_raster_birds(folder, out_dir, pixel_size = 100000, remove_sp_list = []):
    """richness_to_raster() for birds."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    table_landscape = create_array_for_raster(proj_extent('behrmann'), pixel_size = pixel_size)
    for file in os.listdir(folder):
        if file.endswith('.shp'):
            sp_name, wkt_reproj = sp_reproj_birds(folder, file)
            sp_landscape = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size) 
            table_landscape += sp_landscape
            
    out_file = out_dir + '/birds_richness_' + str(int(pixel_size)) + '.tif'
    convert_array_to_raster(table_landscape, [xmin, ymax], out_file, pixel_size)
    return None
 
def richness_to_raster_ver2(sp_array, out_dir, out_name, pixel_size = 100000, remove_sp_list = []):  
    """Another way to convert richess to raster using existing array of species list, which is much faster."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    table_richness = np.zeros(shape = sp_array.shape)
    for j in range(len(sp_array)):
        for i in range(len(sp_array[0])):
            sp_grid = sp_array[j][i]
            richness = len([sp for sp in sp_grid if sp not in remove_sp_list])
            table_richness[j][i] = richness
    out_file = out_dir + '/' + out_name + '_richness_' + str(pixel_size) + '.tif'
    convert_array_to_raster(table_richness, [xmin, ymax], out_file, pixel_size)
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

def create_sp_range_dic(postgis_cur, table_name, out_file_name):
    """Obtain the range size of species (in m^2), put in a dictionary, and save to file. 
    
    Range sizes are calculated under Behrmann equal area projection.
    
    """
    sp_list_exec = 'SELECT DISTINCT binomial FROM ' + table_name
    postgis_cur.execute(sp_list_exec)
    sp_list = [x[0] for x in postgis_cur.fetchall()] 
    sp_range_dic = {}
    
    for sp in sp_list:
        wkt_reproj = sp_reproj(postgis_cur, table_name, sp)
        sp_range_shape = shapely.wkt.loads(wkt_reproj)
        sp_range_dic[sp] = sp_range_shape.area
    
    out_file = open(out_file_name, 'wb')
    cPickle.dump(sp_range_dic, out_file, protocol = 2)
    out_file.close()
    return None        

def create_sp_range_dic_bird(folder, out_file_name, Attr = None, Attr_filter = None):
    """create_sp_range_dic() for birds, where there is one shapefile per species instead of one shapefile for all species combined."""
    sp_range_dic = {}
    for file in os.listdir(folder):
        if file.endswith('.shp'):
            sp_name, wkt_reproj = sp_reproj_birds(folder, file, Attr = Attr, Attr_filter = Attr_filter)
            if wkt_reproj is not None:
                if sp_name in sp_range_dic: print "Warning: " + sp_name + " has duplicates."
                sp_range_shape = shapely.wkt.loads(wkt_reproj)
                sp_range_dic[sp_name] = sp_range_shape.area
    
    out_file = open(out_file_name, 'wb')
    cPickle.dump(sp_range_dic, out_file, protocol = 2)
    out_file.close()
    return None        

def create_array_sp_list(postgis_cur, table_name, out_file_name, pixel_size = 100000):
    """Create an array with species list in each grid covering the globe."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    sp_list_exec = 'SELECT DISTINCT binomial FROM ' + table_name
    postgis_cur.execute(sp_list_exec)
    sp_list = [x[0] for x in postgis_cur.fetchall()] 
    
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    array_list = np.array([[[] for i in range(x_res)] for j in range(y_res)])
    
    for sp in sp_list:
        wkt_reproj = sp_reproj(postgis_cur, table_name, sp)
        # Convert species range to raster array
        sp_array = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)        
        array_list = np.array([[list(array_list[j][i]) + [sp] if sp_array[j][i] > 0 else list(array_list[j][i]) for i in range(x_res)] for j in range(y_res)])
    
    # Save to file
    out_file = open(out_file_name, 'wb')
    cPickle.dump(array_list, out_file, protocol = 2)
    out_file.close()
    return None

def create_array_sp_list_birds(folder, out_file_name, Attr = None, Attr_filter = None, pixel_size = 100000):
    """create_array_sp_list() for birds."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    array_list = np.array([[[] for i in range(x_res)] for j in range(y_res)])
    
    file_list = os.listdir(folder)
    for file in file_list:
        if file.endswith('.shp'):
            sp, wkt_reproj = sp_reproj_birds(folder, file, Attr = Attr, Attr_filter = Attr_filter)
            # Convert species range to raster array
            sp_array = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)        
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
    """Plot the weighing parameter q against R^2 values (full model, three reduced models, 
    
    and the unique contributions from the three sets of parameters.
    The output is a 1*2 plot, with R^2 of the four models in the left subplot, and unique contributions in the right subplot.
    Inputs:
    results_dir: txt file with the output from multiple regression, generated in weighted_richness.py
    taxon: string, name of the taxon
    out_dir: directory for the output figure
    legend: whether legend is to be included in the figure
    
    """
    results = np.genfromtxt(results_dir, dtype = None, delimiter = '\t')
    results_taxon = results[results['f0'] == taxon]
    q, tot_r2, T_r2, prod_r2, janzen_r2, T_r2_part, prod_r2_part, janzen_r2_part = \
        zip(*results_taxon[['f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8']])
    r2 = [tot_r2, T_r2, prod_r2, janzen_r2]
    unique_r2 = [T_r2_part, prod_r2_part, janzen_r2_part]
    col = ['#000000', '#FF4040', '#1874CD', '#8B7765']
    fig = plt.figure(figsize = (7, 4))
    ax1 = plt.subplot(1, 2, 1) 
    line_list = [None] * 4
    for i, r2_val in enumerate(r2):
        line_list[i],  = plt.plot(q, r2_val, c = col[i], linewidth = 2)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('R-squared', fontsize = 8)
    plt.ylim((0, 1))
    ax2 = plt.subplot(1, 2, 2)
    for j, unique_r2_val in enumerate(unique_r2):
        plt.plot(q, unique_r2_val, c = col[j + 1], linewidth = 2)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 6)
    plt.xlabel('Weighing parameter q', fontsize = 8)
    plt.ylabel('Unique R-squared', fontsize = 8)
    plt.ylim((0, 0.7))
    if legend:
        plt.legend(line_list, ['Full', 'Temperature', 'Productivity', 'Janzen'], loc = 2, prop = {'size': 10})
    
    plt.subplots_adjust(wspace = 0.29)
    plt.savefig(out_dir, dpi = 600)
    return None
 
def obtain_avg_ndvi(in_folder, start_year, end_year, out_dir):
    """Average NDVI across the grid using monthly data from multiple years.
    
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
