from __future__ import division
import os
import numpy as np
from scipy import stats
import cPickle
from osgeo import ogr, osr, gdal
import psycopg2
import shapely.wkt, shapely.ops
import multiprocessing

def import_pickle_file(in_dir):
    """Read in pickled file."""
    in_file = open(in_dir, 'rb')
    in_obj = cPickle.load(in_file)
    in_file.close()
    return in_obj

def import_shapefile(in_dir): 
    """Read in a .shp file and save the geometry each feature in a list of WKT."""
    driver = ogr.GetDriverByName('ESRI Shapefile')
    datasource = driver.Open(in_dir, 0)
    layer = datasource.GetLayer()
    geom_list = []
    for feature in layer:
        geom_list.append(feature.GetGeometryRef().ExportToWkt())
    return geom_list

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

def proj_extent(proj_name):
    """Give the global extent (xmin, xmax, ymin, ymax) of a projection for raster. Only Behrmann is available now."""
    in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' # Lat and long
    if proj_name == 'behrmann':
        out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
    
    in_proj = osr.SpatialReference()
    in_proj.ImportFromProj4(in_proj4)
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(out_proj4)
    transform = reproj(in_proj4 = in_proj4, out_proj4 = out_proj4)
    
    x, y, z = transform.TransformPoint(180, 90)
    return [-x, x, -y, y]

def weighted_sample_range_size(range_size_list, size):
    """Generate a random sample of a certain size (length) from a list of range sizes,
    
    weighted by the range sizes.
    
    """
    sum_size = np.sum(range_size_list)
    rand_unif = np.sort(stats.uniform.rvs(0, sum_size, size = size))
    partial_sum = 0
    selected_range = []
    i, j = 0, 0
    while len(selected_range) < size:
        partial_sum += range_size_list[i]
        while j < len(rand_unif) and partial_sum >= rand_unif[j]:
            selected_range.append(range_size_list[i])
            j += 1
        i += 1
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
    return band.ReadAsArray()

def convert_array_to_raster(array, rasterOrigin, out_file, pixel_size, no_value = 0, \
                            out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Convert an array to raster."""
    cols = array.shape[1]
    rows = array.shape[0]
    xmin, ymax = rasterOrigin
    geotrans = (xmin, pixel_size, 0, ymax, 0, -pixel_size)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(out_proj4)
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

def sp_reproj_birds(folder, file_name):
    """sp_reproj() for birds. This function returns two strings - species binomial, and its reprojected range."""
    in_dir = folder + '/' + file_name
    file_split = file_name.split('_')
    sp_name = file_split[0] + ' ' + file_split[1]
    sp_driver = ogr.GetDriverByName('ESRI Shapefile')
    sp_datasource = sp_driver.Open(folder + '/' + file_name, 0)
    sp_layer = sp_datasource.GetLayer()
    sp_feature = sp_layer[0]
    sp_geom_wkt = sp_feature.GetGeometryRef().ExportToWkt()
    wkt_reproj = reproj_geom(sp_geom_wkt) # Reproject to Behrmann
    return sp_name, wkt_reproj

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

def create_sp_range_dic_bird(folder, out_file_name):
    """create_sp_range_dic() for birds, where there is one shapefile per species instead of one shapefile for all species combined."""
    sp_range_dic = {}
    for file in os.listdir(folder):
        if file.endswith('.shp'):
            sp_name, wkt_reproj = sp_reproj_birds(folder, file)
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

def create_array_sp_list_birds(folder, out_file_name, pixel_size = 100000):
    """create_array_sp_list() for birds."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    array_list = np.array([[[] for i in range(x_res)] for j in range(y_res)])
    
    file_list = os.listdir(folder)
    for file in file_list:
        if file.endswith('.shp'):
            sp, wkt_reproj = sp_reproj_birds(folder, file)
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
    metric: metric  to compute. See docstring for compare_range_size_dists() for details.    full_dist_expc: the expected values from the weighted full distribution of the same length as dist, only needed if metric = "ks".
    
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
    The output file is written to the designated output directory, with the same out_name + "_" + metric + ".pck".
    
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

def reproj_raster_to_match(in_dir, out_dir, match_dir, alg = gdal.GRA_Bilinear, \
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
  
    match_file = gdal.Open(match_dir)
    match_proj_wkt = match_file.GetProjection()
    match_geotrans = match_file.GetGeoTransform()
    wide = match_file.RasterXSize
    high = match_file.RasterYSize
    
    outRaster = write_raster_to_file(out_dir, wide, high, match_geotrans, match_proj_wkt, nodata = in_nodata)
    res = gdal.ReprojectImage(in_file, outRaster, in_proj_wkt, match_proj_wkt, alg)
    outRaster.FlushCache()
    outRaster = None    
    
    return None
