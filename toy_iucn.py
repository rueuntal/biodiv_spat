from __future__ import division
import os
import numpy as np
from scipy import stats
import cPickle
from osgeo import ogr, osr, gdal
import psycopg2
import shapely.wkt, shapely.ops

def import_pickle_file(in_dir):
    """Read in pickled file."""
    in_file = open(in_dir, 'rb')
    in_obj = cPickle.load(in_file)
    in_file.close()
    return in_obj

def reproj(in_geom, in_proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', \
           out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Function to reproject geometries (defined by WKT). 
    
    Default input is unprojected WGS84 (EPSG 4326), default output is Behrmann equal area.
    
    """
    in_geom_ogr = ogr.CreateGeometryFromWkt(in_geom)
    in_proj = osr.SpatialReference()
    in_proj.ImportFromProj4(in_proj4)
    out_proj = osr.SpatialReference()
    out_proj.ImportFromProj4(out_proj4)
    transform = osr.CoordinateTransformation(in_proj, out_proj)
    in_geom_ogr.Transform(transform)
    return in_geom_ogr.ExportToWkt()

def proj_extent(proj_name):
    """Give the extend (xmin, xmax, ymin, ymax) of a projection for raster. Only Behrmann is available now."""
    if proj_name == 'behrmann':
        # Obtained from http://tiles.arcgis.com/tiles/BG6nSlhZSAWtExvp/arcgis/rest/services/coordsys_Behrmann/MapServer?f=pjson
        xmin = -1.73675293395E7
        ymin = -7342230.1365
        xmax = 1.73675293395E7
        ymax = 7342230.136500001
    return [xmin, xmax, ymin, ymax]

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
    geom: Default is empty, leading to an array filled with non-values. Or could be a multipolygon (WKT), leading to 1's on the landscape.
    no_value: value to fill the grid with no data. Default is 0.
    pixel_size: size of each pixel. Default is 100km. 
    
    """
    xmin, xmax, ymin, ymax = extent
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    target_ds = gdal.GetDriverByName('MEM').Create('', x_res, y_res, gdal.GDT_Byte)
    target_ds.SetGeoTransform((xmin, pixel_size, 0, ymax, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(no_value)
    if geom: 
        # Convert geom to layer
        virtual_shp = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('/vsimem/virtual.shp')
        virtual_layer = virtual_shp.CreateLayer('sp_range', geom_type = ogr.wkbPolygon)
        idField = ogr.FieldDefn('id', ogr.OFTInteger)
        virtual_layer.CreateField(idField)
        featureDefn = virtual_layer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(ogr.CreateGeometryFromWkt(geom))
        feature.SetField('id', 1)
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
    outRaster = gdal.GetDriverByName('GTiff').Create(out_file, cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform((xmin, pixel_size, 0, ymax, 0, -pixel_size))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(out_proj4)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())    
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outband.SetNoDataValue(no_value)
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
    wkt_reproj = reproj(sp_geom_wkt)
    return wkt_reproj
    
def richness_to_raster_by_family(postgis_cur, table_name, out_dir, pixel_size = 100000):
    """Convert an IUCN shapefile with range maps of a taxon into rasters of richness, one for each family, 
    
    assuming that the shapefile has already been loaded into a PostGIS database.
    Inputs:
    postgis_cur - cursor pointing to the PostGIS database
    table_name - name of the table in the PostGIS databse
    out_dir - output directory for the resulting raster
    pixel_size - grain size of raster, default is 100km * 100km
    
    """
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    family_list_exec = 'SELECT DISTINCT family_nam FROM ' + table_name
    postgis_cur.execute(family_list_exec)
    family_list = [x[0] for x in postgis_cur.fetchall()]
    
    for family in family_list:
        out_file = out_dir + '/' + table_name + '_' + family + '_' + str(int(pixel_size)) + '.tif'
        # Check to make sure file doesn't exist already - useful to resume from unexpected breaks
        if not os.path.exists(out_file):
            sp_list_exec = "SELECT DISTINCT binomial FROM " + table_name + " WHERE family_nam='" + family + "'"
            postgis_cur.execute(sp_list_exec)
            sp_list = [x[0] for x in postgis_cur.fetchall()]
            family_landscape = create_array_for_raster(proj_extent('behrmann'), pixel_size = pixel_size)
            # Create empty array 
            for sp in sp_list:
                wkt_reproj = sp_reproj(postgis_cur, table_name, sp)
                # Convert species range to raster array
                sp_landscape = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)        
                family_landscape += sp_landscape
            convert_array_to_raster(family_landscape, [xmin, ymax], out_file, pixel_size)
    return None

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
            file_split = file.split('_')
            sp_name = file_split[0] + ' ' + file_split[1]
            if sp_name in sp_range_dic: print "Warning: " + sp_name + " has duplicates."
            sp_driver = ogr.GetDriverByName('ESRI Shapefile')
            sp_datasource = sp_driver.Open(folder + '/' + file, 0)
            sp_layer = sp_datasource.GetLayer()
            sp_feature = sp_layer[0]
            sp_geom_wkt = sp_feature.GetGeometryRef().ExportToWkt()
            wkt_reproj = reproj(sp_geom_wkt) # Reproject to Behrmann
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
    
def compare_range_size_dists(sp_list_array, range_size_dic, out_dir, out_name, Nsample, metric, threshold = 5):
    """Compare with empirical range size distribution in each grid with the expected distribution, 
    
    and save the output into an array.
    
    Inputs:
    sp_list_array: array with species list in each grid, created by create_array_sp_list().
    range_size_dic: dictionary of range size for each species, created by create_sp_range_dic().
    out_dir: output directory
    out_name: name of the output file, e.g., "terrestrial_mammals", "reptiles", etc. 
    Nsample: number of random samples to draw from the full distribution.
    metric: metric used to compare the empirical and randomly generated size distributions, which can take the following values:
    "mean", "sd" (standard deviation), "skew" (skewness). All values are calcuated on log scale.
    threshold: minimal species richness of a grid, below which it is not analyzed.
    
    Output:
    File with an array of the same dimension as sp_list_array, where the value of each grid is the quantile of the empirical size distribution among the 
    randomly generated size distributions, measured with the designated metric. 
    The output file is written to the designated output directory, with the same out_name + "_" + metric + ".pck".
    
    """
    array_out = np.empty([len(sp_list_array), len(sp_list_array[0])], dtype = float)
    range_size_list = range_size_dic.values()
    for j in range(len(sp_list_array)):
        for i in range(len(sp_list_array[0])):
            sp_grid = sp_list_array[j][i]
            richness_grid = len(sp_grid)
            if richness_grid < threshold: array_out[j][i] = -1 # Fill unanalyzed grids with -1
            else:
                emp_range_dist_grid = [range_size_dic[sp] for sp in sp_grid]
                emp_metric = metric_dist(emp_range_dist_grid, metric)
                rand_metric_list = []
                for k in range(Nsample):
                    rand_dist_grid = weighted_sample_range_size(range_size_list, len(emp_range_dist_grid))
                    rand_metric_list.append(metric_dist(rand_dist_grid, metric))
                # Compute quantile
                quan = len([x for x in rand_metric_list if x <= emp_metric]) / Nsample
                array_out[j][i] = quan
    
        # Save to file
        out_file_name = out_dir + '/' + out_name + '_' + metric + '.pck'
        out_file = open(out_file_name, 'wb')
        cPickle.dump(array_out, out_file, protocol = 2)
        out_file.close()
        return None
    
