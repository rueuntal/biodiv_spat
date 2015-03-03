from __future__ import division
import os
import numpy as np
from scipy import stats
from osgeo import ogr, osr, gdal
import psycopg2
import shapely.wkt, shapely.ops

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

def convert_array_to_raster(array, rasterOrigin, out_file, pixel_size, \
                            out_proj4 = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'):
    """Convert an array to raster."""
    cols = array.shape[1]
    rows = array.shape[0]
    xmin, ymax = rasterOrigin
    outRaster = gdal.GetDriverByName('GTiff').Create(out_file, cols, rows, 1, gdal.GDT_Int16)
    outRaster.SetGeoTransform((xmin, pixel_size, 0, ymax, 0, -pixel_size))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(out_proj4)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())    
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outband.SetNoDataValue(0)
    outRaster.FlushCache()
    outRaster = None
    return None
    
def convert_sp_array(postgis_cur, table_name, sp, pixel_size):
    """Subfunction to obtain a raster array of presence/absence for a given species."""
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
    # Convert species range to raster array
    sp_landscape = create_array_for_raster(proj_extent('behrmann'), geom = wkt_reproj, pixel_size = pixel_size)
    return sp_landscape
    
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
                sp_landscape =  convert_sp_array(postgis_cur, table_name, sp, pixel_size)
                family_landscape += sp_landscape
            convert_array_to_raster(family_landscape, [xmin, ymax], out_file, pixel_size)
    return None

def create_array_sp_list(postgis_cur, table_name, out_dir, pixel_size = 100000):
    """Create an array with species list in each grid covering the globe."""
    xmin, xmax, ymin, ymax = proj_extent('behrmann')
    sp_list_exec = 'SELECT DISTINCT binomial FROM ' + table_name
    postgis_cur.execute(sp_list_exec)
    sp_list = [x[0] for x in postgis_cur.fetchall()] 
    
    x_res = int((xmax - xmin) / pixel_size)
    y_res = int((ymax - ymin) / pixel_size)
    array_list = np.array([[[] for i in range(x_res)] for j in range(y_res)])
    
    for sp in sp_list:
        sp_array = convert_sp_array(postgis_cur, table_name, sp, pixel_size)
        array_list = np.array([[list(array_list[j][i]) + [sp] if sp_array[j][i] > 0 else list(array_list[j][i]) for i in range(x_res)] for j in range(y_res)])

    