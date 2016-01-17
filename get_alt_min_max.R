.libPaths(c(.libPaths(), "C:/Users/Xiao/Documents/R/win-library/3.1"))

library(raster)
dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\pred_vars\\'
in_file = paste(dir, 'alt_behrmann_1000.tif', sep = '')

alt_behrmann = raster(in_file)
out_file_min = paste(dir, 'alt_behrmann_min.tif', sep = '')
out_file_max = paste(dir, 'alt_behrmann_max.tif', sep = '')

alt_min = aggregate(alt_behrmann, 100, fun = min, expand = F, na.rm = T)
writeRaster(alt_min, out_file_min, overwrite = T, NAflag = -9999)
alt_max = aggregate(alt_behrmann, 100, fun = max, expand = F, na.rm = T)
writeRaster(alt_max, out_file_max, overwrite = T, NAflag = -9999)
