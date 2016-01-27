
.libPaths(c(.libPaths(), "C:/Users/Xiao/Documents/R/win-library/3.1"))
library(raster)

# The following input is for test purpose
# S_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\weighted_S\\amphibians_weighted_richness_100000_0.0.tif'

# This function reads in the raster from S_dir, and conducts linear regression between 
# S and the environmental variables. 
# The full model is S ~ (temp + PET) + (AET + NDVI) + (alt range * seasonality).
# The path of the environmental rasters are fixed.
# It returns the following 14 values (in this order): 
# r^2 for the full, the temperature, the productivity and Janzen's models
# Unique contribution of temperature, productivity, and Janzen's effect
# p-value of the seven variables in the full model
multilin = function(S_dir){
  dat = c(as.matrix(raster(S_dir)))
  pred_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\pred_vars\\'
  pred_names = c('mean_annual_T', 'PET', 'AET', 'NDVI', 'alt_range', 'seasonality')
  for (pred in pred_names){
    pred_file = paste(pred_dir, pred, '.tif', sep = '')
    pred_raster = c(as.matrix(raster(pred_file)))
    dat = cbind(dat, pred_raster)
  }
  dat = as.data.frame(dat)
  names(dat) = c('S', pred_names)
  # Only keep rows where all variables are available
  index = complete.cases(dat)
  dat = dat[index, ]
  
  out = vector(mode = 'numeric', length = 14)
  full_model = lm(log(S) ~ mean_annual_T + PET + AET + NDVI + alt_range * seasonality, data = dat)
  out[1] = summary(full_model)$r.squared
  out[8:14] = as.numeric(summary(full_model)$coef[2:8, 4])
  T_model = lm(log(S) ~ mean_annual_T + PET, data = dat)
  out[2] = summary(T_model)$r.squared
  no_T_model = lm(log(S) ~ AET + NDVI + alt_range * seasonality, data = dat)
  out[5] = summary(full_model)$r.squared - summary(no_T_model)$r.squared
  prod_model = lm(log(S) ~ AET + NDVI, data = dat)
  out[3] = summary(prod_model)$r.squared
  no_prod_model = lm(log(S) ~ mean_annual_T + PET + alt_range * seasonality, data = dat)
  out[6] = summary(full_model)$r.squared - summary(no_prod_model)$r.squared
  Janzen_model = lm(log(S) ~ alt_range * seasonality, data = dat)
  out[4] = summary(Janzen_model)$r.squared
  no_Janzen_model = lm(log(S) ~ mean_annual_T + PET + AET + NDVI, data = dat)
  out[7] = summary(full_model)$r.squared - summary(no_Janzen_model)$r.squared
  return(out)
}

