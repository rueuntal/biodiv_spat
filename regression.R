.libPaths(c(.libPaths(), "C:/Users/Xiao/Documents/R/win-library/3.2"))
library(raster)
library(quantreg)

pred_folder = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\pred_vars\\'

# The following input is for test purpose
# S_dir = 'C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\Janzen\\weighted_S\\amphibians_weighted_richness_100000_0.0.tif'

# This function reads in the raster from S_dir, and conducts linear regression between 
# S and the environmental variables. 
# The full model is S ~ (AET + NDVI) + (alt range * seasonality).
# The path of the environmental rasters are fixed.
# It returns the following 13 values (in this order): 
# r^2 for the full, the productivity and Janzen's models
# Unique contribution of productivity and Janzen's effect
# Unique contribution of altitudinal range, seasonality, and interaction in Janzen's model
# p-value of the five variables in the full model
multilin = function(S_dir, env = 'annual', transform = 'untransformed'){
  dat = c(as.matrix(raster(S_dir)))
  if (env == 'annual'){
    pred_names = c('AET', 'NDVI', 'alt_range', 'seasonality')
  } else if (env == 'breeding'){
    pred_names = c('mean_T_max', 'PET_max', 'AET_max', 'NDVI_max', 'alt_range', 'seasonality')
  } else {
    pred_names = c('mean_T_min', 'PET_min', 'AET_min', 'NDVI_min', 'alt_range', 'seasonality')
  }
  for (pred in pred_names){
    pred_file = paste(pred_folder, pred, '.tif', sep = '')
    pred_raster = c(as.matrix(raster(pred_file)))
    dat = cbind(dat, pred_raster)
  }
  dat = as.data.frame(dat)
  dat[, 6] = dat[, 4] * dat[, 5]
  names(dat) = c('S', pred_names, 'altbyseas')
  # Only keep rows where all variables are available
  index = complete.cases(dat)
  dat = dat[index, ]
  
  # Optional transformation of S
  if (transform == 'log'){
    dat$S = log(dat$S)
  }
  else if (transform == 'sqrt'){
    dat$S = sqrt(dat$S)
  }
  
  out = vector(mode = 'numeric', length = 13)
  full_model = lm(dat$S ~ dat$AET + dat$NDVI + dat$alt_range*dat$seasonality)
  out[1] = summary(full_model)$r.squared
  out[9:13] = as.numeric(summary(full_model)$coef[2:6, 4])
  prod_model = lm(dat$S ~ dat$AET + dat$NDVI)
  out[2] = summary(prod_model)$r.squared
  Janzen_model = lm(dat$S ~ dat$alt_range*dat$seasonality)
  out[3] = summary(Janzen_model)$r.squared
  out[4] = out[1] - out[3] # Unique r2 of productivity
  out[5] = out[1] - out[2] # Unique r2 of Janzen
  no_alt_model = lm(dat$S ~ dat$seasonality + dat$alt_range:dat$seasonality)
  no_seas_model = lm(dat$S ~ dat$alt_range + dat$alt_range:dat$seasonality)
  no_int_model = lm(dat$S ~ dat$seasonality + dat$alt_range)
  out[6] = out[3] - summary(no_alt_model)$r.squared
  out[7] = out[3] - summary(no_seas_model)$r.squared
  out[8] = out[3] - summary(no_int_model)$r.squared
  return(out)
}

# The function quanreg_constraint() applies quantile regression at quantile = 0.9
# to transformed or untransformed S(q) with respect to each of a given list of 
# predictor variables. In each grid cell, it records the lowest predicted S(q) 
# among the regressions, as well as the predictor associated with this minimum.
# The output is returned as a list with the two columns.
quanreg_constraint = function(S_dir, pred_vars, transform = 'untransformed'){
  dat = c(as.matrix(raster(S_dir)))
  for (var in pred_vars){
    var_dir = paste(pred_folder, var, '.tif', sep = '')
    dat = cbind(dat, c(as.matrix(raster(var_dir))))
  }
  dat = as.data.frame(dat)
  names(dat) = c('S', pred_vars)
  index = complete.cases(dat)
  
  S_keep = dat$S[index]
  for (i in 2:(length(pred_vars) + 1)){
    pred_var = dat[index, i]
    if (transform == 'untransformed'){
      model = rq(S_keep ~ pred_var, tau = 0.9)
      dat[index, i + length(pred_vars)] = fitted(model)
    }
    else if (transform == 'sqrt'){
      model = rq(sqrt(S_keep) ~ pred_var, tau = 0.9)
      dat[index, i + length(pred_vars)] = fitted(model)^2
    }
    else if (transform == 'log'){
      model = rq(log(S_keep) ~ pred_var, tau = 0.9)
      dat[index, i + length(pred_vars)] = exp(fitted(model))
    }
  }
  dat$pred_min = apply(dat[, (length(pred_vars) + 2):(2 * length(pred_vars) + 1)], 1, min)
  dat$limit_var = apply(dat[, (length(pred_vars) + 2):(2 * length(pred_vars) + 1)], 1, 
                        function(x) ifelse(anyNA(x), NA, pred_vars[which(x == min(x))]))
  return(list(dat$pred_min, dat$limit_var))
}
