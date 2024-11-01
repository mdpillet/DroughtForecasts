library(terra)

# Read masking raster
maskRaster <- rast("D:/Research/DroughtForecasts/Data/Predictors/ase_UKESM1-0-LL_SSP585.tif")
maskRaster <- maskRaster[[c(1:26, 28, 31, 34, 37)]] 
getMask <- anyNA(maskRaster)

# List all rasters
predFiles <- list.files("D:/Research/DroughtForecasts/Data/Predictors/", "tif$", full.names = T)

# Mask all rasters
for (i in 1:length(predFiles)) {
  tmp <- rast(predFiles[i])
  tmp <- mask(tmp, getMask, maskvalue = T)
  plot(tmp)
  writeRaster(tmp, predFiles[i], overwrite = T)
}