library(terra)

# List all rasters
predFiles <- list.files("D:/Research/DroughtForecasts/Data/Predictors/", "tif$", full.names = T)

# Subset rasters to discard severe drought variables
for (i in 1:length(predFiles)) {
  tmp <- rast(predFiles[i])
  tmp <- tmp[[1:26]]
  outName <- gsub("ase_", "as_", predFiles[i], fixed = T)
  if (grepl("current", outName)) {
    outName <- gsub("GFDL-ESM4|MPI-ESM1-2-HR|UKESM1-0-LL", "NA", outName)
  }
  writeRaster(tmp, outName, overwrite = T)
}