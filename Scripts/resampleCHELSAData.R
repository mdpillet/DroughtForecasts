library(terra)

# Set directory structure
inPathSoilGrids <- "D:/Research/DroughtForecasts/Data/SoilGrids/Processed/"
inPathCHELSA <- "D:/Research/DroughtForecasts/Data/CHELSA/Raw/"
outPath <- "D:/Research/DroughtForecasts/Data/CHELSA/Processed/Resampled/"

# Read reference grid
refGrid <- rast(list.files(path = inPathSoilGrids, pattern = "cec", full.names = T))

# Resample layers
fileNames <- c(list.files(path = paste0(inPathCHELSA, "Current/"), pattern = "tif", full.names = T),
               list.files(path = paste0(inPathCHELSA, "Future/"), pattern = "tif", full.names = T))

for (i in 1:length(fileNames)) { 
  print(i)
  rawRast <- rast(fileNames[i])
  rawRast <- project(rawRast, refGrid)
  resampledRast <- resample(rawRast, refGrid, method = "bilinear", threads = T)
  resampledRast <- crop(resampledRast, refGrid)
  outName <- paste(strsplit(fileNames[i], "/", fixed = T)[[1]][7:8], collapse = "/")
  writeRaster(resampledRast, paste0(outPath, outName))
}