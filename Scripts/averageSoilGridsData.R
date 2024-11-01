library(terra)

# Set directory structure
inPath <- "D:/Research/DroughtForecasts/Data/SoilGrids/Raw/"
outPath <- "D:/Research/DroughtForecasts/Data/SoilGrids/Processed/"

# Average layers
vars <- c("phh2o", "cec", "clay", "sand", "silt", "cfvo", "soc")
extent <- ext(c(-14000000, -3500000, -61490000, 7000000))
for (i in 1:length(vars)) {
  print(i)
  fileNames <- list.files(path = inPath, pattern = vars[i], full.names = T)
  print(fileNames)
  tempStack <- rast(fileNames)
  tempStack <- crop(tempStack, extent)
  tempStackAvg <- tempStack[[1]] * 5/60 +
    tempStack[[2]] * 15/60 +
    tempStack[[3]] * 30/60 +
    tempStack[[4]] * 10/30
  writeRaster(tempStackAvg, paste0(outPath, vars[i], ".tif"))
}