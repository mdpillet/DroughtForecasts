library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
inPathSoilGrids <- "Data/SoilGrids/Processed/"
inPathCHELSA <- "D:/Research/DroughtForecasts/Data/BIOFI2/Raw/"
outPath <- "Data/BIOFI2/Processed/Resampled/"

# Read reference grid
refGrid <- rast(list.files(path = paste0(relPath, inPathSoilGrids), pattern = "cfvo", full.names = T))

# List files
climFiles <- list.files(inPathCHELSA, pattern = "tif$", recursive = T)

# Resample layers
scenarios <- c("historical", "ssp126", "ssp370", "ssp585")
GCMs <- c("GFDL-ESM4", "MPI-ESM1-2-HR", "UKESM1-0-LL")
chars <- c("hms12", "had12", "has12", "hmd12")

for (i in GCMs) { 
  print(i)
  tempGCM <- grep(i, climFiles, value = T)
  for (j in scenarios) {
    if (j == "historical") {
      tempYears <- grep("198[1-9]|199[0-9]|200[0-9]|2010", tempGCM, value = T)
    }
    else {
      tempYears <- grep("207[1-9]|20[8-9][0-9]|2100", tempGCM, value = T)
    }
    tempScenario <- grep(j, tempYears, value = T)
    for (k in chars) {
      tempChars <- grep(k, tempScenario, value = T)
      tempStack <- rast(paste0(inPathCHELSA, tempChars))
      layerAvg <- mean(tempStack, na.rm = T)
      layerQuant <- quantile(tempStack, probs = c(0.5, 0.75), na.rm = T)
      layerSumm <- c(layerAvg, layerQuant)
      names(layerSumm) <- paste0(k, "_", names(layerSumm))
      rawRast <- project(layerSumm, refGrid)
      resampledRast <- resample(rawRast, refGrid, method = "bilinear", threads = T)
      resampledRast <- crop(resampledRast, refGrid)
      if (j == "historical") {
        outName <- paste0(relPath, outPath, "Current/", k, "_", i, ".tif")
      }
      else {
        outName <- paste0(relPath, outPath, "Future/", k, "_", i, "_", j, ".tif")
      }
      writeRaster(resampledRast, outName, overwrite = T)
    }
  }
}