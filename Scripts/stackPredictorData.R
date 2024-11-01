library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
CHELSAPath <- "Data/CHELSA/Processed/Resampled/"
soilPath <- "Data/SoilGrids/Processed/"
BIOFI2Path <- "Data/BIOFI2/Processed/Resampled/"
outPath <- "Data/Predictors/"

# Stack SoilGrids
soilLayers <- rast(list.files(path = paste0(relPath, soilPath), pattern = "tif", full.names = T))

# Stack current data (A)
fileNamesCHELSACurrent <- list.files(path = paste0(relPath, CHELSAPath, "Current/"), pattern = "tif", full.names = T)
CHELSAcurrent <- rast(fileNamesCHELSACurrent)
writeRaster(CHELSAcurrent, paste0(relPath, outPath, "a_current_current.tif"), overwrite = T)

# Stack current data (AS)
writeRaster(c(CHELSAcurrent, soilLayers), paste0(relPath, outPath, "as_current_current.tif"), overwrite = T)

# Stack current data (AE and ASE)
parametersBIOFI2 <- c("GFDL-ESM4", "MPI-ESM1-2-HR", "UKESM1-0-LL")
fileNamesBIOFI2Current <- list.files(path = paste0(relPath, BIOFI2Path, "Current/"), pattern = "tif", full.names = T)
for (i in 1:length(parametersBIOFI2)) { 
  print(i)
  tmpStack <- rast(fileNamesBIOFI2Current[grepl(parametersBIOFI2[i], fileNamesBIOFI2Current)])
  writeRaster(c(CHELSAcurrent, tmpStack), paste0(relPath, outPath, "ae_", parametersBIOFI2[i], "_current.tif"), overwrite = T)
  writeRaster(c(CHELSAcurrent, soilLayers, tmpStack), paste0(relPath, outPath, "ase_", parametersBIOFI2[i], "_current.tif"), overwrite = T)
}

# Stack future data (A and AS)
fileNamesCHELSAFuture <- list.files(path = paste0(relPath, CHELSAPath, "Future/"), pattern = "tif", full.names = T)
parametersCHELSA <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll")
scenariosCHELSA <- c("ssp126", "ssp370", "ssp585")
for (i in 1:length(parametersCHELSA)) { 
  print(i)
  for (j in 1:length(scenariosCHELSA)) {
    tmpFiles <- grep(parametersCHELSA[i], fileNamesCHELSAFuture, value = T)
    tmpFiles <- grep(scenariosCHELSA[j], tmpFiles, value = T)
    tmpStack <- rast(tmpFiles)
    writeRaster(tmpStack, paste0(relPath, outPath, "a_", toupper(parametersCHELSA[i]), "_", toupper(scenariosCHELSA[j]), ".tif"), overwrite = T)  
    writeRaster(c(tmpStack, soilLayers), paste0(relPath, outPath, "as_", toupper(parametersCHELSA[i]), "_", toupper(scenariosCHELSA[j]), ".tif"), overwrite = T)  
  }
}

# Stack future data (AE and ASE)
fileNamesBIOFI2Future <- list.files(path = paste0(relPath, BIOFI2Path, "Future/"), pattern = "tif", full.names = T)
for (i in 1:length(parametersCHELSA)) { 
  print(i)
  for (j in 1:length(scenariosCHELSA)) {
    tmpFilesCHELSA <- grep(parametersCHELSA[i], fileNamesCHELSAFuture, value = T)
    tmpFilesCHELSA <- grep(scenariosCHELSA[j], tmpFilesCHELSA, value = T)
    tmpFilesBIOFI2 <- grep(toupper(parametersCHELSA[i]), fileNamesBIOFI2Future, value = T)
    tmpFilesBIOFI2 <- grep(scenariosCHELSA[j], tmpFilesBIOFI2, value = T)
    writeRaster(c(rast(tmpFilesCHELSA), rast(tmpFilesBIOFI2)), paste0(relPath, outPath, "ae_", toupper(parametersCHELSA[i]), "_", toupper(scenariosCHELSA[j]), ".tif"), overwrite = T)  
    writeRaster(c(rast(tmpFilesCHELSA), soilLayers, rast(tmpFilesBIOFI2)), paste0(relPath, outPath, "ase_", toupper(parametersCHELSA[i]), "_", toupper(scenariosCHELSA[j]), ".tif"), overwrite = T)  
  }
}