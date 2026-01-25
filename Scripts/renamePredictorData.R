# List predictor files
tmp <- list.files("E:/Research/DroughtForecasts/Data/Predictors/", "tif$", full.names = T)

# Rename variables
for (i in 1:length(tmp)) {
  tmpRast <- rast(tmp[i])
  if (grepl("ase_", tmp[i], fixed = T)) {
    varNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC",
                                                      "HAD_mean", "HAD_median", "HAD_q0.75",
                                                      "HAS_mean", "HAS_median", "HAS_q0.75",
                                                      "HMD_mean", "HMD_median", "HMD_q0.75",
                                                      "HMS_mean", "HMS_median", "HMS_q0.75")
  } else if (grepl("as_", tmp[i], fixed = T)) {
    varNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC")
  }
  names(tmpRast) <- varNames
  terra::writeRaster(tmpRast, gsub(".tif", "_renamed.tif", tmp[i], fixed = T), overwrite = T)
}

# Delete original files
file.remove(tmp)

# Rename files back to original names
tmp <- list.files("E:/Research/DroughtForecasts/Data/Predictors/", "tif$", full.names = T)
for (i in 1:length(tmp)) {
  outName <- gsub("_renamed", "", tmp[i], fixed = T)
  file.rename(tmp[i], outName)
}