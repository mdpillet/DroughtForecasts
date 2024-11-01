library(RCurl)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
outPath <- "Data/SoilGrids/Raw/"

# Download files
options(timeout = 36000)
urlPrefix <- "https://files.isric.org/soilgrids/latest/data_aggregated/1000m/"
vars <- c("phh2o", "cec", "clay", "sand", "silt", "cfvo", "soc")
depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm")
for (i in 1:length(vars)) {
  print(i)
  for (j in 1:length(depths)) {
    inName <- paste0(urlPrefix, vars[i], "/", vars[i], "_", depths[j], "_mean_1000.tif")
    outName <- paste0(vars[i], "_", depths[j], "_mean_1000.tif")
    if (!file.exists(paste0(relPath, outPath, outName))) {
      download.file(inName, paste0(relPath, outPath, outName), mode = "wb")  
    }
  }
}