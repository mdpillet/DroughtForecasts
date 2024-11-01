library(terra)
library(XML)

# Set paths
relDir <- "D:/Research/DroughtForecasts/Data/"
occPath <- "Occurrences/combined_final_nativesOnly.csv"
outPath <- "Occurrences/BySpecies/All/"
out10Path <- "Occurrences/BySpecies/Over10/"

# Read occurrence data
occ <- read.csv(paste0(relDir, occPath), header = T)

# Create shapefile
spOcc <- vect(occ, crs = "+proj=longlat")

# Split by species
uniqSp <- unique(spOcc$FinalSpecies)
for (i in 1:length(uniqSp)) {
  print(i/length(uniqSp)*100)
  tmp <- spOcc[spOcc$FinalSpecies == uniqSp[i],]
  writeVector(tmp, paste0(relDir, outPath, uniqSp[i], ".kml"), overwrite = T)
  writeVector(tmp, paste0(relDir, outPath, uniqSp[i], ".shp"), overwrite = T)
  if (nrow(tmp) >= 10) {
    writeVector(tmp, paste0(relDir, out10Path, uniqSp[i], ".kml"), overwrite = T)
    writeVector(tmp, paste0(relDir, out10Path, uniqSp[i], ".shp"), overwrite = T)
  }
}