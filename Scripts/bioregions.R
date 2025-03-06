library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtPredictions/"
summaryPath <- "Data/DiversityMaps/"
changePath <- "Data/DiversityChangeMaps/"
bioregionPath <- "Data/Bioregions/"

# Read bioregions file
bioregions <- vect(paste0(relPath, bioregionPath, "bioregions_nonhierarchical.shp"))
bioregions <- aggregate(bioregions, by = "bioregio")

# Create change maps
d0_current <- rast(paste0(relPath, summaryPath, "0km_current.tif"))
d0_SSP126 <- rast(paste0(relPath, summaryPath, "0km_SSP126.tif"))
d0_SSP370 <- rast(paste0(relPath, summaryPath, "0km_SSP370.tif"))
d0_SSP585 <- rast(paste0(relPath, summaryPath, "0km_SSP585.tif"))
d100_current <- rast(paste0(relPath, summaryPath, "100km_current.tif"))
d100_SSP126 <- rast(paste0(relPath, summaryPath, "100km_SSP126.tif"))
d100_SSP370 <- rast(paste0(relPath, summaryPath, "100km_SSP370.tif"))
d100_SSP585 <- rast(paste0(relPath, summaryPath, "100km_SSP585.tif"))
writeRaster(d0_SSP126 - d0_current, paste0(relPath, changePath, "abs_d0_SSP126.tif"))
writeRaster(d0_SSP370 - d0_current, paste0(relPath, changePath, "abs_d0_SSP370.tif"))
writeRaster(d0_SSP585 - d0_current, paste0(relPath, changePath, "abs_d0_SSP585.tif"))
writeRaster(d100_SSP126 - d100_current, paste0(relPath, changePath, "abs_d100_SSP126.tif"))
writeRaster(d100_SSP370 - d100_current, paste0(relPath, changePath, "abs_d100_SSP370.tif"))
writeRaster(d100_SSP585 - d100_current, paste0(relPath, changePath, "abs_d100_SSP585.tif"))
writeRaster(d0_SSP126 / d0_current, paste0(relPath, changePath, "rel_d0_SSP126.tif"))
writeRaster(d0_SSP370 / d0_current, paste0(relPath, changePath, "rel_d0_SSP370.tif"))
writeRaster(d0_SSP585 / d0_current, paste0(relPath, changePath, "rel_d0_SSP585.tif"))
writeRaster(d100_SSP126 / d100_current, paste0(relPath, changePath, "rel_d100_SSP126.tif"))
writeRaster(d100_SSP370 / d100_current, paste0(relPath, changePath, "rel_d100_SSP370.tif"))
writeRaster(d100_SSP585 / d100_current, paste0(relPath, changePath, "rel_d100_SSP585.tif"))

# List change files
absChangeFiles <- list.files(paste0(relPath, changePath), "abs")
relChangeFiles <- list.files(paste0(relPath, changePath), "rel")

# Extract statistics by bioregion (mean for absolute changes per bioregion, median for relative changes per bioregion)
absRankings <- matrix(nrow = nrow(bioregions), ncol = length(absChangeFiles))
for (i in 1:length(absChangeFiles)) {
  print(i)
  tmp <- rast(paste0(relPath, changePath, absChangeFiles[i]))
  for (j in 1:nrow(bioregions)) {
    extr <- terra::extract(tmp, bioregions[j], ID = F)
    extrMean <- mean(extr$sum, na.rm = T)
    absRankings[j, i] <- extrMean
  }
}
absRankings <- as.data.frame(absRankings)
rownames(absRankings) <- bioregions$bioregio
colnames(absRankings) <- absChangeFiles
write.csv(absRankings, paste0(relPath, bioregionPath, "rankings_changeAbs.csv"))

relRankings <- matrix(nrow = nrow(bioregions), ncol = length(relChangeFiles))
for (i in 1:length(relChangeFiles)) {
  print(i)
  tmp <- rast(paste0(relPath, changePath, relChangeFiles[i]))
  for (j in 1:nrow(bioregions)) {
    extr <- terra::extract(tmp, bioregions[j], ID = F)
    extrMean <- median(extr$sum, na.rm = T)
    relRankings[j, i] <- extrMean
  }
}
relRankings <- as.data.frame(relRankings)
rownames(relRankings) <- bioregions$bioregio
colnames(relRankings) <- relChangeFiles
write.csv(relRankings, paste0(relPath, bioregionPath, "rankings_changeRel.csv"))