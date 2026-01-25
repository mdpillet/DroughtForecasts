library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
summaryPath <- "Data/DiversityMaps/"
changePath <- "Data/DiversityChangeMaps/"
bioregionPath <- "Data/Bioregions/"

# Read bioregions file
bioregions <- vect(paste0(relPath, bioregionPath, "bioregions_nonhierarchical.shp"))
bioregions <- aggregate(bioregions, by = "bioregio")

# Create change maps
current <- rast(paste0(relPath, summaryPath, "Current/avg.tif"))
future <- rast(paste0(relPath, summaryPath, "Future/avg.tif"))
writeRaster(future - current, paste0(relPath, changePath, "abs_avg.tif"))
writeRaster(future / current, paste0(relPath, changePath, "rel_avg.tif"))

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