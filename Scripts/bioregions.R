library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
summaryPath <- "Data/DiversityMaps/Summary/"
bioregionPath <- "Data/Bioregions/"

# Read bioregions file
bioregions <- vect(paste0(relPath, bioregionPath, "bioregions_nonhierarchical.shp"))
bioregions <- aggregate(bioregions, by = "bioregio")

# List summary files
absChangeFiles <- list.files(paste0(relPath, summaryPath), "changeAbs.tif")
relChangeFiles <- list.files(paste0(relPath, summaryPath), "changeRel.tif")

# Extract statistics by bioregion (mean for absolute changes per bioregion, median for relative changes per bioregion)
absRankings <- matrix(nrow = nrow(bioregions), ncol = length(absChangeFiles))
for (i in 1:length(absChangeFiles)) {
  print(i)
  tmp <- rast(paste0(relPath, summaryPath, absChangeFiles[i]))
  for (j in 1:nrow(bioregions)) {
    extr <- extract(tmp, bioregions[j], ID = F)
    extrMean <- mean(extr$mean, na.rm = T)
    absRankings[j, i] <- extrMean
  }
}
absRankings <- as.data.frame(absRankings)
rownames(absRankings) <- bioregions$bioregio
colnames(absRankings) <- gsub("_changeAbs.tif", "", absChangeFiles)
write.csv(absRankings, paste0(relPath, bioregionPath, "rankings_changeAbs.csv"))

relRankings <- matrix(nrow = nrow(bioregions), ncol = length(relChangeFiles))
for (i in 1:length(relChangeFiles)) {
  print(i)
  tmp <- rast(paste0(relPath, summaryPath, relChangeFiles[i]))
  for (j in 1:nrow(bioregions)) {
    extr <- extract(tmp, bioregions[j], ID = F)
    extrMean <- median(extr$mean, na.rm = T)
    relRankings[j, i] <- extrMean
  }
}
relRankings <- as.data.frame(relRankings)
rownames(relRankings) <- bioregions$bioregio
colnames(relRankings) <- gsub("_changeRel.tif", "", relChangeFiles)
write.csv(relRankings, paste0(relPath, bioregionPath, "rankings_changeRel.csv"))