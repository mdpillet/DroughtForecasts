library(terra)
library(spThin)

# Set directory structure and load workflow functions
spPath <- "E:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Filtered/Over10/"
thinPath <- "E:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Thinned/"
summaryFile <- "E:/Research/DroughtForecasts/Data/Occurrences/BySpecies/thinningSummary.csv"

# Set seed number
seedNo <- 2024

# Set parameters
thinRadius <- c(0, 5, 10) # Kilometers for spatial thinning

# List species to be modeled
spList <- list.files(spPath, "shp$", full.names = T)

# Thinning summary
thinSumm <- data.frame(Taxon = character(),
                       OccRaw = integer())
thinSumm <- cbind(thinSumm, as.data.frame(matrix(ncol = length(thinRadius), nrow = 0,
                                                 dimnames = list(NULL, paste0("OccThin_", thinRadius)))))

# Spatially thin occurrences
for (currentSpecies in 1:length(spList)) {
  # Get species name
  sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
  sp <- sp[length(sp)]
  sp <- strsplit(sp, ".", fixed = T)[[1]][1]
  thinSumm[currentSpecies, "Taxon"] <- sp
  
  # Read occurrences
  spOcc <- vect(spList[currentSpecies])
  thinSumm[currentSpecies, "OccRaw"] <- nrow(spOcc)
  
  # Thinning
  for (i in 1:length(thinRadius)) {
    if (thinRadius[i] != 0) {
      thinDF <- crds(spOcc, df = T)
      thinDF$sp <- sp
      thinnedOcc <- thin(thinDF, lat.col = "y", long.col = "x", spec.col = "sp", thin.par = thinRadius[i], reps = 1, 
                         locs.thinned.list.return = T, write.files = F, write.log.file = F, verbose = F)[[1]]
      spOccThinned <- spOcc[as.integer(row.names(thinnedOcc))]  
      thinSumm[currentSpecies, paste0("OccThin_", thinRadius[i])] <- nrow(spOccThinned)
      writeVector(spOccThinned, paste0(thinPath, sp, "_thin", thinRadius[i], ".shp"), overwrite = T)
    } else {
      thinSumm[currentSpecies, paste0("OccThin_", thinRadius[i])] <- nrow(spOcc)
      writeVector(spOcc, paste0(thinPath, sp, "_thin", thinRadius[i], ".shp"), overwrite = T)
    }
  }
  print(thinSumm[currentSpecies,])
}

# Export summary
write.csv(thinSumm, summaryFile, row.names = F)