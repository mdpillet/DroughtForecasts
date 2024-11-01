# Set directory structure
relPath <- "D:/Research/Chapter3/"
inPath <- "Data/iNaturalist/Raw/"
outPath <- "Data/iNaturalist/Processed/observationsCombined.csv"

# Get list of observation files
obsFiles <- list.files(paste0(relPath, inPath), "*observations*")

# Combine observation files
obs <- read.csv(paste0(relPath, inPath, obsFiles[1]), header = T)
for (i in 2:length(obsFiles)) {
  tmp <- read.csv(paste0(relPath, inPath, obsFiles[i]), header = T)
  obs <- rbind(obs, tmp)
}

# Add in records for species flagged as geoprivate by iNaturalist
obscured <- read.csv(paste0(relPath, inPath, "TaxaGeoprivate_10-24-2023.csv"), header = T)
for (i in 1:nrow(obscured)) {
  tmpID <- obscured[i, "obs_id"]
  obs[obs$id == tmpID, "positional_accuracy"] <- obscured[i, "positional_accuracy"]
  obs[obs$id == tmpID, "latitude"] <- obscured[i, "latitude"]
  obs[obs$id == tmpID, "longitude"] <- obscured[i, "longitude"]
}

# Export data
write.csv(obs, paste0(relPath, outPath), row.names = F)