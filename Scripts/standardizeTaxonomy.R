library(TNRS)
library(dplyr)

# Set paths
relDir <- "D:/Research/DroughtForecasts/Data/"
inPath <- "Occurrences/combined_coordsCleaned.csv"
inPathWorking <- "TNRS/namesMatched_working.csv"
outPath <- "TNRS/namesMatched_raw.csv"
outPathFinal <- "TNRS/namesMatched_final.csv"
occPathFinal <- "Occurrences/combined_final.csv"

# Read occurrence data
occ <- read.csv(paste0(relDir, inPath), header = T)

# Get unique names
uniqNames <- unique(occ$Taxon)

# Check taxonomy
occNames <- TNRS(taxonomic_names = uniqNames,
     sources = "cact",
     url = "http://vegbiendev.nceas.ucsb.edu:9975/tnrs_api.php")

# Export names
write.csv(occNames, paste0(relDir, outPath), row.names = F)

# Import names
taxonNames <- read.csv(paste0(relDir, inPathWorking), header = T)

# Discard non-resolved names
taxonNames <- subset(taxonNames, !is.na(Name_infraspecific_matched_override))
# Create final names
for (i in 1:nrow(taxonNames)) {
  if (taxonNames[i, "Name_infraspecific_matched_override"] == "" | taxonNames[i, "Name_infraspecific_matched_override"] == "Accept") {
    taxonNames[i, "Final_name_infraspecific"] <- taxonNames[i, "Accepted_name"]
    taxonNames[i, "Final_name_species"] <- taxonNames[i, "Accepted_species"]
  }
  else {
    taxonNames[i, "Final_name_infraspecific"] <- taxonNames[i, "Name_infraspecific_matched_override"]
    taxonNames[i, "Final_name_species"] <- taxonNames[i, "Name_species_matched_override"]
  }
}

# Export names
write.csv(taxonNames, paste0(relDir, outPathFinal), row.names = F)

# Update taxonomy for occurrence data
for (i in 1:nrow(taxonNames)) {
  print(i)
  tmpInputName <- taxonNames[i, "Name_submitted"]
  tmpOutputNameInfra <- taxonNames[i, "Final_name_infraspecific"]
  tmpOutputNameSpecies <- taxonNames[i, "Final_name_species"]
  occ[occ$Taxon == tmpInputName, "FinalSpecies"] <- tmpOutputNameSpecies
  occ[occ$Taxon == tmpInputName, "FinalInfraspecific"] <- tmpOutputNameInfra
}

# Only retain occurrences with resolved name
occFinal <- subset(occ, !is.na(FinalSpecies))

# Get summary statistics
length(unique(occFinal$FinalSpecies))
occBySpecies <- occFinal %>%
  group_by(FinalSpecies) %>%
  summarize(OccCount = n())
species10 <- subset(occBySpecies, OccCount >= 10)
write.csv(occBySpecies, paste0(relDir, "Occurrences/occCountAll.csv"), row.names = F)
write.csv(species10, paste0(relDir, "Occurrences/occCount10.csv"), row.names = F)

# Export final occurrences
write.csv(occFinal, paste0(relDir, occPathFinal), row.names = F)