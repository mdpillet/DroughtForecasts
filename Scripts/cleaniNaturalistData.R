# Set directory structure
relPath <- "D:/Research/Chapter3/"
inPath <- "Data/iNaturalist/Processed/observationsCombined.csv"
outPath <- "Data/iNaturalist/Processed/observationsCombined_working.csv"
outPathFinal <- "Data/iNaturalist/Processed/INAT.csv"

# Read observations
obs <- read.csv(paste0(relPath, inPath), header = T)

# Remove extraneous columns
obsOut <- obs[, c("id", "observed_on", "quality_grade", "description", "captive_cultivated", 
                  "latitude", "longitude", "positional_accuracy", "private_latitude", "private_longitude", "public_positional_accuracy",
                  "geoprivacy", "taxon_geoprivacy", "coordinates_obscured", "scientific_name", "taxon_genus_name",
                  "taxon_genushybrid_name", "taxon_species_name", "taxon_hybrid_name", "taxon_subspecies_name",
                  "taxon_variety_name", "taxon_form_name")]

# Filter observations
obsOut <- subset(obsOut, quality_grade == "research")
obsOut <- subset(obsOut, captive_cultivated == "false")
obsOut <- subset(obsOut, !is.na(latitude) & !is.na(longitude))
obsOut <- subset(obsOut, geoprivacy != "obscured" & geoprivacy != "private")
obsOut <- subset(obsOut, !(taxon_genus_name == "" & taxon_genushybrid_name == ""))
obsOut <- subset(obsOut, !(taxon_species_name == "" & taxon_hybrid_name == ""))

# Fix aberrant observations
obsOut[obsOut$id == 67429082, "scientific_name"] <- "Eriosyce laui"

# Export observations
write.csv(obsOut, paste0(relPath, outPath), row.names = F)

# Read observations (manually cleaned: extraneous columns removed, hybrids without nothospecies removed, hybrid names formatted, occurrences before 1945 removed)
obs <- read.csv(paste0(relPath, outPath), header = T)

# Remove observations with low positional accuracy
obsOut <- subset(obs, positional_accuracy < 1000 | is.na(positional_accuracy))

# Add ID and remove extraneous columns
obsOut$AccNoGlobal <- paste0("INAT", obsOut$id)
obsOut <- obsOut[, c("observed_on", "description", "latitude", "longitude", "scientific_name", "AccNoGlobal")]

# Export observations
write.csv(obsOut, paste0(relPath, outPathFinal), row.names = F)