library(rnaturalearth)
library(rnaturalearthdata)
library(GNRS)
library(NSR)
library(dplyr)

# Set paths
relDir <- "D:/Research/DroughtForecasts/Data/"
occPath <- "Occurrences/combined_final.csv"
nsrPath <- "NSR/checklistRaw.csv"
checklistPath <- "NSR/checklistFinal.csv"
nativePath <- "Occurrences/combined_final_nativesOnly.csv"

# Read occurrence data
occ <- read.csv(paste0(relDir, occPath), header = T)

# Get country information
occCoords <- data.frame(lon = occ$lon, lat = occ$lat)
occCoords$country <- extract(vect(countries50), occCoords)$name
rownames(occCoords) <- NULL

# Standardize country information
# Matches should be checked manually at this point!
query <- GNRS_template(nrow = length(unique(occCoords$country)))
query$country <- unique(occCoords$country)
results <- GNRS(query)
query$matchedCountry <- results$country

for (i in 1:nrow(query)) {
  print(query[i, "country"])
  occCoords[occCoords$country == query[i, "country"], "country"] <- query[i, "matchedCountry"] 
}

# Check for invasive species
occ$Country <- occCoords$country
counts10 <- occ %>% group_by(FinalSpecies) %>% mutate(count = n())
counts10 <- subset(counts10, count >= 10)
nsr <- NSR_template(nrow = nrow(counts10))
nsr$species <- counts10$FinalSpecies
nsr$country <- counts10$Country
queryNSR <- distinct(nsr)
resultsNSR <- NSR(queryNSR)
resultsNSR <- resultsNSR[, c("species", "country", "native_status")]
write.csv(resultsNSR, paste0(relDir, nsrPath), row.names = F)

# Remove records thought to be invasive, belonging to species with highly uncertain native range, and Rhipsalis baccifera
invasiveRecords <- read.csv(paste0(relDir, checklistPath), header = T)
invasiveRecords <- subset(invasiveRecords, native_status_override == "C" | native_status_override == "I")
nativeRecords <- occ
for (i in 1:nrow(invasiveRecords)) {
  print(i)
  nativeRecords <- nativeRecords[!(nativeRecords$FinalSpecies == invasiveRecords[i, "species"] &
                  nativeRecords$Country == invasiveRecords[i, "country"]),]
}
write.csv(nativeRecords[, 1:6], paste0(relDir, nativePath), row.names = F)