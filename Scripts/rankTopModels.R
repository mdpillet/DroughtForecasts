# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
outputPath <- "Outputs/"
modelSetPath <- "/ModelSetComparison/modelSetComparison.csv"

# List species
speciesDirs <- list.dirs(paste0(relPath, outputPath), recursive = F)

# Create ranking data frame
rankings <- data.frame(VariableSet = character(length(speciesDirs)),
                       GCM = character(length(speciesDirs)))

# Find top models for each species
for (i in 1:length(speciesDirs)) {
  comp <- read.csv(paste0(speciesDirs[i], modelSetPath), header = T)
  rankings[i, "VariableSet"] <- subset(comp, DeltaAICc == 0)$VariableSet
  rankings[i, "GCM"] <- subset(comp, DeltaAICc == 0)$GCM
}

# Rank models...
rankings$ModelNameFull <- paste0(rankings$VariableSet, "_", rankings$GCM)
# ...by variable set and GCM
table(rankings$ModelNameFull)
table(rankings$ModelNameFull) / length(speciesDirs) * 100
# ...by variable set only
table(rankings$VariableSet)
table(rankings$VariableSet) / length(speciesDirs) * 100
