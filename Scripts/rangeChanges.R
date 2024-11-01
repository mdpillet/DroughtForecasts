library(lmerTest)
library(ggplot2)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
outputPath <- "Outputs/"
rangeSizePath <- "/RangeSizes/"
figPath <- "Manuscript/Figures/"
tablePath <- "Manuscript/Tables/"

# List species
speciesDirs <- list.dirs(paste0(relPath, outputPath), recursive = F)

# Number of range change calculations per species
offset <- 36

# Create data frame for range changes
rangeChanges <- data.frame(Species = character(offset * length(speciesDirs)),
                           VariableSet = character(offset * length(speciesDirs)),
                           DroughtGCM = character(offset * length(speciesDirs)),
                           VariableSetDroughtGCM = character(offset * length(speciesDirs)),
                           FutureGCM = character(offset * length(speciesDirs)),
                           SSP = character(offset * length(speciesDirs)),
                           RangeChange = numeric(offset * length(speciesDirs)))

# Loop over species and calculate range changes
counter <- 0
for (i in 1:length(speciesDirs)) {
  modelFiles <- list.files(paste0(speciesDirs[i], rangeSizePath))
  sp <- strsplit(speciesDirs[i], "/", fixed = T)[[1]]
  sp <- sp[length(sp)]
  print(sp)
  for (j in 1:length(modelFiles)) {
    variableSet <- toupper(strsplit(modelFiles[j], "_")[[1]][1])
    modelGCM <- strsplit(modelFiles[j], "_")[[1]][2]
    tmp <- read.csv(paste0(speciesDirs[i], rangeSizePath, modelFiles[j]), header = T)
    for (k in 2:nrow(tmp)) {
      rangeChange <- tmp[k, "RangeSize"] / tmp[1, "RangeSize"]
      forecastGCM <- strsplit(tmp[k, "Map"], "_")[[1]][2]
      SSP <- strsplit(tmp[k, "Map"], "_")[[1]][3]
      counter <- counter + 1
      rangeChanges[counter, "Species"] <- sp
      rangeChanges[counter, "VariableSet"] <- variableSet
      rangeChanges[counter, "DroughtGCM"] <- modelGCM
      rangeChanges[counter, "VariableSetDroughtGCM"] <- paste0(variableSet, "_", modelGCM)
      rangeChanges[counter, "FutureGCM"] <- forecastGCM
      rangeChanges[counter, "SSP"] <- SSP
      rangeChanges[counter, "RangeChange"] <- rangeChange
    }
  }
}
rangeChanges[rangeChanges$DroughtGCM == "current", "DroughtGCM"] <- "NA"

# Regress range changes
# lmodel <- lm(RangeChange ~ VariableSet + SSP + DroughtGCM + FutureGCM, data = rangeChanges)
lmodel <- lm(RangeChange ~ VariableSetDroughtGCM + SSP + FutureGCM, data = rangeChanges)
summary(lmodel)
confint(lmodel)
# lmmodel <- lmer(RangeChange ~ VariableSet + SSP + DroughtGCM + FutureGCM + (1|Species), data = rangeChanges)
lmmodel <- lmer(RangeChange ~ VariableSetDroughtGCM + SSP + FutureGCM + (1|Species), data = rangeChanges)
coef(summary(lmmodel))
confint(lmmodel)

# Boxplots of range size changes
rangeChanges[rangeChanges$VariableSet == "AE", "VariableSet"] <- "AD"
rangeChanges[rangeChanges$VariableSet == "ASE", "VariableSet"] <- "ASD"
rangeChanges[rangeChanges$SSP == "SSP126", "SSP"] <- "SSP1-2.6"
rangeChanges[rangeChanges$SSP == "SSP370", "SSP"] <- "SSP3-7.0"
rangeChanges[rangeChanges$SSP == "SSP585", "SSP"] <- "SSP5-8.5"
png(paste0(relPath, figPath, "Fig2.png"), units = "mm", width = 180, height = 180, res = 300)
ggplot(data = rangeChanges, aes(x = VariableSet, y = RangeChange)) +
  geom_boxplot(outlier.size = 0.5) + ylim(0, 2) + facet_grid(FutureGCM ~ SSP) + geom_hline(yintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Variable set") + ylab("Predicted change in suitable climatic area")
dev.off()

# Export range change table for Supplementary Information
rangeChangesTable <- rangeChanges[, c("Species", "VariableSetDroughtGCM", "FutureGCM", "SSP", "RangeChange")]
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "A_current", "VariableSetDroughtGCM"] <- "A"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "AS_current", "VariableSetDroughtGCM"] <- "AS"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "AE_GFDL-ESM4", "VariableSetDroughtGCM"] <- "AD (GFDL-ESM4)"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "AE_MPI-ESM1-2-HR", "VariableSetDroughtGCM"] <- "AD (MPI-ESM1-2-HR)"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "AE_UKESM1-0-LL", "VariableSetDroughtGCM"] <- "AD (UKESM1-0-LL)"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "ASE_GFDL-ESM4", "VariableSetDroughtGCM"] <- "ASD (GFDL-ESM4)"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "ASE_MPI-ESM1-2-HR", "VariableSetDroughtGCM"] <- "ASD (MPI-ESM1-2-HR)"
rangeChangesTable[rangeChangesTable$VariableSetDroughtGCM == "ASE_UKESM1-0-LL", "VariableSetDroughtGCM"] <- "ASD (UKESM1-0-LL)"
write.csv(rangeChangesTable, paste0(relPath, tablePath, "SupplementaryTable2.csv"), row.names = F)

# Summarize directional range changes by variable set, future GCM, and SSP
rangeChanges$ChangeDirection <- ifelse(rangeChanges$RangeChange >= 1, "+", "-")
table(subset(rangeChanges, VariableSet == "A")$ChangeDirection, subset(rangeChanges, VariableSet == "A")$SSP, subset(rangeChanges, VariableSet == "A")$FutureGCM) / length(speciesDirs)
table(subset(rangeChanges, VariableSet == "A")$ChangeDirection, subset(rangeChanges, VariableSet == "A")$SSP, subset(rangeChanges, VariableSet == "A")$FutureGCM)
table(subset(rangeChanges, VariableSet == "AS")$ChangeDirection, subset(rangeChanges, VariableSet == "AS")$SSP, subset(rangeChanges, VariableSet == "AS")$FutureGCM) / length(speciesDirs)
table(subset(rangeChanges, VariableSet == "AS")$ChangeDirection, subset(rangeChanges, VariableSet == "AS")$SSP, subset(rangeChanges, VariableSet == "AS")$FutureGCM)
table(subset(rangeChanges, VariableSet == "AD")$ChangeDirection, subset(rangeChanges, VariableSet == "AD")$SSP, subset(rangeChanges, VariableSet == "AD")$FutureGCM) / length(speciesDirs)
table(subset(rangeChanges, VariableSet == "AD")$ChangeDirection, subset(rangeChanges, VariableSet == "AD")$SSP, subset(rangeChanges, VariableSet == "AD")$FutureGCM)
table(subset(rangeChanges, VariableSet == "ASD")$ChangeDirection, subset(rangeChanges, VariableSet == "ASD")$SSP, subset(rangeChanges, VariableSet == "ASD")$FutureGCM) / length(speciesDirs)
table(subset(rangeChanges, VariableSet == "ASD")$ChangeDirection, subset(rangeChanges, VariableSet == "ASD")$SSP, subset(rangeChanges, VariableSet == "ASD")$FutureGCM)