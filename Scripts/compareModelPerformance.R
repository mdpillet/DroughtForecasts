library(ggplot2)
library(tidyr)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
outputPath <- "Outputs/"
modelSetPath <- "/ModelSetComparison/modelSetComparison.csv"
figPath <- "Manuscript/Figures/"
tablePath <- "Manuscript/Tables/"

# List species
speciesDirs <- list.dirs(paste0(relPath, outputPath), recursive = F)

# Create model quality table for Supplementary Information
modelQuality <- lapply(paste0(speciesDirs, modelSetPath), read.csv, header = T)
for (i in 1:length(speciesDirs)) {
  modelQuality[[i]]$Species <- strsplit(speciesDirs[i], "/", fixed = T)[[1]][length(strsplit(speciesDirs[i], "/", fixed = T)[[1]])]
}
modelQuality <- do.call(rbind, modelQuality)
modelQuality <- modelQuality[, c("Species", "VariableSet", "GCM", "AICc", "DeltaAICc", "AverageValidationAUC", "AverageOmissionRate")]
modelQuality[modelQuality$VariableSet == "AE", "VariableSet"] <- "AD"
modelQuality[modelQuality$VariableSet == "ASE", "VariableSet"] <- "ASD"
write.csv(modelQuality, paste0(relPath, tablePath, "SupplementaryTable1.csv"), row.names = F)

# Create ranking data frame
rankings <- data.frame(A = numeric(length(speciesDirs)),
                       AS = numeric(length(speciesDirs)),
                       AE_GFDL = numeric(length(speciesDirs)),
                       AE_MPI = numeric(length(speciesDirs)),
                       AE_UKESM = numeric(length(speciesDirs)),
                       ASE_GFDL = numeric(length(speciesDirs)),
                       ASE_MPI = numeric(length(speciesDirs)),
                       ASE_UKESM = numeric(length(speciesDirs)))

# Get AUC for each species
for (i in 1:length(speciesDirs)) {
  comp <- read.csv(paste0(speciesDirs[i], modelSetPath), header = T)
  rankings[i, "A"] <- subset(comp, VariableSet == "A")$AverageValidationAUC
  rankings[i, "AS"] <- subset(comp, VariableSet == "AS")$AverageValidationAUC
  rankings[i, "AE_GFDL"] <- subset(comp, VariableSet == "AE" & GCM == "GFDL-ESM4")$AverageValidationAUC
  rankings[i, "AE_MPI"] <- subset(comp, VariableSet == "AE" & GCM == "MPI-ESM1-2-HR")$AverageValidationAUC
  rankings[i, "AE_UKESM"] <- subset(comp, VariableSet == "AE" & GCM == "UKESM1-0-LL")$AverageValidationAUC
  rankings[i, "ASE_GFDL"] <- subset(comp, VariableSet == "ASE" & GCM == "GFDL-ESM4")$AverageValidationAUC
  rankings[i, "ASE_MPI"] <- subset(comp, VariableSet == "ASE" & GCM == "MPI-ESM1-2-HR")$AverageValidationAUC
  rankings[i, "ASE_UKESM"] <- subset(comp, VariableSet == "ASE" & GCM == "UKESM1-0-LL")$AverageValidationAUC
}

# Create boxplot
rankingsTrans <- pivot_longer(rankings, cols = everything())
rankingsTrans$name <- as.factor(rankingsTrans$name)
levels(rankingsTrans$name) <- c("A", 
                                "AD \n(GFDL)",
                                "AD \n(MPI)",
                                "AD \n(UKESM1)",
                                "AS",
                                "ASD \n(GFDL)",
                                "ASD \n(MPI)",
                                "ASD \n(UKESM1)")

png(paste0(relPath, figPath, "Fig1.png"), units = "mm", width = 180, height = 120, res = 300)
ggplot(data = rankingsTrans, aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylab("AUROC") + ylim(0.5, 1) + xlab("Variable set")
dev.off()

# Get summary statistics
colMeans(rankings)
median(rankings$A)
median(rankings$AS)
median(rankings$ASE_GFDL)
median(rankings$ASE_MPI)
median(rankings$ASE_UKESM)
median(rankings$AE_GFDL)
median(rankings$AE_MPI)
median(rankings$AE_UKESM)

# Paired t-tests between groups
t.test(rankings$ASE_UKESM, rankings$A, paired = T, alternative = "greater")
t.test(rankings$ASE_GFDL, rankings$A, paired = T, alternative = "greater")
t.test(rankings$ASE_MPI, rankings$A, paired = T, alternative = "greater")
t.test(rankings$AE_UKESM, rankings$A, paired = T, alternative = "greater")
t.test(rankings$AE_GFDL, rankings$A, paired = T, alternative = "greater")
t.test(rankings$AE_MPI, rankings$A, paired = T, alternative = "greater")
t.test(rankings$AS, rankings$A, paired = T, alternative = "greater")
t.test(rankings$ASE_UKESM, rankings$AS, paired = T, alternative = "greater")
t.test(rankings$ASE_GFDL, rankings$AS, paired = T, alternative = "greater")
t.test(rankings$ASE_MPI, rankings$AS, paired = T, alternative = "greater")