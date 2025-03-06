library(dplyr)
library(tidyr)
library(ggplot2)

# Set directory structure
relPath <- "D:/Research/DroughtPredictions/"
outputPath <- "Outputs/"
figPath <- "Manuscript/Figures/"
tablePath <- "Manuscript/Tables/"

# Summarize range change results
species <- list.dirs(paste0(relPath, outputPath), recursive = F)
results <- data.frame(species = species,
                      bestModel = character(length(species)),
                      AUC = numeric(length(species)),
                      variables = character(length(species)),
                      droughtsIncluded = logical(length(species))) 

for (i in 1:length(species)) {
  print(i)
  if (file.exists(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"))) {
    tmp <- read.csv(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"), header = T)
    bestModel <- subset(tmp, AverageValidationAUC > 0.5 & AUC_50000 > 0.7)
    if (nrow(bestModel) > 0) {
      changes <- read.csv(paste0(species[i], "/RangeSizesBestModel/rangeSizes.csv"))
      results[i, "bestModel"] <- changes[1, "GeneratingModel"]
      results[i, "changeSSP1-2.6"] <- changes[2, "RangeSizeChange"]
      results[i, "changeSSP3-7.0"] <- changes[3, "RangeSizeChange"]
      results[i, "changeSSP5-8.5"] <- changes[4, "RangeSizeChange"]
      changes <- read.csv(paste0(species[i], "/RangeSizesBestModel/rangeSizesNoDisp.csv"))
      results[i, "changeSSP1-2.6NoDisp"] <- changes[2, "RangeSizeChange"]
      results[i, "changeSSP3-7.0NoDisp"] <- changes[3, "RangeSizeChange"]
      results[i, "changeSSP5-8.5NoDisp"] <- changes[4, "RangeSizeChange"]
      results[i, "AUC"] <- bestModel[1, "AverageValidationAUC"]
      results[i, "AUC_0"] <- bestModel[1, "AUC_0"]
      results[i, "AUC_25"] <- bestModel[1, "AUC_25000"]
      results[i, "AUC_50"] <- bestModel[1, "AUC_50000"]
      results[i, "AUC_75"] <- bestModel[1, "AUC_75000"]
      results[i, "variables"] <- bestModel[1, "Variables"]
      if (grepl("H(A[DS]|M[SD])", results[i, "variables"])) {
        results[i, "droughtsIncluded"] <- T
      }
      else {
        results[i, "droughtsIncluded"] <- F
      }
    }
    else {
      results[i, "bestModel"] <- NA
    }
  }
  else {
    results[i, "bestModel"] <- NA
  }
}

# Plot variable inclusion rate
results <- subset(results, !is.na(bestModel))
results$BIO1 <- grepl("BIO1", results$variables)
results$BIO2 <- grepl("BIO2", results$variables)
results$BIO3 <- grepl("BIO3", results$variables)
results$BIO4 <- grepl("BIO4", results$variables)
results$BIO5 <- grepl("BIO5", results$variables)
results$BIO6 <- grepl("BIO6", results$variables)
results$BIO7 <- grepl("BIO7", results$variables)
results$BIO8 <- grepl("BIO8", results$variables)
results$BIO9 <- grepl("BIO9", results$variables)
results$BIO10 <- grepl("BIO10", results$variables)
results$BIO11 <- grepl("BIO11", results$variables)
results$BIO12 <- grepl("BIO12", results$variables)
results$BIO13 <- grepl("BIO13", results$variables)
results$BIO14 <- grepl("BIO14", results$variables)
results$BIO15 <- grepl("BIO15", results$variables)
results$BIO16 <- grepl("BIO16", results$variables)
results$BIO17 <- grepl("BIO17", results$variables)
results$BIO18 <- grepl("BIO18", results$variables)
results$BIO19 <- grepl("BIO19", results$variables)
results$pH <- grepl("pH", results$variables)
results$CEC <- grepl("CEC", results$variables)
results$Sand <- grepl("sand", results$variables)
results$Silt <- grepl("silt", results$variables)
results$Clay <- grepl("clay", results$variables)
results$CFVO <- grepl("CFVO", results$variables)
results$SOC <- grepl("SOC", results$variables)
results$HAD_mean <- grepl("HAD_mean", results$variables)
results$HAD_q0.5 <- grepl("HAD_q0.5", results$variables)
results$HAD_q0.75 <- grepl("HAD_q0.75", results$variables)
results$HAS_mean <- grepl("HAS_mean", results$variables)
results$HAS_q0.5 <- grepl("HAS_q0.5", results$variables)
results$HAS_q0.75 <- grepl("HAS_q0.75", results$variables)
results$HMD_mean <- grepl("HMD_mean", results$variables)
results$HMD_q0.5 <- grepl("HMD_q0.5", results$variables)
results$HMD_q0.75 <- grepl("HMD_q0.75", results$variables)
results$HMS_mean <- grepl("HMS_mean", results$variables)
results$HMS_q0.5 <- grepl("HMS_q0.5", results$variables)
results$HMS_q0.75 <- grepl("HMS_q0.75", results$variables)

df_long <- results[, 16:53] %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Included") %>%
  filter(Included)

df_fraction <- df_long %>%
  count(Variable) %>%
  mutate(Fraction = n / nrow(results),
         Color = case_when(grepl("BIO", Variable) ~ "Average", 
                           Variable == "pH" ~ "Soil",
                           TRUE ~ "SDE"))

custom_labels <- c("HAS_q0.75" = "HAS (75th percentile)",
                   "HMS_q0.75" = "HMS (75th percentile)",
                   "HAD_q0.75" = "HAD (75th percentile)",
                   "HMD_q0.75" = "HMD (75th percentile)",
                   "HAS_mean" = "HAS (mean)",
                   "HMS_mean" = "HMS (mean)",
                   "HAD_mean" = "HAD (mean)",
                   "HMD_mean" = "HMD (mean)")
png(paste0(relPath, figPath, "Fig1.png"), units = "mm", width = 180, height = 180, res = 300)
ggplot(df_fraction, aes(x = reorder(Variable, Fraction), y = Fraction, fill = Color)) +
  geom_bar(stat = "identity", col = "black") +
  coord_flip() +
  xlab("Variable") +
  ylab("Fraction of models") +
  theme_bw() +
  scale_fill_manual(values = c("black", "white", "grey")) + theme(legend.position = "none") +
  scale_x_discrete(labels = custom_labels)
dev.off()

# Compare model performance
results <- data.frame(species = species,
                      bestModel = character(length(species)),
                      AUC = numeric(length(species)),
                      variables = character(length(species)),
                      droughtsIncluded = logical(length(species)))
for (i in 1:length(species)) {
  print(i)
  if (file.exists(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"))) {
    tmp <- read.csv(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"), header = T)
    bestModel <- subset(tmp, AverageValidationAUC > 0.5 & AUC_50000 > 0.7)
    if (nrow(bestModel) > 0) {
      results[i, "bestModel"] <- paste0(tolower(bestModel[1, "VariableSet"]), "_", bestModel[1, "GCM"], "_current.csv")
      results[i, "AUC"] <- bestModel[1, "AverageValidationAUC"]
      results[i, "AUC_0"] <- bestModel[1, "AUC_0"]
      results[i, "AUC_25"] <- bestModel[1, "AUC_25000"]
      results[i, "AUC_50"] <- bestModel[1, "AUC_50000"]
      results[i, "AUC_75"] <- bestModel[1, "AUC_75000"]
      results[i, "variables"] <- bestModel[1, "Variables"]
    }
    else {
      results[i, "bestModel"] <- NA
      results[i, "change"] <- NA
    }
  }
  else {
    results[i, "bestModel"] <- NA
    results[i, "change"] <- NA
  }
}
results <- subset(results, !is.na(bestModel))
resultsAUCPlot <- results[, c("AUC", "AUC_0", "AUC_25", "AUC_50", "AUC_75")]
names(resultsAUCPlot) <- c("Cross-validated", "Full model (unfiltered)", "Full model (25 km filter)", "Full model (50 km filter)", "Full model (75 km filter)")

df_long <- resultsAUCPlot %>%
  pivot_longer(cols = everything(), names_to = "AUC_Type", values_to = "Value")
df_long$AUC_Type <- factor(df_long$AUC_Type, levels = c("Cross-validated", "Full model (unfiltered)", "Full model (25 km filter)", "Full model (50 km filter)", "Full model (75 km filter)"))

# Create side-by-side boxplots
png(paste0(relPath, figPath, "ExtendedDataFig1.png"), units = "mm", width = 270, height = 180, res = 300)
ggplot(df_long, aes(x = AUC_Type, y = Value)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "AUC type", y = "AUC") +
  theme(legend.position = "none")
dev.off()

# Create model quality table
results <- data.frame(species = species,
                      bestModel = character(length(species)),
                      AUC = numeric(length(species)))
counter <- 1
for (i in 1:length(species)) {
  print(i)
  if (file.exists(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"))) {
    tmp <- read.csv(paste0(species[i], "/ModelSetComparison/modelSetComparison.csv"), header = T)
    bestModel <- tmp
    if (nrow(bestModel) > 0) {
      for (j in 1:nrow(bestModel)) {
        results[counter, "species"] <- species[i]
        results[counter, "bestModel"] <- bestModel[j, "GCM"]
        results[counter, "AUC"] <- bestModel[j, "AverageValidationAUC"]
        results[counter, "AUC_0"] <- bestModel[j, "AUC_0"]
        results[counter, "AUC_25"] <- bestModel[j, "AUC_25000"]
        results[counter, "AUC_50"] <- bestModel[j, "AUC_50000"]
        results[counter, "AUC_75"] <- bestModel[j, "AUC_75000"]        
        results[counter, "BIC"] <- bestModel[j, "BIC"]
        results[counter, "DeltaBIC"] <- bestModel[j, "DeltaBIC"]
        results[counter, "LambdaRule"] <- bestModel[j, "LambdaRule"]
        results[counter, "AverageOmissionRate"] <- bestModel[j, "AverageOmissionRate"]
        results[counter, "Variables"] <- bestModel[j, "Variables"]
        counter <- counter + 1
      }
    }
  }
}

results$species <- gsub(paste0(relPath, outputPath), "", results$species, fixed = T)
write.csv(results, paste0(relPath, tablePath, "SupplementaryTable2.csv"), row.names = F)