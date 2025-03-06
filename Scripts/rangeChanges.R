library(lmerTest)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(broom.mixed)
library(dotwhisker)
library(ggrepel)

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

# Clean results
results <- results[, c("species", "bestModel", "droughtsIncluded", "changeSSP1-2.6", "changeSSP3-7.0", "changeSSP5-8.5", "changeSSP1-2.6NoDisp",
                       "changeSSP3-7.0NoDisp", "changeSSP5-8.5NoDisp")]

# Convert from wide to long format
df_long <- results %>%
  pivot_longer(cols = c("changeSSP1-2.6", "changeSSP3-7.0", "changeSSP5-8.5", "changeSSP1-2.6NoDisp",
               "changeSSP3-7.0NoDisp", "changeSSP5-8.5NoDisp"), names_to = "SSP_Disp", values_to = "RangeChange")

for (i in 1:nrow(df_long)) {
  if (df_long[i, "SSP_Disp"] %in% c("changeSSP1-2.6NoDisp", "changeSSP3-7.0NoDisp", "changeSSP5-8.5NoDisp")) df_long[i, "Dispersal"] <- "Dispersal: none"
  else df_long[i, "Dispersal"] <- "Dispersal: 100 km"
}

for (i in 1:nrow(df_long)) {
  if (df_long[i, "SSP_Disp"] %in% c("changeSSP1-2.6NoDisp", "changeSSP1-2.6")) df_long[i, "SSP"] <- "SSP1-2.6"
  if (df_long[i, "SSP_Disp"] %in% c("changeSSP3-7.0NoDisp", "changeSSP3-7.0")) df_long[i, "SSP"] <- "SSP3-7.0"
  if (df_long[i, "SSP_Disp"] %in% c("changeSSP5-8.5NoDisp", "changeSSP5-8.5")) df_long[i, "SSP"] <- "SSP5-8.5"
}

df_long <- df_long[, c("species", "bestModel", "droughtsIncluded", "RangeChange", "Dispersal", "SSP")]
df_long <- subset(df_long, !is.na(bestModel) & !is.na(RangeChange))

# Regress range changes
tmodel <- glmmTMB(RangeChange ~ bestModel + Dispersal + droughtsIncluded + SSP + (1 | species),
                  family = tweedie(link = "log"), data = df_long)
tmodelSumm <- summary(tmodel)
tmodelCI <- confint(tmodel)
exp(summary(tmodel)$coefficients$cond)
exp(tmodelCI)

# Create coefficient plot
tidymodel <- tidy(tmodel, effects = "fixed", conf.int = T) %>%
  mutate(estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
png(paste0(relPath, figPath, "Fig3.png"), units = "mm", width = 180, height = 180, res = 300)
dwplot(tidymodel, dot_args = list(color = "black"), whisker_args = list(color = "black")) + 
  theme_bw() + ylab("Independent variable") + xlab("Exponentiated coefficients") +
  scale_y_discrete(labels = c("SSP5-8.5", "SSP3-7.0", "With severe drought", "No dispersal", "GCM: UKESM1-0-LL", "GCM: MPI-ESM1-2-HR")) + theme(legend.position = "none") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + scale_x_continuous(breaks = seq(0, max(tidymodel$conf.high), by = 0.5),
                                                                                        limits = c(0, max(tidymodel$conf.high)))
dev.off()


# Boxplots of range size changes
df_long[df_long$droughtsIncluded, "Droughts"] <- "Yes"
df_long[!df_long$droughtsIncluded, "Droughts"] <- "No"
png(paste0(relPath, figPath, "Fig2.png"), units = "mm", width = 270, height = 180, res = 300)
ggplot(data = df_long, aes(y = RangeChange, x = Droughts, fill = Droughts)) +
  geom_violin(size = 1) + ylim(0, 5) + facet_grid(Dispersal ~ SSP) + geom_hline(yintercept = 1, linetype = "dashed", color = "black") + theme_bw() + 
  ylab("Predicted change in suitable climate area") + xlab("Severe drought variables in model") + theme(legend.position = "none")
dev.off()

# Summary statistics
tmp <- subset(df_long, Dispersal == "Dispersal: none" & SSP == "SSP1-2.6")
table(tmp$droughtsIncluded) / nrow(tmp)
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)
tmp <- subset(df_long, Dispersal == "Dispersal: 100 km" & SSP == "SSP1-2.6")
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)
tmp <- subset(df_long, Dispersal == "Dispersal: none" & SSP == "SSP3-7.0")
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)
tmp <- subset(df_long, Dispersal == "Dispersal: 100 km" & SSP == "SSP3-7.0")
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)
tmp <- subset(df_long, Dispersal == "Dispersal: none" & SSP == "SSP5-8.5")
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)
tmp <- subset(df_long, Dispersal == "Dispersal: 100 km" & SSP == "SSP5-8.5")
table(tmp$RangeChange < 1) / nrow(tmp)
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange < 0.20) / nrow(tmp)

# Extinction statistics
tmp <- subset(df_long, Dispersal == "Dispersal: 100 km" & SSP == "SSP1-2.6")
table(tmp$RangeChange < 0.75) / nrow(tmp)
table(tmp$RangeChange <= 0.5)
table(subset(tmp, droughtsIncluded)$RangeChange < 0.1)
table(subset(tmp, !droughtsIncluded)$RangeChange < 0.1)

# Export range change table for Supplementary Information
rangeChangesTable <- df_long[, c("species", "bestModel", "Droughts", "Dispersal", "SSP", "RangeChange")]
rangeChangesTable$species <- gsub(paste0(relPath, outputPath), "", df_long$species, fixed = T)
write.csv(rangeChangesTable, paste0(relPath, tablePath, "SupplementaryTable1.csv"), row.names = F)

# Summarize directional range changes by variable set, future GCM, and SSP
df_long$ChangeDirection <- ifelse(df_long$RangeChange >= 1, "+", "-")
table(subset(df_long, Droughts == "Yes")$ChangeDirection, subset(df_long, Droughts == "Yes")$SSP, subset(df_long, Droughts == "Yes")$Dispersal) / 550
table(subset(df_long, Droughts == "Yes")$ChangeDirection, subset(df_long, Droughts == "Yes")$SSP, subset(df_long, Droughts == "Yes")$Dispersal)
table(subset(df_long, Droughts == "No")$ChangeDirection, subset(df_long, Droughts == "No")$SSP, subset(df_long, Droughts == "No")$Dispersal) / 80
table(subset(df_long, Droughts == "No")$ChangeDirection, subset(df_long, Droughts == "No")$SSP, subset(df_long, Droughts == "No")$Dispersal)