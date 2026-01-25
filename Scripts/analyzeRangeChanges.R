library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(broom.mixed)
library(dotwhisker)
library(ggrepel)

# Set base directory
relDir <- "D:/Research/DroughtForecasts/Outputs/"

# Find all files
dirs <- list.dirs(relDir, full.names = TRUE, recursive = FALSE)
files <- file.path(dirs, "RangeSizes", "rangeChanges.csv")
files <- files[file.exists(files)]

# Combine all files
df_list <- lapply(files, read.csv)
result <- do.call(rbind, df_list)

# Split parameters
result <- result %>%
  separate(
    col = basename.outNameBinary.,
    into = c("variableSet", "GCM", "time", "ignore1", "clamp", "cropDistance", "thresholdType"),
    sep = "_",
    remove = FALSE
  ) %>%
  select(-ignore1)
result <- result %>%
  separate(
    col = sp,
    into = c("taxon", "envFilter", "extFilter", "thinDistance"),
    sep = "_",
    remove = F
  ) %>%
  select(!c(envFilter, extFilter))

# Plot range changes (Clamp: yes; thinning: 0 km; threshold: maxTSS; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin0" & thresholdType == "th-maxTSS.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/Fig2.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: no; thinning: 0 km; threshold: maxTSS; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp0" & thinDistance == "thin0" & thresholdType == "th-maxTSS.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig1.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: yes; thinning: 5 km; threshold: maxTSS; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin5" & thresholdType == "th-maxTSS.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig2.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: yes; thinning: 10 km; threshold: maxTSS; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin10" & thresholdType == "th-maxTSS.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig3.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: yes; thinning: 0 km; threshold: maxTSS; EPV > 10)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin0" & thresholdType == "th-maxTSS.tif" & EPV > 10)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig4.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: yes; thinning: 0 km; threshold: eq. sens. spec.; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin0" & thresholdType == "th-eqsensspec.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig5.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Plot range changes (Clamp: yes; thinning: 0 km; threshold: OR10; EPV > 0)
resultsPlot <- subset(result, time != "current" & clamp == "clamp1" & thinDistance == "thin0" & thresholdType == "th-OR10.tif" & EPV > 0)
resultsPlot$variableSet <- as.factor(resultsPlot$variableSet)
levels(resultsPlot$variableSet) <- c("No severe droughts", "With severe droughts")
resultsPlot$cropDistance <- as.factor(resultsPlot$cropDistance)
levels(resultsPlot$cropDistance) <- c("Dispersal: none", "Dispersal: 100 km", "Dispersal: 500 km")
resultsPlot$time <- as.factor(resultsPlot$time)
levels(resultsPlot$time) <- c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig6.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(resultsPlot, aes(x = rangeChange_sameCrop, linetype = variableSet, col = GCM)) +
  geom_density(size = 0.75, adjust = 0.4) + facet_grid(cropDistance ~ time) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + theme_bw() +
  xlab("Predicted fraction of suitable climatic area remaining in the future") + ylab("Density") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 12, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 12, alpha = 0.05, fill = "green") +
  theme(strip.text = element_text(size = 11), axis.title = element_text(size = 12)) + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(linetype = "Model type") +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      "No severe droughts"   = "dashed",
      "With severe droughts" = "solid"
    )
  ) + coord_cartesian(xlim = c(0, 2), ylim = c(0, 12)) + scale_x_continuous(trans = "sqrt", breaks = c(0, 0.1, 0.25, 0.5, 1, 1.5, 2))
dev.off()

# Calculate summary statistics
summary_df_thin0_clamp1 <- result %>%
  filter(
    time != "current",
    thinDistance == "thin0",
    clamp == "clamp1"
  ) %>%
  group_by(
    taxon,
    time,
    cropDistance, variableSet
  ) %>%
  summarize(
    rangeChange = mean(rangeChange_sameCrop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    time,
    cropDistance, variableSet
  ) %>%
  summarize(
    frac_lt_1    = mean(rangeChange < 1,    na.rm = TRUE),
    frac_le_0.2  = mean(rangeChange <= 0.2, na.rm = TRUE),
    frac_le_0.75 = mean(rangeChange <= 0.75, na.rm = TRUE),
    median_range = median(rangeChange, na.rm = TRUE),
    n_taxa       = n(),
    .groups = "drop"
  )
write.csv(summary_df_thin0_clamp1, "D:/Research/DroughtForecasts/Manuscript/Submission2/Tables/Table1.csv", row.names = F)

summary_df_all <- result %>%
  filter(time != "current") %>%
  group_by(
    taxon,
    time,
    cropDistance,
    thinDistance,
    clamp,
    variableSet
  ) %>%
  summarize(
    rangeChange = mean(rangeChange_sameCrop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    time,
    cropDistance,
    thinDistance,
    clamp,
    variableSet
  ) %>%
  summarize(
    frac_lt_1    = mean(rangeChange < 1,     na.rm = TRUE),
    frac_le_0.2  = mean(rangeChange <= 0.2,  na.rm = TRUE),
    frac_le_0.75 = mean(rangeChange <= 0.75, na.rm = TRUE),
    median_range = median(rangeChange, na.rm = TRUE),
    n_taxa       = n(),
    .groups = "drop"
  ) %>%
  group_by(variableSet) %>%
  summarize(
    frac_lt_1    = mean(frac_lt_1,    na.rm = TRUE),
    frac_le_0.2  = mean(frac_le_0.2,  na.rm = TRUE),
    frac_le_0.75 = mean(frac_le_0.75, na.rm = TRUE),
    median_range = mean(median_range, na.rm = TRUE),
    n_taxa       = mean(n_taxa),
    .groups = "drop"
  )

# Export range changes
write.csv(subset(result, time != "current"), "D:/Research/DroughtForecasts/Manuscript/Submission2/Tables/SupplementaryTable3.csv", row.names = F)

# Fit regression
resultsRegression <- subset(result, time != "current")
tmodel <- glmmTMB(rangeChange_sameCrop ~ variableSet + cropDistance + clamp + time + thresholdType + thinDistance + GCM
                  + (1 | taxon),
                  family = tweedie(link = "log"), data = resultsRegression)
tmodelSumm <- summary(tmodel)
tmodelCI <- confint(tmodel)
exp(summary(tmodel)$coefficients$cond)
exp(tmodelCI)

# Create coefficient plot
tidymodel <- tidy(tmodel, effects = "fixed", conf.int = T) %>%
  mutate(estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/Fig3.png", units = "mm", width = 180, height = 120, res = 300)
dwplot(tidymodel, dot_args = list(color = "black"), whisker_args = list(color = "black")) + 
  theme_bw() + ylab("Independent variable") + xlab("Exponentiated coefficients") + xlim(0.75, 2) +
  scale_y_discrete(labels = c("GCM: UKESM1-0-LL", "GCM: MPI-ESM1-2-HR", "Thinning: 5 km", "Thinning: 10 km", "Threshold: 10th percentile omission rate", "Threshold: maximum true skill statistic", "SSP5-8.5", "SSP3-7.0", "With clamping", "Dispersal: 500 km", "Dispersal: 100 km", "With severe drought")) + 
  theme(legend.position = "none", axis.title = element_text(size = 12)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + scale_x_continuous(breaks = seq(0.75, 2, by = 0.25),
                                                                                        limits = c(0.75, 2)) +
  annotate("rect", xmin = 0.75, xmax = 1, ymin = 0.5, ymax = 12.5, alpha = 0.05, fill = "red") + 
  annotate("rect", xmin = 1, xmax = 2, ymin = 0.5, ymax = 12.5, alpha = 0.05, fill = "green")
dev.off()