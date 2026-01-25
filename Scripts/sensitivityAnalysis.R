library(terra)
library(data.table)
library(tidyr)
library(dplyr)
library(ENMeval)
library(glmmTMB)
library(ggplot2)

# Set map directory
mapDir <- "D:/Research/DroughtForecasts/OutputsSensitivityAnalysis/"

# List summary output files
range_dirs <- list.dirs(mapDir, recursive = TRUE, full.names = TRUE)
range_dirs <- range_dirs[basename(range_dirs) == "RangeSizes"]
files <- file.path(range_dirs, "rangeSizes.csv")
files <- files[file.exists(files)]

# Combine output files
all_data <- as.data.frame(rbindlist(lapply(files, function(f) { fread(f) }), fill = TRUE))

# Create parameter columns
all_data$basename.outNameBinary. <- sub("\\.tif$", "", all_data$basename.outNameBinary)
all_data <- separate(
  all_data,
  col = basename.outNameBinary.,
  into = c("VariableSet", "GCM", "Time", "BufferSize", "Clamp", "CropDist", "Threshold"),
  sep = "_",
  remove = F
)
all_data <- separate(
  all_data,
  col = sp,
  into = c("Species", "EnvOutliers1", "EnvOutliers2", "Thinning"),
  sep = "_",
  remove = F
)
all_data$EnvOutliers <- paste0(all_data$EnvOutliers1, "_", all_data$EnvOutliers2)

# Export model summary
write.csv(all_data, "D:/Research/DroughtForecasts/Manuscript/Submission2/Tables/SupplementaryTable4.csv",
          row.names = F)

# Get AUC and Boyce index summary statistics by buffer size
bufferCombo <- unique(all_data$sp)
bufferStats <- data.frame()
for (i in 1:length(bufferCombo)) {
  print(i)
  tmp <- subset(all_data, sp == bufferCombo[i])
  tmp <- subset(tmp, Time == "current" & Clamp == "clamp0" & CropDist == "crop0" & Threshold == "th-maxTSS")
  tmp <- tmp[, c("sp", "VariableSet", "GCM", "BufferSize", "AUC", "boyce")]
  wide_data <- tmp %>%
    pivot_wider(
      id_cols = c("sp", "VariableSet", "GCM"),
      names_from = BufferSize,
      values_from = c(AUC, boyce)
    )
  bufferStats <- rbind(bufferStats, wide_data)
}
bufferSummary <- bufferStats %>%
  summarize(
    mean_Boyce_100k = mean(boyce_buffer100000, na.rm = TRUE),
    mean_Boyce_500k = mean(boyce_buffer500000, na.rm = TRUE),
    median_Boyce_100k = median(boyce_buffer100000, na.rm = TRUE),
    median_Boyce_500k = median(boyce_buffer500000, na.rm = TRUE),
    sd_Boyce_100k = sd(boyce_buffer100000, na.rm = TRUE),
    sd_Boyce_500k = sd(boyce_buffer500000, na.rm = TRUE)
  )

# Regress Boyce index on buffer size
bufferData$boyceTrans <- ((bufferData$boyce + 1) / 2 * (nrow(bufferData) - 1) + 0.5) / nrow(bufferData) # Transform Boyce index to be strictly positive, followed by Smithson & Verkuilen (2006) transformation
bufferModelBoyce <- glmmTMB(boyceTrans ~ BufferSize + (1|Species) + (1|GCM) + (1|Thinning) + (1|EnvOutliers),
                 family = beta_family(), data = bufferData)
summary(bufferModelBoyce)
confint(bufferModelBoyce)

# Test sensitivity to outlier filtering
all_data$tmpNameOutlier <- paste(all_data$Species, all_data$Thinning, all_data$VariableSet, all_data$GCM, all_data$Clamp, sep = "_")
outlierCombo <- unique(all_data$tmpNameOutlier)
outlierCombo <- grepv("_ase_|_as_NA", outlierCombo)
nicheOverlaps <- data.frame(Comparison = character(),
                            D = numeric(),
                            I = numeric()) 
for (i in 1:length(outlierCombo)) {
  print(i)
  nicheOverlaps[i, "Comparison"] <- outlierCombo[i]
  tmp <- subset(all_data, tmpNameOutlier == outlierCombo[i])
  tmp <- subset(tmp, Time == "current" & CropDist == "crop500000" & BufferSize == "buffer500000" & Threshold == "th-maxTSS")
  if (nrow(tmp) != 0) {
    files <- paste0(mapDir, tmp$sp, "/SuitabilityMaps/Current/", tmp$VariableSet, "_", tmp$GCM, "_", tmp$Time, "_", tmp$BufferSize, "_", tmp$Clamp, ".tif")
    r <- tryCatch({
      rast(files) / 10000
    }, error = function(e) {
      if (grepl("extents do not match", e$message)) {
        tmpRast1 <- rast(files[1]) / 10000
        tmpRast2 <- rast(files[2]) / 10000
        tmpRast2 <- resample(tmpRast2, tmpRast1, method = "bilinear")
        c(tmpRast1, tmpRast2)
      } else {
        stop(e)
      }
    })
    overlapD <- calc.niche.overlap(r, overlapStat = "D", quiet = TRUE)[2, 1]
    overlapI <- calc.niche.overlap(r, overlapStat = "I", quiet = TRUE)[2, 1]
    nicheOverlaps[i, "D"] <- overlapD
    nicheOverlaps[i, "I"] <- overlapI
  }
}

# Get summary statistics
summary(nicheOverlaps$D)
sd(nicheOverlaps$D, na.rm = T)
summary(nicheOverlaps$I)
sd(nicheOverlaps$I, na.rm = T)

# Plot niche overlaps
df_long <- nicheOverlaps %>%
  pivot_longer(cols = c(D, I), names_to = "Metric", values_to = "Value")
df_long$Metric <- as.factor(df_long$Metric)
levels(df_long$Metric) <- c("Schoener's *D*", "*I* similarity statistic")
png("D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/ExtendedDataFig8.png", units = "mm", width = 270, height = 180, res = 300)
ggplot(df_long, aes(x = Value, y = after_stat(count))) +
  geom_histogram(bins = 24, fill = "grey70", color = "black") +
  facet_wrap(~ Metric, nrow = 1) +
  theme_bw() +
  theme(strip.text = element_text(size = 11), legend.title = element_blank(), axis.title = element_text(size = 12)) +
  xlab("Niche overlap") +
  ylab("Count") +
  xlim(0.4, NA) +
  theme(strip.text = ggtext::element_markdown())
dev.off()