library(ggplot2)
library(tidyverse)

# Set base directory
relDir <- "D:/Research/DroughtForecasts/Outputs/"

# Find all files
dirs <- list.dirs(relDir, full.names = TRUE, recursive = FALSE)
files <- file.path(dirs, "RangeSizes", "rangeSizes.csv")
files <- files[file.exists(files)]

# Calculate range size changes
for (i in 1:length(files)) {
  print(i)
  tmpFile <- read.csv(files[i])
  for (j in 1:nrow(tmpFile)) {
    tmpRow <- tmpFile[j, ]
    tmpTime <- strsplit(tmpRow$basename.outNameBinary., "_", fixed = T)[[1]][3]
    if (tmpTime == "current") next # Skip if not a future projection
    # Find range size for present
    tmpClamp <- strsplit(tmpRow$basename.outNameBinary., "_", fixed = T)[[1]][5]
    tmpCrop <- strsplit(tmpRow$basename.outNameBinary., "_", fixed = T)[[1]][6]
    tmpThreshold <- strsplit(tmpRow$basename.outNameBinary., "_", fixed = T)[[1]][7]
    tmpSubset <- gsub(".tif", "", tmpRow$model, fixed = T)
    tmpSubset_crop500 <- paste0(tmpSubset, "_", tmpClamp, "_crop500000_", tmpThreshold)
    tmpSubset_sameCrop <- paste0(tmpSubset, "_", tmpClamp, "_", tmpCrop, "_", tmpThreshold)
    tmpBaseline_crop500 <- subset(tmpFile, basename.outNameBinary. == tmpSubset_crop500)$rangeSize
    tmpBaseline_sameCrop <- subset(tmpFile, basename.outNameBinary. == tmpSubset_sameCrop)$rangeSize
    # Calculate change
    tmpFile[j, "currentRangeSize_crop500"] <- tmpBaseline_crop500
    tmpFile[j, "currentRangeSize_sameCrop"] <- tmpBaseline_sameCrop
    tmpFile[j, "rangeChange_crop500"] <- tmpRow$rangeSize / tmpBaseline_crop500
    tmpFile[j, "rangeChange_sameCrop"] <- tmpRow$rangeSize / tmpBaseline_sameCrop
  }
  outName <- gsub("rangeSizes.csv", "rangeChanges.csv", files[i], fixed = T)
  write.csv(tmpFile, outName, row.names = F)
}