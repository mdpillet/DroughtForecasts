library(dplyr)
library(tidyr)

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

# Decompose variance
anov <- summary(aov(rangeChange_sameCrop ~ 
                           variableSet + 
                           GCM +
                           time +
                           clamp + 
                           cropDistance + 
                           thresholdType +
                           thinDistance + variableSet:cropDistance, 
                         data = subset(result, time != "current")))
(anov[[1]]$`Sum Sq` / sum(anov[[1]]$`Sum Sq`[1:8]) * 100)[1:8]