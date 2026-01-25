library(data.table)
library(glmmTMB)
library(ggplot2)

# Set model directory
modelDir <- "D:/Research/DroughtForecasts/Outputs/"

# List summary output files
range_dirs <- list.dirs(modelDir, recursive = TRUE, full.names = TRUE)
range_dirs <- range_dirs[basename(range_dirs) == "ModelSetComparison"]
files <- file.path(range_dirs, "modelComparison.csv")
files <- files[file.exists(files)]

# Combine output files
all_data <- as.data.frame(rbindlist(lapply(files, function(f) { fread(f) }), fill = TRUE))
all_data <- subset(all_data, ncoeffit > 0)

# Transform Boyce index to be strictly positive, followed by Smithson & Verkuilen (2006) transformation
all_data$boyceTrans <- ((all_data$boyce + 1) / 2 * (nrow(all_data) - 1) + 0.5) / nrow(all_data)

# Regress AUC on variable set and thinning radius
all_data$thinningRadius <- unlist(lapply(strsplit(all_data$sp, "_", fixed = T), "[", 4))
all_data$taxon <- unlist(lapply(strsplit(all_data$sp, "_", fixed = T), "[", 1))
AUCmodel <- glmmTMB(AUC ~ model + thinningRadius + (1 | taxon), family = beta_family(), data = all_data)
summary(AUCmodel)
round(confint(AUCmodel), 3)

# Regress Boyce index on variable set and thinning radius
boyceModel <- glmmTMB(boyceTrans ~ model + thinningRadius + (1 | taxon), family = beta_family(), data = all_data)
summary(boyceModel)
round(confint(boyceModel), 3)

# Subset models
ntaxaPre <- length(unique(all_data$taxon))
# Export summary table
write.csv(all_data[, c("taxon", "model", "thinningRadius", "noPres", "lambdaRule", "nvarraw", "ncoefraw", "nvarfit", "ncoeffit", "EPV", "AUC", "boyce", "BIC")],
          "D:/Research/DroughtForecasts/Manuscript/Submission2/Tables/SupplementaryTable1.csv", row.names = F)
all_data <- subset(all_data, AUC > 0.7 & boyce > 0.5)
ntaxaPost <- length(unique(all_data$taxon))
length(unique(subset(all_data, thinningRadius == "thin0")$taxon))
length(unique(subset(all_data, thinningRadius == "thin5")$taxon))
length(unique(subset(all_data, thinningRadius == "thin10")$taxon))

# Get summary statistics for AUC and Boyce index
summary(all_data$AUC)
sd(all_data$AUC, na.rm = T)
tapply(all_data$AUC, all_data$model, summary)
tapply(all_data$AUC, all_data$model, sd, na.rm = T)
summary(all_data$boyce)
sd(all_data$boyce, na.rm = T)
tapply(all_data$boyce, all_data$model, summary)
tapply(all_data$boyce, all_data$model, sd, na.rm = T)

# Get summary statistics for number of coefficients
tapply(all_data$ncoeffit, all_data$model, mean, na.rm = T)
tapply(all_data$ncoeffit, all_data$model, sd, na.rm = T)
tapply(all_data$ncoeffit / all_data$ncoefraw, all_data$model, mean, na.rm = T)
tapply(all_data$ncoeffit / all_data$ncoefraw, all_data$model, sd, na.rm = T)