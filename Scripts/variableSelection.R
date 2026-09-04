library(ggplot2)
library(tidyverse)
library(ggtext)

# Set base directory
relDir <- "D:/Research/DroughtForecasts/Outputs/"

# Set thresholds
aucThreshold <- 0.7 # Threshold for discarding models based on AUC
boyceThreshold <- 0.5 # Threshold for discarding models based on Boyce index

# Find all files
dirs <- list.dirs(relDir, full.names = TRUE, recursive = FALSE)
files <- file.path(dirs, "ModelSetComparison", "modelComparison.csv")
files <- files[file.exists(files)]

# Combine all files
df_list <- lapply(files, read.csv)
result <- do.call(rbind, df_list)

# Isolate thinning distance
result$thin <- as.integer(sub(".*_thin([0-9]+)$", "\\1", result$sp))

# Get best model statistics
resultSubset <- subset(result, AUC > aucThreshold & boyce > boyceThreshold)
best_models <- resultSubset %>%
  group_by(sp) %>%
  slice_min(BIC, n = 1, with_ties = FALSE) %>%
  ungroup()
model_fractions <- best_models %>%
  count(model) %>%
  mutate(
    fraction = n / sum(n)
  )

best_models_thin <- resultSubset %>%
  group_by(sp, thin) %>%
  slice_min(BIC, n = 1, with_ties = FALSE) %>%
  ungroup()
model_fractions_thin <- best_models_thin %>%
  count(thin, model) %>%
  group_by(thin) %>%
  mutate(
    fraction = n / sum(n)
  ) %>%
  ungroup()

# Count frequencies by model type
varNames <- c(paste0("BIO", 1:19),
              "HAS_mean", "HAD_mean", "HMS_mean", "HMD_mean",
              "HAS_median", "HAD_median", "HMS_median", "HMD_median",
              "HAS_q0.75", "HAD_q0.75", "HMS_q0.75", "HMD_q0.75",
              "SOC", "pH", "silt", "sand", "clay", "CEC", "CFVO")
freqs <- data.frame(model = unique(result$model))
freqs[varNames] <- 0

for (i in 1:nrow(freqs)) {
  tmpModel <- freqs[i, "model"]
  modelSubset <- subset(result, model == tmpModel)
  modelSubset <- subset(modelSubset, AUC > aucThreshold)
  modelSubset <- subset(modelSubset, boyce > boyceThreshold)
  freqs[i, "noModels"] <- nrow(modelSubset)
  for (j in 1:nrow(modelSubset)) {
    tmpVars <- modelSubset[j, "varnamesfit"]
    tmpVars <- unlist(strsplit(tmpVars, " - "))
    for (k in 1:length(tmpVars)) {
      freqs[i, tmpVars[k]] <- freqs[i, tmpVars[k]] + 1
    }
  }
}

# Normalize by number of models
freqs[, 2:39] <- freqs[, 2:39] / freqs[, 40]

# Plot frequencies

# Convert to long format
freqs_long <- freqs %>%
  pivot_longer(
    cols = -c(model, noModels),
    names_to = "variable",
    values_to = "frequency"
  )

# Order variables by mean frequency
freqs_long <- freqs_long %>%
  mutate(variable = fct_reorder(variable, frequency))

# Rename model types
freqs_long$model <- as.factor(freqs_long$model)
levels(freqs_long$model) <- c("Without severe drought",
                              "With severe drought (GFDL-ESM4)",
                              "With severe drought (MPI-ESM1-2-HR)",
                              "With severe drought (UKESM1-0-LL)")

# Add variable types
freqs_long <- freqs_long %>%
  mutate(
    varType = case_when(
      str_detect(variable, "BIO") ~ "Bioclimatic",
      str_detect(variable, "HMS|HMD|HAS|HAD") ~ "Severe drought",
      TRUE ~ "Soil"
    )
  )

# Rename variables
lvls <- levels(freqs_long$variable)
lvls <- gsub("_mean", " (mean)", lvls)
lvls <- gsub("_median", " (median)", lvls)
lvls <- gsub("_q0.75", " (75th percentile)", lvls)
lvls <- gsub("sand", "Sand", lvls)
lvls <- gsub("silt", "Silt", lvls)
lvls <- gsub("clay", "Clay", lvls)
levels(freqs_long$variable) <- lvls

# Color drought variables
label_colors <- freqs_long %>%
  group_by(variable) %>%
  summarize(is_severe = any(varType == "Severe drought")) %>%
  mutate(color = ifelse(is_severe, "red", "black")) %>%
  pull(color)

# Exclude drought variable rows for models without drought variables
freqs_long <- subset(freqs_long, !(model == "Without severe drought" & varType == "Severe drought"))

# Single-panel Cleveland dot plot, colored by model
cairo_pdf("D:/Research/DroughtForecasts/Manuscript/Revision2/Figures/Fig1.pdf",
          width = 180/25.4, height = 170/25.4)
ggplot(freqs_long, aes(x = frequency, y = variable, color = model, shape = varType)) +
  geom_point(size = 2) +
  labs(
    x = "Fraction of models",
    y = "Variable",
    color = "Model type",
    shape = "Variable type"
  ) +
  geom_richtext(inherit.aes = F,
    x = 0.425, y = levels(freqs_long$variable)[33],
    label = 
      "<b>BIO3</b>: isothermality<br>
      <b>BIO2</b>: mean diurnal temperature range<br>
      <b>HMS</b>: mean drought severity<br>
    <b>HAS</b>: accumulated drought severity<br>
    <b>BIO8</b>: mean temperature wettest quarter<br>
    <b>BIO9</b>: mean temperature driest quarter<br>
    <b>HMD</b>: mean drought duration<br>
    <b>HAD</b>: accumulated drought duration<br>
    <b>pH</b>: soil pH<br>
    <b>BIO7</b>: annual temperature range<br>
    <b>BIO14</b>: precipitation driest month<br>
    <b>BIO15</b>: precipitation seasonality<br>
    <b>BIO5</b>: maximum temperature warmest month<br>
    <b>CEC</b>: cation exchange capacity<br>
    <b>BIO19</b>: precipitation coldest quarter<br>
    <b>CFVO</b>: coarse fragments fraction<br>
    <b>BIO1</b>: mean annual temperature<br>
    <b>BIO4</b>: temperature seasonality<br>
    <b>BIO6</b>: minimum temperature coldest month<br>
    <b>Silt</b>: soil silt fraction<br>
    <b>Clay</b>: soil clay fraction<br>
    <b>BIO11</b>: mean temperature coldest quarter<br>
    <b>SOC</b>: soil organic carbon<br>
    <b>BIO10</b>: mean temperature warmest quarter<br>
    <b>BIO18</b>: precipitation warmest quarter<br>
    <b>BIO13</b>: precipitation wettest month<br>
    <b>Sand</b>: soil sand fraction<br>
    <b>BIO17</b>: precipitation driest quarter<br>
    <b>BIO16</b>: precipitation wettest quarter<br>
    <b>BIO12</b>: annual precipitation",
    hjust = 0, vjust = 1, size = 2.5,
    fill = "white", label.color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.5, 0.120),
    legend.box = "horizontal",
    axis.text.y = element_text(color = label_colors)
  )
dev.off()

# Export data on included variables
for (i in 1:nrow(resultSubset)) {
  print(i)
  tmp <- strsplit(resultSubset[i, "model"], "_", fixed = T)
  if (tmp[[1]][1] == "as") resultSubset[i, "Variable set"] <- "Without drought"
  if (tmp[[1]][1] == "ase") resultSubset[i, "Variable set"] <- "With drought"
  resultSubset[i, "GCM"] <- tmp[[1]][2]
  resultSubset[i, "Taxon"] <- strsplit(resultSubset[i, "sp"], "_", fixed = T)[[1]][1]
}
write.csv(resultSubset, "D:/Research/DroughtForecasts/Manuscript/Submission2/Tables/SupplementaryTable3.csv", row.names = F)