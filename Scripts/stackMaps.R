library(terra)
library(dplyr)
library(readr)
library(purrr)
library(scico)
library(fields)

# Set directory structure and load workflow functions
relDir <- "D:/Research/DroughtForecasts/Outputs/"
outPath <- "D:/Research/DroughtForecasts/Data/DiversityMaps/"
countryPath <- "D:/Research/DroughtForecasts/Data/CountriesOutline/world-administrative-boundaries.shp"
figPath <- "D:/Research/DroughtForecasts/Manuscript/Submission2/Figures/"

# Set thresholds
aucThreshold <- 0.7 # Threshold for discarding models based on AUC
boyceThreshold <- 0.5 # Threshold for discarding models based on Boyce index

# Find all model comparison files
dirs <- list.dirs(relDir, full.names = TRUE, recursive = FALSE)
files <- file.path(dirs, "ModelSetComparison", "modelComparison.csv")
files <- files[file.exists(files)]
files <- files[grepl("_thin0", files, fixed = T)]

# Find best models for each species exceeding AUC and Boyce index thresholds
best_models <- map_dfr(files, function(f) {
  df <- read_csv(f, show_col_types = FALSE) %>%
    filter(
      AUC   > aucThreshold,
      boyce > boyceThreshold
    )
  if (nrow(df) == 0) return(NULL)
  df %>%
    arrange(BIC, desc(AUC)) %>%
    slice(1) %>%
    select(sp, model)
})

# Create current diversity maps
paramCurrent <- c(
  "_clamp0_crop0_th-maxTSS",
  "_clamp0_crop0_th-eqsensspec",
  "_clamp0_crop0_th-OR10",
  "_clamp0_crop100000_th-maxTSS",
  "_clamp0_crop100000_th-eqsensspec",
  "_clamp0_crop100000_th-OR10",
  "_clamp0_crop500000_th-maxTSS",
  "_clamp0_crop500000_th-eqsensspec",
  "_clamp0_crop500000_th-OR10",
  "_clamp1_crop0_th-maxTSS",
  "_clamp1_crop0_th-eqsensspec",
  "_clamp1_crop0_th-OR10",
  "_clamp1_crop100000_th-maxTSS",
  "_clamp1_crop100000_th-eqsensspec",
  "_clamp1_crop100000_th-OR10",
  "_clamp1_crop500000_th-maxTSS",
  "_clamp1_crop500000_th-eqsensspec",
  "_clamp1_crop500000_th-OR10"
)
stackFiles <- NULL
for (i in 1:nrow(best_models)) {
  tmpDir <- paste0(relDir, best_models[i, "sp"], "/BinaryMaps/Current/")
  tmpFiles <- list.files(tmpDir, "tif$", full.names = T)
  tmpModel <- strsplit(as.character(best_models[i, "model"]), ".tif", fixed = T)[[1]][1]
  tmpFiles <- tmpFiles[grepl(tmpModel, tmpFiles, fixed = T)]
  stackFiles <- c(stackFiles, tmpFiles)
}
for (j in 1:length(paramCurrent)) {
  tmpStack <- stackFiles[grepl(paramCurrent[j], stackFiles, fixed = T)]
  tmp <- lapply(tmpStack, rast)
  tmp <- sprc(tmp)
  stacked <- mosaic(tmp, fun = "sum")
  writeRaster(stacked, paste0(outPath, "Current/", paramCurrent[j], ".tif"), overwrite = T)
  png(filename = paste0(outPath, "Current/", paramCurrent[j], ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(stacked, range = c(0, 150))
  dev.off()
}

# Calculate average current diversity
tmpCurrent <- list.files(paste0(outPath, "Current/"), "tif$", full.names = T)
tmpCurrent <- lapply(tmpCurrent, rast)
tmpCurrent <- sprc(tmpCurrent)
stackedCurrent <- mosaic(tmpCurrent, fun = "mean")
writeRaster(stackedCurrent, paste0(outPath, "Current/avg.tif"), overwrite = T)
png(filename = paste0(outPath, "Current/avg.png"), width = 10, height = 10, units = "in", res = 320)
plot(stackedCurrent, range = c(0, 100))
dev.off()

# Create future diversity maps
dir.create(paste0(outPath, "Temp/"))
paramFuture <- c(
  "_SSP126_buffer500000_clamp0_crop0_th-maxTSS", "_SSP126_buffer500000_clamp0_crop0_th-eqsensspec", "_SSP126_buffer500000_clamp0_crop0_th-OR10",
  "_SSP126_buffer500000_clamp0_crop100000_th-maxTSS", "_SSP126_buffer500000_clamp0_crop100000_th-eqsensspec", "_SSP126_buffer500000_clamp0_crop100000_th-OR10",
  "_SSP126_buffer500000_clamp0_crop500000_th-maxTSS", "_SSP126_buffer500000_clamp0_crop500000_th-eqsensspec", "_SSP126_buffer500000_clamp0_crop500000_th-OR10",
  "_SSP126_buffer500000_clamp1_crop0_th-maxTSS", "_SSP126_buffer500000_clamp1_crop0_th-eqsensspec", "_SSP126_buffer500000_clamp1_crop0_th-OR10",
  "_SSP126_buffer500000_clamp1_crop100000_th-maxTSS", "_SSP126_buffer500000_clamp1_crop100000_th-eqsensspec", "_SSP126_buffer500000_clamp1_crop100000_th-OR10",
  "_SSP126_buffer500000_clamp1_crop500000_th-maxTSS", "_SSP126_buffer500000_clamp1_crop500000_th-eqsensspec", "_SSP126_buffer500000_clamp1_crop500000_th-OR10",
  "_SSP370_buffer500000_clamp0_crop0_th-maxTSS", "_SSP370_buffer500000_clamp0_crop0_th-eqsensspec", "_SSP370_buffer500000_clamp0_crop0_th-OR10",
  "_SSP370_buffer500000_clamp0_crop100000_th-maxTSS", "_SSP370_buffer500000_clamp0_crop100000_th-eqsensspec", "_SSP370_buffer500000_clamp0_crop100000_th-OR10",
  "_SSP370_buffer500000_clamp0_crop500000_th-maxTSS", "_SSP370_buffer500000_clamp0_crop500000_th-eqsensspec", "_SSP370_buffer500000_clamp0_crop500000_th-OR10",
  "_SSP370_buffer500000_clamp1_crop0_th-maxTSS", "_SSP370_buffer500000_clamp1_crop0_th-eqsensspec", "_SSP370_buffer500000_clamp1_crop0_th-OR10",
  "_SSP370_buffer500000_clamp1_crop100000_th-maxTSS", "_SSP370_buffer500000_clamp1_crop100000_th-eqsensspec", "_SSP370_buffer500000_clamp1_crop100000_th-OR10",
  "_SSP370_buffer500000_clamp1_crop500000_th-maxTSS", "_SSP370_buffer500000_clamp1_crop500000_th-eqsensspec", "_SSP370_buffer500000_clamp1_crop500000_th-OR10",
  "_SSP585_buffer500000_clamp0_crop0_th-maxTSS", "_SSP585_buffer500000_clamp0_crop0_th-eqsensspec", "_SSP585_buffer500000_clamp0_crop0_th-OR10",
  "_SSP585_buffer500000_clamp0_crop100000_th-maxTSS", "_SSP585_buffer500000_clamp0_crop100000_th-eqsensspec", "_SSP585_buffer500000_clamp0_crop100000_th-OR10",
  "_SSP585_buffer500000_clamp0_crop500000_th-maxTSS", "_SSP585_buffer500000_clamp0_crop500000_th-eqsensspec", "_SSP585_buffer500000_clamp0_crop500000_th-OR10",
  "_SSP585_buffer500000_clamp1_crop0_th-maxTSS", "_SSP585_buffer500000_clamp1_crop0_th-eqsensspec", "_SSP585_buffer500000_clamp1_crop0_th-OR10",
  "_SSP585_buffer500000_clamp1_crop100000_th-maxTSS", "_SSP585_buffer500000_clamp1_crop100000_th-eqsensspec", "_SSP585_buffer500000_clamp1_crop100000_th-OR10",
  "_SSP585_buffer500000_clamp1_crop500000_th-maxTSS", "_SSP585_buffer500000_clamp1_crop500000_th-eqsensspec", "_SSP585_buffer500000_clamp1_crop500000_th-OR10"
)

stackFiles <- NULL
for (i in 1:nrow(best_models)) {
  tmpDir <- paste0(relDir, best_models[i, "sp"], "/BinaryMaps/Future/")
  tmpFiles <- list.files(tmpDir, "tif$", full.names = T)
  tmpModel <- strsplit(as.character(best_models[i, "model"]), "_current_", fixed = T)[[1]][1]
  if (tmpModel == "as_NA") tmpModel <- "/as_"
  tmpFiles <- tmpFiles[grepl(tmpModel, tmpFiles, fixed = T)]
  if (tmpModel == "/as_") {
    stackFilesAS <- NULL
    for (j in 1:length(paramFuture)) {
      tmpGCM <- tmpFiles[grepl(paramFuture[j], tmpFiles, fixed = T)]
      tmpGCM <- lapply(tmpGCM, rast)
      tmpGCM <- sprc(tmpGCM)
      stackedGCM <- mosaic(tmpGCM, fun = "mean")
      writeRaster(stackedGCM, paste0(outPath, "Temp/", i, paramFuture[j], ".tif"), overwrite = T)
      stackFilesAS <- c(stackFilesAS, paste0(outPath, "Temp/", i, paramFuture[j], ".tif"))
    }
    stackFiles <- c(stackFiles, stackFilesAS)
    
  }
  else {
    stackFiles <- c(stackFiles, tmpFiles)
  }
}
for (j in 1:length(paramFuture)) {
  print(j)
  tmpStack <- stackFiles[grepl(paramFuture[j], stackFiles, fixed = T)]
  tmp <- lapply(tmpStack, rast)
  tmp <- sprc(tmp)
  stacked <- mosaic(tmp, fun = "sum")
  writeRaster(stacked, paste0(outPath, "Future/", paramFuture[j], ".tif"), overwrite = T)
  png(filename = paste0(outPath, "Future/", paramFuture[j], ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(stacked, range = c(0, 150))
  dev.off()
}
unlink(paste0(outPath, "Temp/"), recursive = TRUE)

# Calculate average future diversity
tmpFuture <- list.files(paste0(outPath, "Future/"), "tif$", full.names = T)
tmpFuture <- lapply(tmpFuture, rast)
tmpFuture <- sprc(tmpFuture)
stackedFuture <- mosaic(tmpFuture, fun = "mean")
writeRaster(stackedFuture, paste0(outPath, "Future/avg.tif"), overwrite = T)
png(filename = paste0(outPath, "Future/avg.png"), width = 10, height = 10, units = "in", res = 320)
plot(stackedFuture, range = c(0, 100))
dev.off()

# Create overview figure
richnessCurrent <- rast(paste0(outPath, "Current/avg.tif"))
richnessFuture <- rast(paste0(outPath, "Future/avg.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "Fig4.png"), units = "mm", width = 180, height = 180, res = 300)
# --- COLOR PREP FOR A & B---
my_cols <- c(
  "white",   "#001959", "#031D5A", "#06215B", "#07245A", "#0A285B", "#0A2C5C", 
  "#0C2F5C", "#0D335D", "#0D365D", "#0E395E", "#0F3C5F", "#0F3F60", "#11415F", 
  "#114460", "#114761", "#134961", "#144C62", "#144E62", "#165061", "#165261", 
  "#185561", "#1A5761", "#1B5961", "#1E5C62", "#1F5D60", "#226060", "#246160", 
  "#27635F", "#2B655D", "#2E665C", "#31695A", "#346959", "#386B57", "#3B6C55", 
  "#406E53", "#447051", "#47704F", "#4C724C", "#4F734B", "#537548", "#587645", 
  "#5C7744", "#607941", "#64793F", "#687A3D", "#6C7B3B", "#727D38", "#767F36", 
  "#7A7F34", "#7F8132", "#838231", "#88842F", "#8D842E", "#92862C", "#97882C", 
  "#9C892A", "#A18A2B", "#A68A2C", "#AB8B2D", "#B08C2E", "#B68D30", "#BB8E33", 
  "#BF9035", "#C49138", "#C9913B", "#CD923F", "#D29243", "#D79347", "#DB944B", 
  "#DF9550", "#E39654", "#E69759", "#EA995E", "#EE9B64", "#F09C6A", "#F39E70", 
  "#F59F76", "#F8A17B", "#F8A281", "#FAA587", "#FCA78E", "#FBA894", "#FCAB9A", 
  "#FCAC9F", "#FCAFA5", "#FDB0AA", "#FCB2AF", "#FDB4B6", "#FDB6BB", "#FDB8C0", 
  "#FCB9C6", "#FDBCCB", "#FBBDD0", "#FCBFD7", "#FBC2DC", "#FBC3E2", "#FAC6E8", 
  "#FAC7EE", "#FBCAF4", "#F9CCF9"
)
my_breaks <- c(-0.00001, 0.00001, seq(1, 100, by = 1))
# --- GLOBAL SETTINGS ---
par(mfrow = c(2, 2), mar = c(1, 0, 1, 0), oma = c(1, 2, 1, 1))
# --- PLOT A: CURRENT RICHNESS ---
plot(richnessCurrent, 
     box = FALSE, axes = FALSE, col = my_cols, breaks = my_breaks,
     legend = FALSE, mar = c(1, 1, 1, 4), reset = TRUE) 
plot(countries, add = TRUE)
mtext(expression(bold("a")), side = 3, adj = 0, line = 0.5)
fields::image.plot(
  zlim = c(0, 100), legend.only = TRUE, col = my_cols,
  smallplot = c(0.80, 0.82, 0.2, 0.8), 
  axis.args = list(at = seq(0, 100, 20), cex.axis = 0.8),
  legend.args = list(text = "Present richness\n(species)", side = 3, font = 1, line = 0.5, cex = 0.7)
)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
# --- PLOT B: FUTURE RICHNESS ---
plot(richnessFuture, 
     box = FALSE, axes = FALSE, col = my_cols, breaks = my_breaks,
     legend = FALSE, mar = c(1, 0, 1, 5), reset = TRUE)
plot(countries, add = TRUE)
mtext(expression(bold("b")), side = 3, adj = 0, line = 0.5)
fields::image.plot(
  zlim = c(0, 100), legend.only = TRUE, col = my_cols,
  smallplot = c(0.75, 0.77, 0.2, 0.8), 
  axis.args = list(at = seq(0, 100, 20), cex.axis = 0.8),
  legend.args = list(text = "Future richness\n(species)", side = 3, font = 1, line = 0.5, cex = 0.7)
)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
# --- COLOR PREP FOR C & D---
breaks_C <- c(seq(-50, 0, length.out = 51), seq(0.001, 25, length.out = 25))
cols_C   <- c(scico(50, palette = "roma", begin = 0, end = 0.45), # The Loss side
              scico(25, palette = "roma", begin = 0.55, end = 1)) # The Gain side
breaks_D <- c(seq(0, 1, length.out = 51), seq(1.001, 2, length.out = 50))
cols_D   <- c(scico(50, palette = "roma", begin = 0, end = 0.45),
              scico(50, palette = "roma", begin = 0.55, end = 1))
# --- PLOT C: ABSOLUTE CHANGE ---
plot(richnessAbs, 
     box = FALSE, axes = FALSE, col = cols_C, breaks = breaks_C,
     legend = FALSE, mar = c(1, 1, 1, 4), reset = TRUE)
plot(countries, add = TRUE)
mtext(expression(bold("c")), side = 3, adj = 0, line = 0.5)
fields::image.plot(
  zlim = c(-50, 25), legend.only = TRUE, col = cols_C,
  smallplot = c(0.80, 0.82, 0.2, 0.8), 
  axis.args = list(at = c(-50, -25, 0, 25), cex.axis = 0.8),
  legend.args = list(text = "Richness\nchange\n(species)", side = 3, font = 1, line = 0.5, cex = 0.7)
)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
# --- PLOT D: RELATIVE CHANGE ---
plot(richnessRel, 
     box = FALSE, axes = FALSE, col = cols_D, breaks = breaks_D,
     legend = FALSE, mar = c(1, 0, 1, 5), reset = TRUE)
plot(countries, add = TRUE)
mtext(expression(bold("d")), side = 3, adj = 0, line = 0.5)
fields::image.plot(
  zlim = c(0, 2), legend.only = TRUE, col = cols_D,
  smallplot = c(0.75, 0.77, 0.2, 0.8),
  axis.args = list(at = seq(0, 2, 0.5), cex.axis = 0.8),
  legend.args = list(text = "Richness\nchange\n(relative)", side = 3, font = 1, line = 0.5, cex = 0.7)
)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
dev.off()