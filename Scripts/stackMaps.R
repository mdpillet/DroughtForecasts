library(terra)

# Set directory structure and load workflow functions
spPath <- "D:/Research/DroughtForecasts/Outputs/"
outPath <- "D:/Research/DroughtForecasts/Data/DiversityMaps/"
countryPath <- "D:/Research/DroughtForecasts/Data/CountriesOutline/world-administrative-boundaries.shp"
figPath <- "D:/Research/DroughtForecasts/Manuscript/Figures/"

# List species to be modeled
files <- list.files(path = "D:/Research/DroughtForecasts/Outputs/",
                    pattern = "tif$", 
                    recursive = TRUE, 
                    full.names = TRUE)

# Filter for files in subdirectories matching "BinaryMapsBestModel"
files100 <- files[grepl("BinaryMapsBestModel", files)]
files0 <- files[grepl("BinaryMapsNoDispBestModel", files)]
maps <- c("current", "SSP126", "SSP370", "SSP585")

# Create diversity maps
for (i in 1:length(maps)) {
  tmp100 <- files100[grepl(maps[i], files100)]
  print(i)
  tmp <- lapply(tmp100, rast)
  tmp <- sprc(tmp)
  stacked <- mosaic(tmp, fun = "sum")
  writeRaster(stacked, paste0(outPath, "100km_", maps[i], ".tif"), overwrite = T)
  png(filename = paste0(outPath, "100km_", maps[i], ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(stacked, range = c(0, 100))
  dev.off()
}

for (i in 1:length(maps)) {
  tmp0 <- files0[grepl(maps[i], files0)]
  print(i)
  tmp <- lapply(tmp0, rast)
  tmp <- sprc(tmp)
  stacked <- mosaic(tmp, fun = "sum")
  writeRaster(stacked, paste0(outPath, "0km_", maps[i], ".tif"), overwrite = T)
  png(filename = paste0(outPath, "0km_", maps[i], ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(stacked, range = c(0, 100))
  dev.off()
}

# Create publication figures
richnessCurrent <- rast(paste0(outPath, "100km_current.tif"))
richnessFuture <- rast(paste0(outPath, "100km_SSP585.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "Fig4.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()

richnessCurrent <- rast(paste0(outPath, "100km_current.tif"))
richnessFuture <- rast(paste0(outPath, "100km_SSP370.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "ExtendedDataFig2.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()

richnessCurrent <- rast(paste0(outPath, "100km_current.tif"))
richnessFuture <- rast(paste0(outPath, "100km_SSP126.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "ExtendedDataFig3.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()

richnessCurrent <- rast(paste0(outPath, "0km_current.tif"))
richnessFuture <- rast(paste0(outPath, "0km_SSP585.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "ExtendedDataFig4.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()

richnessCurrent <- rast(paste0(outPath, "0km_current.tif"))
richnessFuture <- rast(paste0(outPath, "0km_SSP370.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "ExtendedDataFig5.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()

richnessCurrent <- rast(paste0(outPath, "0km_current.tif"))
richnessFuture <- rast(paste0(outPath, "0km_SSP126.tif"))
richnessAbs <- richnessFuture - richnessCurrent 
richnessRel <- richnessFuture / richnessCurrent
countries <- vect(countryPath)
countries <- project(countries, richnessCurrent)

png(paste0(figPath, "ExtendedDataFig6.png"), units = "mm", width = 180, height = 180, res = 300)
options(terra.pal = map.pal("viridis", 100))
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))
plotA <- plot(richnessCurrent, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("a")), side = 3, adj = 0)
plotB <- plot(richnessFuture, box = F, axes = F, range = c(0, 100), plg = list(title = "Species", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("b")), side = 3, adj = 0)
options(terra.pal = hcl.colors(100, palette = "Red-Green"))
plotC <- plot(richnessAbs, box = F, axes = F, range = c(-50, 50), plg = list(title = "Richness\nchange\n(species)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("c")), side = 3, adj = 0)
plotD <- plot(richnessRel, box = F, axes = F, range = c(0, 2), plg = list(title = "Richness\nchange\n(relative)", title.cex = 1, cex = 1, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
plot(countries, add = T)
mtext(expression(bold("d")), side = 3, adj = 0)
dev.off()