library(terra)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
summaryPath <- "Data/DiversityMaps/Summary/"
figPath <- "Manuscript/Figures/"
expertPath <- "Data/ExpertMaps/"
mapPath <- "Data/DiversityMaps/"
countryPath <- "Data/CountriesOutline/world-administrative-boundaries.shp"

# Stack expert maps
tmp <- lapply(list.files(paste0(relPath, expertPath), "shp$", full.names = T), vect)
tmp <- do.call(rbind, tmp)
tmp$Richness <- 1
tmp <- project(tmp, crs(rast(paste0(relPath, summaryPath, "a_current_mean.tif"))))

# Calculate diversity
tmpRast <- rasterize(tmp, rast(paste0(relPath, summaryPath, "a_current_mean.tif")), field = "Richness", fun = "sum")

# Calculate correlations with expert map
richnessA <- rast(paste0(relPath, mapPath, "Current/a_current_current.tif"))
richnessB <- rast(paste0(relPath, mapPath, "Current/as_current_current.tif"))
richnessC <- rast(paste0(relPath, summaryPath, "ae_current_mean.tif"))
richnessD <- rast(paste0(relPath, summaryPath, "ase_current_mean.tif"))
cor.test(values(richnessA), values(tmpRast))
cor.test(values(richnessB), values(tmpRast))
cor.test(values(richnessC), values(tmpRast))
cor.test(values(richnessD), values(tmpRast))

# Create manuscript figure
countries <- vect(paste0(relPath, countryPath))
countries <- project(countries, crs(rast(paste0(relPath, summaryPath, "a_current_mean.tif"))))
png(paste0(relPath, figPath, "ExtendedDataFig2.png"), units = "mm", width = 90, height = 90, res = 300)
plot <- plot(tmpRast, box = F, axes = F, range = c(0, 80), plg = list(title = "Species", title.cex = 0.5, cex = 0.5, size = c(0.75, 0.5)), mar = c(1.5, 1.5, 1.5, 1.5), buffer = F, cex.main = 0.75)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.5)
plot(countries, add = T)
dev.off()