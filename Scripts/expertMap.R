library(terra)
library(scico)
library(fields)

# Set directory structure
relPath <- "D:/Research/DroughtForecasts/"
figPath <- "Manuscript/Submission2/Figures/"
expertPath <- "Data/ExpertMaps/"
mapPath <- "Data/DiversityMaps/Current/"
countryPath <- "Data/CountriesOutline/world-administrative-boundaries.shp"

# Stack expert maps
tmp <- lapply(list.files(paste0(relPath, expertPath), "shp$", full.names = T), vect)
tmp <- do.call(rbind, tmp)
tmp$Richness <- 1
tmp <- project(tmp, crs(rast(paste0(relPath, mapPath, "avg.tif"))))

# Calculate diversity
tmpRast <- rasterize(tmp, rast(paste0(relPath, mapPath, "avg.tif")), field = "Richness", fun = "sum")

# Calculate correlations with expert map
richness <- rast(paste0(relPath, mapPath, "avg.tif"))
cor.test(values(richness), values(tmpRast))

# Create manuscript figure
my_cols <- c(
  "#001959", "#031D5A", "#06215B", "#07245A", "#0A285B", "#0A2C5C", 
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
countries <- vect(paste0(relPath, countryPath))
countries <- project(countries, crs(rast(paste0(relPath, mapPath, "avg.tif"))))
png(paste0(relPath, figPath, "ExtendedDataFig7.png"), units = "mm", width = 180, height = 180, res = 300)
plot(tmpRast, 
     box = FALSE, axes = FALSE, col = c("white", rev(my_cols)), breaks = my_breaks,
     legend = FALSE, mar = c(1, 1, 1, 4), reset = TRUE)
plot(countries, add = TRUE)
fields::image.plot(
  zlim = c(0, 100), legend.only = TRUE, col = c("white", rev(my_cols)),
  smallplot = c(0.85, 0.87, 0.2, 0.8), 
  axis.args = list(at = seq(0, 100, 20), cex.axis = 0.8),
  legend.args = list(text = "Richness\n(species)", side = 3, font = 1, line = 0.5, cex = 0.7)
)
sbar(1000000, "bottomleft", type = "bar", labels = c("", "1,000 km", ""), cex = 0.75)
dev.off()