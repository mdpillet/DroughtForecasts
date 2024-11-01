library(terra)
library(occTest)
library(ggpubr)

# Set paths
occPath <- "D:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Over10/"
envRaster <- "D:/Research/DroughtForecasts/Data/Predictors/ase_UKESM1-0-LL_SSP585.tif"
outPath <- "D:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Filtered/" 

# Read occurrence data
species <- list.files(occPath, "shp$", full.names = T)

# Read environmental raster
env <- rast(envRaster)
envProj <- project(env, crs(vect(species[1])))
names(envProj) <- c("bio1", paste0("bio", 10:19), paste0("bio", 2:9), "cec", "cfvo", "clay", "pH", "sand", "silt", "SOC", "HAD_mean", "HAD_median", "HAD_q0.75",
                          "HAS_mean", "HAS_median", "HAS_q0.75",
                          "HMD_mean", "HMD_median", "HMD_q0.75",
                          "HMS_mean", "HMS_median", "HMS_q0.75")

# Set settings
customSettings <- defaultSettings()
customSettings$analysisSettings$filterAtlas <- F
# ncores <- 7

# Note that currently occTest fails for:
# i = 418 (Eriosyce jussieui) - acceptable failure: all points in single grid cell
# Loop over species
for (i in 1:length(species)) {
  occ <- vect(species[i])
  print(unique(occ$FinalSpeci))
  coords <- as.data.frame(crds(occ))
  names(coords) <- c("decimalLongitude", "decimalLatitude") 
  occTest <- occTest(sp.name = unique(occ$FinalSpeci),
              habitat = "terrestrial",
              sp.table = coords,
              r.env = envProj, 
              interactiveMode = F,
              verbose = F,
              analysisSettings = customSettings$analysisSettings,
              doParallel = F,
              mc.cores = ncores)
  # Strict filtering
  occFilter_strict <- occFilter(df = occTest, errorAcceptance = "strict")
  list_of_plots <- plot(x = occTest, occFilter_list = occFilter_strict, show_plot = F)
  print(paste0("Retained/total: ", length(occFilter_strict$filteredDataset$taxonobservationID), "/", nrow(occ)))
  if (length(occFilter_strict$filteredDataset$taxonobservationID) >= 10) {
    ggarrange(list_of_plots[[1]], list_of_plots[[2]], list_of_plots[[3]], list_of_plots[[4]])
    ggsave(paste0(outPath, "Over10/", unique(occ$FinalSpeci), ".jpg"))
    writeVector(occ[occFilter_strict$filteredDataset$taxonobservationID], paste0(outPath, "Over10/", unique(occ$FinalSpeci), ".kml"), overwrite = T)
    writeVector(occ[occFilter_strict$filteredDataset$taxonobservationID], paste0(outPath, "Over10/", unique(occ$FinalSpeci), ".shp"), overwrite = T)
    save(occTest, occFilter_strict, file = paste0(outPath, "Over10/", unique(occ$FinalSpeci), ".rda"))
  }
  else {
    ggarrange(list_of_plots[[1]], list_of_plots[[2]], list_of_plots[[3]], list_of_plots[[4]])
    ggsave(paste0(outPath, "Other/", unique(occ$FinalSpeci), ".jpg"))
    writeVector(occ[occFilter_strict$filteredDataset$taxonobservationID], paste0(outPath, "Other/", unique(occ$FinalSpeci), ".kml"), overwrite = T)
    writeVector(occ[occFilter_strict$filteredDataset$taxonobservationID], paste0(outPath, "Other/", unique(occ$FinalSpeci), ".shp"), overwrite = T)
    save(occTest, occFilter_strict, file = paste0(outPath, "Other/", unique(occ$FinalSpeci), ".rda"))
  }
  # Majority filtering
  # occFilter_maj <- occFilter(df = occTest, errorAcceptance = "majority")
  # list_of_plots <- plot(x = occTest, occFilter_list = occFilter_maj, show_plot = F)
  # ggarrange(list_of_plots[[1]], list_of_plots[[2]], list_of_plots[[3]], list_of_plots[[4]])
  # ggsave(paste0(outPath, unique(occ$FinalSpeci), "_maj.jpg"))
  # writeVector(occ[occFilter_maj$filteredDataset$taxonobservationID], paste0(outPath, unique(occ$FinalSpeci), "_maj.kml"), overwrite = T)
  # writeVector(occ[occFilter_maj$filteredDataset$taxonobservationID], paste0(outPath, unique(occ$FinalSpeci), "_maj.shp"), overwrite = T)
  # Relaxed filtering
  # occFilter_relax <- occFilter(df = occTest, errorAcceptance = "relaxed")
  # list_of_plots <- plot(x = occTest, occFilter_list = occFilter_relax, show_plot = F)
  # ggarrange(list_of_plots[[1]], list_of_plots[[2]], list_of_plots[[3]], list_of_plots[[4]])
  # ggsave(paste0(outPath, unique(occ$FinalSpeci), "_relax.jpg"))
  # writeVector(occ[occFilter_relax$filteredDataset$taxonobservationID], paste0(outPath, unique(occ$FinalSpeci), "_relax.kml"), overwrite = T)
  # writeVector(occ[occFilter_relax$filteredDataset$taxonobservationID], paste0(outPath, unique(occ$FinalSpeci), "_relax.shp"), overwrite = T)
  # Save objects
  # save(occTest, occFilter_strict, occFilter_maj, occFilter_relax, file = paste0(outPath, unique(occ$FinalSpeci), ".rda"))
  closeAllConnections()
}