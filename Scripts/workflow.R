library(terra)
library(ENMeval)
library(rmarkdown)
library(plyr)

# Set directory structure and load workflow functions
predPath <- "D:/Research/DroughtForecasts/Data/Predictors/"
spPath <- "D:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Filtered/Over10/"
scriptPath <- "D:/Research/DroughtForecasts/Scripts/functions.R"
outPath <- "D:/Research/DroughtForecasts/Outputs/"

# Create directories if needed
if (!dir.exists(outPath)) dir.create(outPath)

# Set seed number
seedNo <- 2024

# Set number of cores for modeling
numCores <- 7

# Set modeling parameters
crs <- "+proj=longlat" # Projection for occurrence data
bgNo <- 10000 # Number of background points to be sampled
minNo <- 10 # Minimum number of occurrences for modeling
bufferSize <- 100000 # Buffer size
# tune.args <- list(fc = c("LQ", "H", "LQH"), rm = seq(1, 5, 0.5)) # Set feature complexity and regularization
tune.args <- list(fc = c("LQ"), rm = seq(1, 5, 0.5)) # Set feature complexity and regularization
nullModels <- F # Create null models for model selection
centile <- 0.1 # Quantile for thresholding based on omission rate

# Set cleaning options
rmCroppedPreds <- T
rmTransformedPreds <- T

# Load modeling functions
source(scriptPath)

# List species to be modeled
spList <- list.files(spPath, "shp$", full.names = T)

# Loop workflow for all species
# NOTE: when the workflow throws an error, this is most likely due to model selection not succeeding when certain evaluation statistics (e.g. AICc) cannot be calculated.
# NOTE: when this occurs, the workflow should be restarted after the species in question.
for (currentSpecies in 1:length(spList)) {
  # Get species name
  sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
  sp <- sp[length(sp)]
  sp <- strsplit(sp, ".", fixed = T)[[1]][1]
  
  # Print species name
  print(paste0("Starting workflow for species: ", sp))
  
  # Create output directories
  if (!dir.exists(paste0(outPath, sp))) {
    dir.create(paste0(outPath, sp))
    dir.create(paste0(outPath, sp, "/Background/"))
    dir.create(paste0(outPath, sp, "/Buffer/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/ENM/"))
    dir.create(paste0(outPath, sp, "/ModelsPCA/"))
    dir.create(paste0(outPath, sp, "/Presences/"))
    dir.create(paste0(outPath, sp, "/TransformedPredictors/"))
    dir.create(paste0(outPath, sp, "/TransformedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/StatusReport/"))
    dir.create(paste0(outPath, sp, "/BestModels/"))
    if (nullModels) dir.create(paste0(outPath, sp, "/NullModels/"))
    dir.create(paste0(outPath, sp, "/ModelSetComparison/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Current/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Future/"))
    dir.create(paste0(outPath, sp, "/TransformedPredictors/Future/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Future/"))
    dir.create(paste0(outPath, sp, "/MarginalEffects/"))
    dir.create(paste0(outPath, sp, "/Thresholds/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Current/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Future/"))
    dir.create(paste0(outPath, sp, "/RangeSizes/"))
    dir.create(paste0(outPath, sp, "/ConcordanceMaps/"))
  }
  
  # Match occurrence data projection to environmental data
  spOcc <- vect(spList[currentSpecies])
  predFiles <- list.files(predPath, "current.tif$", full.names = T)
  preds <- rast(predFiles[1])
  spOcc <- project(spOcc, crs(preds))
  occProj <- crds(spOcc, df = T)
  names(occProj) <- c("lon", "lat")
  occProj$sp <- sp
  occProj$ID <- spOcc$AccNo
  
  # Prepare data frame to track status
  statusReport <- data.frame(Stage = character(), Status = character())
  
  # Check sample size
  print(paste0("Checking sample size requirements for species: ", sp))
  sampleStatus <- logical(length(predFiles))
  for (i in 1:length(predFiles)) {
    preds <- predFiles[i]
    sampleStatus[i] <- checkSampleSize(occProj, preds, n = minNo)
  }
  statusReport[1, 1] <- "Sample size requirements"
  if (all(sampleStatus)) {
    print(paste0("Sample size requirements met for species: ", sp))
    statusReport[1, 2] <- "Passed"
    write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  } else {
    print(paste0("Sample size requirements NOT met for species: ", sp))
    statusReport[1, 2] <- "Failed"
    write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
    next
  }
  
  # Get coordinates for final presence records
  print(paste0("Retrieving final presence records for species: ", sp))
  mostPreds <- list.files(predPath, "ase.*tif$", full.names = T)[1]
  finalPres <- getPresenceRecords(occProj, mostPreds)
  write.csv(finalPres, paste0(outPath, sp, "/Presences/finalPres.csv"), row.names = F)
  finalPresSp <- finalPres
  names(finalPresSp) <- c("lon", "lat")
  writeVector(vect(finalPresSp, crs = crs(spOcc)), paste0(outPath, sp, "/Presences/finalPres.shp"), overwrite = T)
  statusReport[2, 1] <- "Final number of presences"
  statusReport[2, 2] <- nrow(finalPres)
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Build buffer around unprojected presences
  # NOTE: the buffer is built around all points that passed occTest, even those for which environmental variables can't be extracted.
  # NOTE: this is proper behavior, as those points are likely true presences.
  # NOTE: such points may not be caught by occTest, as occTest uses an unprojected environmental raster for extractions, 
  # NOTE: whereas in this workflow we use a projected environmental raster for extractions.
  print(paste0("Construct buffered convex hull for species: ", sp))
  statusReport[3, 1] <- "Buffer constructed"
  statusReport[3, 2] <- buildBuffer(vect(spList[currentSpecies], crs = crs), bufferSize, predFiles[1], paste0(outPath, sp, "/Buffer/buffer.shp"), plot = F)
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Crop predictor files using buffer
  print(paste0("Cropping predictor files for species: ", sp))
  buffer <- vect(paste0(outPath, sp, "/Buffer/buffer.shp"))
  sampleStatus <- logical(length(predFiles))
  for (i in 1:length(predFiles)) {
    preds <- predFiles[i]
    sampleStatus[i] <- cropPredictors(preds, buffer, paste0(outPath, sp, "/CroppedPredictors/Current/"))
  }
  statusReport[4, 1] <- "Cropped predictor files"
  if (all(sampleStatus)) {
    print(paste0("Predictor files cropped for species: ", sp))
    statusReport[4, 2] <- "Passed"
  } else {
    print(paste0("Predictor files NOT cropped for species: ", sp))
    statusReport[4, 2] <- "Failed"
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create background samples
  # NOTE: if fewer than the requested number of background points can be sampled (in cases of restricted ranges),
  # NOTE: then the function will proceed with a warning. The number of background points is always properly recorded.
  set.seed(seedNo)
  print(paste0("Sampling background points for species: ", sp))
  croppedPreds <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "ase.*tif$", full.names = T)[1]
  bgPoints <- sampleBackground(croppedPreds, occProj, bgNo)
  write.csv(bgPoints, paste0(outPath, sp, "/Background/background.csv"), row.names = F)
  bgSp <- bgPoints
  names(bgSp) <- c("lon", "lat")
  writeVector(vect(bgSp, crs = crs(spOcc)), paste0(outPath, sp, "/Background/background.shp"), overwrite = T)
  statusReport[5, 1] <- "Number of background points"
  statusReport[5, 2] <- nrow(bgPoints)
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Fit PCA models
  print(paste0("Fitting PCA models for species: ", sp))
  croppedPredFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  modelStatus <- integer(length(croppedPredFiles))
  for (i in 1:length(croppedPredFiles)) {
    preds <- croppedPredFiles[i]
    modelStatus[i] <- fitPCA(preds, vect(bgSp, crs = crs(spOcc)), paste0(outPath, sp, "/ModelsPCA/"), seedNo = seedNo)
  }
  statusReport[6, 1] <- "Minimum no. of PCs"
  statusReport[7, 1] <- "Maximum no. of PCs"
  statusReport[6, 2] <- min(modelStatus)
  statusReport[7, 2] <- max(modelStatus)
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Transform cropped predictors using PCA models
  print(paste0("Transforming cropped predictors using PCA models for species: ", sp))
  sampleStatus <- logical(length(croppedPredFiles))
  for (i in 1:length(croppedPredFiles)) {
    preds <- croppedPredFiles[i]
    sampleStatus[i] <- transformPCA(preds, paste0(outPath, sp, "/ModelsPCA/"), paste0(outPath, sp, "/TransformedPredictors/Current/"))
  }
  statusReport[8, 1] <- "Transformed predictor files"
  if (all(sampleStatus)) {
    print(paste0("Predictor files transformed for species: ", sp))
    statusReport[8, 2] <- "Passed"
  } else {
    print(paste0("Predictor files NOT transformed for species: ", sp))
    statusReport[8, 2] <- "Failed"
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Fit models
  print(paste0("Fitting models for species: ", sp))
  transPredFiles <- list.files(paste0(outPath, sp, "/TransformedPredictors/Current/"), "tif$", full.names = T)
  for (i in 1:length(transPredFiles)) {
    preds <- transPredFiles[i]
    fitModel(finalPres, preds, bgPoints, sp, outDir = paste0(outPath, sp, "/ENM/"), tune.args, numCores, plot = T)
  }
  statusReport[9, 1] <- "Built ENMs"
  print(paste0("ENMs built for species: ", sp))
  statusReport[9, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Select best model
  print(paste0("Selecting best models for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  for (i in 1:length(modelFiles)) {
    print(paste0("Model selection for set ", i))
    finalModel <- selectModel(modelFiles[i], 
                              outDirModel = paste0(outPath, sp, "/BestModels/"),
                              numCores = numCores,
                              nullModels = nullModels,
                              outDirNulls = paste0(outPath, sp, "/NullModels/"))
  }
  statusReport[10, 1] <- "Selected best ENMs"
  print(paste0("Best ENMs selected for species: ", sp))
  statusReport[10, 2] <- "Passed"
  statusReport[11, 1] <- "Null models for model selection"
  statusReport[11, 2] <- nullModels
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Compare model sets
  print(paste0("Comparing model sets for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/BestModels/"), "rda$", full.names = T)
  nullModelFiles <- list.files(paste0(outPath, sp, "/NullModels/"), "rda$", full.names = T)
  compareModelSets(modelFiles,
                   nullModels = nullModelFiles,
                   outDir = paste0(outPath, sp, "/ModelSetComparison/"))
  statusReport[12, 1] <- "Model set comparison"
  statusReport[12, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create suitability maps for current environments
  print(paste0("Creating current suitability maps for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/BestModels/"), "rda$", full.names = T)
  predFilesTrans <- list.files(paste0(outPath, sp, "/TransformedPredictors/Current/"), "tif$", full.names = T)
  for (i in 1:length(modelFiles)) {
    tmp <- strsplit(modelFiles[i], "/", fixed = T)
    tmp <- tmp[[1]][length(tmp[[1]])]
    tmp <- strsplit(tmp, ".", fixed = T)
    tmp <- tmp[[1]][1]
    preds <- grep(tmp, predFilesTrans, value = T)
    predictSuitability(modelFiles[i], preds, outDir = paste0(outPath, sp, "/SuitabilityMaps/Current/"))
  }
  statusReport[13, 1] <- "Create suitability maps for current environments"
  statusReport[13, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create suitability maps for future environments
  print(paste0("Cropping future predictor files for species: ", sp))
  predFiles <- list.files(predPath, "SSP....tif$", full.names = T)
  buffer <- vect(paste0(outPath, sp, "/Buffer/buffer.shp"))
  sampleStatus <- logical(length(predFiles))
  for (i in 1:length(predFiles)) {
    preds <- predFiles[i]
    sampleStatus[i] <- cropPredictors(preds, buffer, paste0(outPath, sp, "/CroppedPredictors/Future/"))
  }
  statusReport[14, 1] <- "Cropped future predictor files"
  if (all(sampleStatus)) {
    print(paste0("Future predictor files cropped for species: ", sp))
    statusReport[14, 2] <- "Passed"
  } else {
    print(paste0("Future predictor files NOT cropped for species: ", sp))
    statusReport[14, 2] <- "Failed"
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Transforming cropped future predictors using PCA models for species: ", sp))
  croppedPredFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Future/"), "tif$", full.names = T)
  sampleStatus <- logical(length(croppedPredFiles))
  for (i in 1:length(croppedPredFiles)) {
    preds <- croppedPredFiles[i]
    tmp <- strsplit(croppedPredFiles[i], "/", fixed = T)
    tmp <- tmp[[1]][length(tmp[[1]])]
    tmp <- strsplit(tmp, ".", fixed = T)
    tmp <- tmp[[1]][1]
    tmpVar <- strsplit(tmp, "_", fixed = T)[[1]][1]
    tmpGCM <- strsplit(tmp, "_", fixed = T)[[1]][2]
    if (tmpVar %in% c("a", "as", "ac", "asc")) {
      sampleStatus[i] <- transferPCA(preds, paste0(outPath, sp, "/ModelsPCA/", tmpVar, "_current_current.rda"), paste0(outPath, sp, "/TransformedPredictors/Future/"))  
    } else {
      sampleStatus[i] <- transferPCA(preds, paste0(outPath, sp, "/ModelsPCA/", tmpVar, "_", tmpGCM, "_current.rda"), paste0(outPath, sp, "/TransformedPredictors/Future/"))
    }
  }
  statusReport[15, 1] <- "Transformed future predictor files"
  if (all(sampleStatus)) {
    print(paste0("Future predictor files transformed for species: ", sp))
    statusReport[15, 2] <- "Passed"
  } else {
    print(paste0("Future predictor files NOT transformed for species: ", sp))
    statusReport[15, 2] <- "Failed"
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Creating future suitability maps for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/BestModels/"), "rda$", full.names = T)
  predFiles <- list.files(paste0(outPath, sp, "/TransformedPredictors/Future/"), "tif$", full.names = F)
  for (i in 1:length(modelFiles)) {
    tmp <- strsplit(modelFiles[i], "/", fixed = T)
    tmp <- tmp[[1]][length(tmp[[1]])]
    tmp <- strsplit(tmp, ".", fixed = T)
    tmp <- tmp[[1]][1]
    tmpVar <- strsplit(tmp, "_", fixed = T)[[1]][1]
    tmpGCM <- strsplit(tmp, "_", fixed = T)[[1]][2]
    if (tmpVar %in% c("a", "as", "ac", "asc")) {
      tmpPreds <- grep(paste0(tmpVar, "_"), predFiles, fixed = T, value = T)
    } else {
      tmpPreds <- grep(paste0(tmpVar, "_", tmpGCM), predFiles, fixed = T, value = T)
    }
    for (j in 1:length(tmpPreds)) {
      predictSuitability(modelFiles[i], paste0(outPath, sp, "/TransformedPredictors/Future/", tmpPreds[j]), outDir = paste0(outPath, sp, "/SuitabilityMaps/Future/"))
    }
  }
  statusReport[16, 1] <- "Create suitability maps for future environments"
  statusReport[16, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create concordance maps
  print(paste0("Create concordance maps for: ", sp))
  suitMapsCur <- list.files(paste0(outPath, sp, "/SuitabilityMaps/Current/"), "tif$", full.names = T)
  suitMapsFut <- list.files(paste0(outPath, sp, "/SuitabilityMaps/Future/"), "tif$", full.names = T)
  ## Calculate model weights
  suitMapsCur_noE <- grep("/a_|/ac_|/as_|/asc_", suitMapsCur, value = T)
  suitMapsFut_noE <- grep("/a_|/ac_|/as_|/asc_", suitMapsFut, value = T)
  modelWeightsCur_noE <- calculateWeights(suitMapsCur_noE, paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))
  modelWeightsCur_all <- calculateWeights(suitMapsCur, paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))
  modelWeightsFut_noE <- calculateWeights(suitMapsFut_noE, paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))
  modelWeightsFut_all <- calculateWeights(suitMapsFut, paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))
  ## Current - no extreme droughts
  tmp1 <- app(rast(suitMapsCur_noE), weighted.mean, w = modelWeightsCur_noE$weight, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/noE_current.tif"), overwrite = T)
  sd1 <- app(rast(suitMapsCur_noE), sd, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/noE_current_sd.tif"), overwrite = T)
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/noE_current.png"), width = 10, height = 10, units = "in", res = 320)
  plot(tmp1, range = c(0, 1))
  dev.off()
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/noE_current_sd.png"), width = 10, height = 10, units = "in", res = 320)
  plot(sd1, range = c(0, 0.5))
  dev.off()
  ## Current - all
  tmp2 <- app(rast(suitMapsCur), weighted.mean, w = modelWeightsCur_all$weight, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/all_current.tif"), overwrite = T)
  sd2 <- app(rast(suitMapsCur), sd, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/all_current_sd.tif"), overwrite = T)
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all_current.png"), width = 10, height = 10, units = "in", res = 320)
  plot(tmp2, range = c(0, 1))
  dev.off()
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all_current_sd.png"), width = 10, height = 10, units = "in", res = 320)
  plot(sd2, range = c(0, 0.5))
  dev.off()
  ## Future - no extreme droughts
  tmp3 <- app(rast(suitMapsFut_noE), weighted.mean, w = modelWeightsFut_noE$weight, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/noE_future.tif"), overwrite = T)
  sd3 <- app(rast(suitMapsFut_noE), sd, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/noE_future_sd.tif"), overwrite = T)
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/noE_future.png"), width = 10, height = 10, units = "in", res = 320)
  plot(tmp3, range = c(0, 1))
  dev.off()
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/noE_future_sd.png"), width = 10, height = 10, units = "in", res = 320)
  plot(sd3, range = c(0, 0.5))
  dev.off()
  ## Future - all
  tmp4 <- app(rast(suitMapsFut), weighted.mean, w = modelWeightsFut_all$weight, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/all_future.tif"), overwrite = T)
  sd4 <- app(rast(suitMapsFut), sd, na.rm = T, filename = paste0(outPath, sp, "/ConcordanceMaps/all_future_sd.tif"), overwrite = T)
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all_future.png"), width = 10, height = 10, units = "in", res = 320)
  plot(tmp4, range = c(0, 1))
  dev.off()
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all_future_sd.png"), width = 10, height = 10, units = "in", res = 320)
  plot(sd4, range = c(0, 0.5))
  dev.off()
  ## Create comparison plot
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all.png"), width = 10, height = 10, units = "in", res = 320)
  par(mfrow = c(2, 2))
  plot(tmp1, range = c(0, 1), main = "No drought variables - historical", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(tmp2, range = c(0, 1), main = "All variables - historical", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(tmp3, range = c(0, 1), main = "No drought variables - 2071-2100", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(tmp4, range = c(0, 1), main = "All variables - 2071-2100", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  dev.off()
  ## Create SD plot
  png(filename = paste0(outPath, sp, "/ConcordanceMaps/all_sd.png"), width = 10, height = 10, units = "in", res = 320)
  par(mfrow = c(2, 2))
  upperLim <- round_any(max(minmax(sd1)[2], minmax(sd2)[2], minmax(sd3)[2], minmax(sd4)[2]), 0.25, ceiling)
  plot(sd1, range = c(0, upperLim), main = "No drought variables - historical", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(sd2, range = c(0, upperLim), main = "All variables - historical", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(sd3, range = c(0, upperLim), main = "No drought variables - 2071-2100", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  plot(sd4, range = c(0, upperLim), main = "All variables - 2071-2100", axes = F)
  points(spOcc, cex = 0.25, col = "red")
  dev.off()
  
  # Create marginal effect plots
  print(paste0("Extracting summary statistics for predictor files for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/BestModels/"), "rda$", full.names = T)
  croppedPredFilesCurrent <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  for (i in 1:length(croppedPredFilesCurrent)) {
    extractPredSummary(croppedPredFilesCurrent[i], paste0(outPath, sp, "/CroppedPredictors/Future/"), spOcc, paste0(outPath, sp, "/MarginalEffects/"))  
  }
  statusReport[17, 1] <- "Extracted summary statistics for predictor files"
  statusReport[17, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Create marginal effect plots for: ", sp))
  summFiles <- list.files(paste0(outPath, sp, "/MarginalEffects/"), "csv$", full.names = T)
  pcaFiles <- list.files(paste0(outPath, sp, "/ModelsPCA/"), "rda$", full.names = T)
  modelFiles <- list.files(paste0(outPath, sp, "/BestModels/"), "rda$", full.names = T)
  for (i in 1:length(summFiles)) {
    plotMarginal(summFiles[i], pcaFiles[i], modelFiles[i], paste0(outPath, sp, "/MarginalEffects/"))
  }
  statusReport[18, 1] <- "Plotted marginal effects"
  statusReport[18, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Threshold maps
  print(paste0("Find thresholds for: ", sp))
  suitMaps <- list.files(paste0(outPath, sp, "/SuitabilityMaps/Current/"), "tif$", full.names = T)
  pres <- vect(list.files(paste0(outPath, sp, "/Presences/"), "shp$", full.names = T))
  getThresholds(suitMaps, pres, centile, outDir = paste0(outPath, sp, "/Thresholds/"))
  statusReport[19, 1] <- "Thresholds calculated"
  statusReport[19, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Thresholding maps for: ", sp))
  suitMaps <- list.files(paste0(outPath, sp, "/SuitabilityMaps/"), "tif$", full.names = T, recursive = T)
  thresholds <- read.csv(paste0(outPath, sp, "/Thresholds/thresholds.csv"), header = T)
  threshold(suitMaps, thresholds, outDir = paste0(outPath, sp, "/RangeSizes/"))
  statusReport[20, 1] <- "Maps thresholded"
  statusReport[20, 2] <- "Passed"
  
  # Post-modeling cleaning
  if (rmCroppedPreds) unlink(paste0(outPath, sp, "/CroppedPredictors/"), recursive = T)
  if (rmTransformedPreds) unlink(paste0(outPath, sp, "/TransformedPredictors/"), recursive = T)
  
  # Remove temporary files
  tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}