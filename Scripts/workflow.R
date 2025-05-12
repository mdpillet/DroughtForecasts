library(terra)
library(remotes)
# Custom version of ENMeval package needs to be installed.
# install_github("mdpillet/ENMeval", force = T)
library(ENMeval)
library(rmarkdown)
library(plyr)
library(spThin)
library(sf)
library(pROC)
library(glmnet)

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
thinRadius <- 10 # Kilometers for spatial thinning
bgNo <- 10000 # Number of background points to be sampled
minNo <- 10 # Minimum number of occurrences for modeling
bufferSize <- 100000 # Buffer size
tune.args <- list(fc = c("LQP"), rm = 1) # Set feature complexity and regularization
nullModels <- F # Create null models for model selection
centile <- 0.1 # Quantile for thresholding based on omission rate
corrThreshold <- 0.7 # Threshold for removing correlated predictor variables
aucThreshold <- 0.5 # Threshold for discarding models based on cross-validated AUC
aucDistanceThreshold <- 0.7 # Threshold for discarding model based on full-model AUC with a background point distance filter of 50 km

# Set cleaning options
rmCroppedPreds <- T
rmTransformedPreds <- T

# Load modeling functions
source(scriptPath)

# List species to be modeled
spList <- list.files(spPath, "shp$", full.names = T)

# Loop workflow for all species
# NOTE: when the workflow throws an error, this is most likely due to model selection not succeeding.
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
    dir.create(paste0(outPath, sp, "/BufferNull/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/UncorrelatedPredictors/"))
    dir.create(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/ENM/"))
    dir.create(paste0(outPath, sp, "/Presences/"))
    dir.create(paste0(outPath, sp, "/StatusReport/"))
    if (nullModels) dir.create(paste0(outPath, sp, "/NullModels/"))
    dir.create(paste0(outPath, sp, "/ModelSetComparison/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Current/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Future/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Future/"))
    dir.create(paste0(outPath, sp, "/MarginalEffects/"))
    dir.create(paste0(outPath, sp, "/Thresholds/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Current/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Future/"))
    dir.create(paste0(outPath, sp, "/RangeSizes/"))
    dir.create(paste0(outPath, sp, "/ConcordanceMaps/"))
  }
  
  # Spatially thin occurrences and match occurrence data projection to environmental data
  spOcc <- vect(spList[currentSpecies])
  thinDF <- crds(spOcc, df = T)
  thinDF$sp <- sp
  thinnedOcc <- thin(thinDF, lat.col = "y", long.col = "x", spec.col = "sp", thin.par = thinRadius, reps = 1, 
                     locs.thinned.list.return = T, write.files = F, write.log.file = F, verbose = F)[[1]]
  spOcc <- spOcc[as.integer(row.names(thinnedOcc))]
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
  statusReport[3, 2] <- buildBuffer(vect(spList[currentSpecies], crs = crs), 0, predFiles[1], paste0(outPath, sp, "/BufferNull/buffer.shp"), plot = F)
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
  bgPoints <- sampleBackground(croppedPreds, bgNo)
  write.csv(bgPoints, paste0(outPath, sp, "/Background/background.csv"), row.names = F)
  bgSp <- bgPoints
  names(bgSp) <- c("lon", "lat")
  writeVector(vect(bgSp, crs = crs(spOcc)), paste0(outPath, sp, "/Background/background.shp"), overwrite = T)
  statusReport[5, 1] <- "Number of background points"
  statusReport[5, 2] <- nrow(bgPoints)
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Remove correlated predictors
  print(paste0("Removing variables with correlations above threshold: ", corrThreshold))
  croppedPredFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  croppedPredFilesShort <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = F)
  for (i in 1:length(croppedPredFiles)) {
    preds <- croppedPredFiles[i]
    tossed <- removeCorrelatedPredictors(rast(preds), corrThreshold, paste0(outPath, sp, "/UncorrelatedPredictors/Current/", croppedPredFilesShort[i]), seedNo = seedNo)
    write.csv(tossed, paste0(outPath, sp, "/UncorrelatedPredictors/Current/", croppedPredFilesShort[i], ".csv"), row.names = F)
  }
  statusReport[6, 1] <- "Removed correlated predictors"
  statusReport[6, 2] <- "Passed"
  
  # Fit models
  print(paste0("Fitting models for species: ", sp))
  uncorrPredFiles <- list.files(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"), "tif$", full.names = T)
  folds <- spatialStratify(finalPres, seedNo)$folds
  modelAlgs <- NULL
  for (i in 1:length(uncorrPredFiles)) {
    preds <- uncorrPredFiles[i]
    alg <- fitModel(finalPres, preds, bgPoints, sp, outDir = paste0(outPath, sp, "/ENM/"), tune.args, numCores, plot = F, folds)
    modelAlgs <- c(modelAlgs, alg)
  }
  statusReport[7, 1] <- "Built ENMs"
  print(paste0("ENMs built for species: ", sp))
  statusReport[7, 2] <- paste0(modelAlgs, collapse = ",")
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Compare model sets
  print(paste0("Comparing model sets for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  nullModelFiles <- list.files(paste0(outPath, sp, "/NullModels/"), "rda$", full.names = T)
  compareModelSets(modelFiles,
                   nullModels = nullModelFiles,
                   outDir = paste0(outPath, sp, "/ModelSetComparison/"))
  statusReport[8, 1] <- "Model set comparison"
  statusReport[8, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create suitability maps for current environments
  print(paste0("Creating current suitability maps for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  predFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  for (i in 1:length(modelFiles)) {
    tmp <- strsplit(modelFiles[i], "/", fixed = T)
    tmp <- tmp[[1]][length(tmp[[1]])]
    tmp <- strsplit(tmp, ".", fixed = T)
    tmp <- tmp[[1]][1]
    preds <- grep(tmp, predFiles, value = T)
    predictSuitability(modelFiles[i], preds, outDir = paste0(outPath, sp, "/SuitabilityMaps/Current/"))
  }
  statusReport[9, 1] <- "Create suitability maps for current environments"
  statusReport[9, 2] <- "Passed"
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
  statusReport[10, 1] <- "Cropped future predictor files"
  if (all(sampleStatus)) {
    print(paste0("Future predictor files cropped for species: ", sp))
    statusReport[10, 2] <- "Passed"
  } else {
    print(paste0("Future predictor files NOT cropped for species: ", sp))
    statusReport[10, 2] <- "Failed"
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Creating future suitability maps for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  predFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Future/"), "tif$", full.names = F)
  for (i in 1:length(modelFiles)) {
    tmp <- strsplit(modelFiles[i], "/", fixed = T)
    tmp <- tmp[[1]][length(tmp[[1]])]
    tmp <- strsplit(tmp, ".", fixed = T)
    tmp <- tmp[[1]][1]
    tmpVar <- strsplit(tmp, "_", fixed = T)[[1]][1]
    tmpGCM <- strsplit(tmp, "_", fixed = T)[[1]][2]
    tmpPreds <- grep(paste0(tmpVar, "_"), predFiles, fixed = T, value = T)
    if (tmpVar %in% c("a", "as", "ac", "asc")) {
      tmpPreds <- grep(paste0(tmpVar, "_"), predFiles, fixed = T, value = T)
    } else {
      tmpPreds <- grep(paste0(tmpVar, "_", tmpGCM), predFiles, fixed = T, value = T)
    }
    for (j in 1:length(tmpPreds)) {
      predictSuitability(modelFiles[i], paste0(outPath, sp, "/CroppedPredictors/Future/", tmpPreds[j]), outDir = paste0(outPath, sp, "/SuitabilityMaps/Future/"))
    }
  }
  statusReport[11, 1] <- "Create suitability maps for future environments"
  statusReport[11, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create marginal effect plots
  print(paste0("Extracting summary statistics for predictor files for species: ", sp))
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  croppedPredFilesCurrent <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  for (i in 1:length(croppedPredFilesCurrent)) {
    extractPredSummary(croppedPredFilesCurrent[i], paste0(outPath, sp, "/CroppedPredictors/Future/"), spOcc, paste0(outPath, sp, "/MarginalEffects/"))  
  }
  statusReport[12, 1] <- "Extracted summary statistics for predictor files"
  statusReport[12, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Create marginal effect plots for: ", sp))
  summFiles <- list.files(paste0(outPath, sp, "/MarginalEffects/"), "csv$", full.names = T)
  modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), "rda$", full.names = T)
  for (i in 1:length(summFiles)) {
    plotMarginal(summFiles[i], modelFiles[i], paste0(outPath, sp, "/MarginalEffects/"))
  }
  statusReport[13, 1] <- "Plotted marginal effects"
  statusReport[13, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Threshold maps
  print(paste0("Find thresholds for: ", sp))
  suitMaps <- list.files(paste0(outPath, sp, "/SuitabilityMaps/Current/"), "tif$", full.names = T)
  pres <- vect(list.files(paste0(outPath, sp, "/Presences/"), "shp$", full.names = T))
  getThresholds(suitMaps, pres, centile, outDir = paste0(outPath, sp, "/Thresholds/"))
  statusReport[14, 1] <- "Thresholds calculated"
  statusReport[14, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  print(paste0("Thresholding maps for: ", sp))
  suitMaps <- list.files(paste0(outPath, sp, "/SuitabilityMaps/"), "tif$", full.names = T, recursive = T)
  thresholds <- read.csv(paste0(outPath, sp, "/Thresholds/thresholds.csv"), header = T)
  threshold(suitMaps, thresholds, outDir = paste0(outPath, sp, "/RangeSizes/"))
  statusReport[15, 1] <- "Maps thresholded"
  statusReport[15, 2] <- "Passed"
  
  # Post-modeling cleaning
  if (rmCroppedPreds) unlink(paste0(outPath, sp, "/CroppedPredictors/"), recursive = T)
  if (rmTransformedPreds) unlink(paste0(outPath, sp, "/TransformedPredictors/"), recursive = T)
  
  # Remove temporary files
  tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}

# Calculate distance-based AUC by subsampling background points, only retaining points a minimum distance away from all presences
# NOTE: code may fail for large distances and will need to be restarted
distances <- c(0, 25000, 50000, 75000)
for (j in 1:length(distances)) {
  for (currentSpecies in 1:length(spList)) {
    # Get species name
    sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
    sp <- sp[length(sp)]
    sp <- strsplit(sp, ".", fixed = T)[[1]][1]
    print(sp)
    
    # Read in presences and absences
    if (file.exists(paste0(outPath, sp, "/Presences/finalPres.shp"))) {
      presences <- vect(paste0(outPath, sp, "/Presences/finalPres.shp"))
      background <- vect(paste0(outPath, sp, "/Background/background.shp"))  
    } else { next }
    
    # Calculate distances
    distance_matrix <- distance(background, presences)
    min_distances <- apply(distance_matrix, 1, min)
    background_filtered <- background[min_distances >= distances[j], ]
    
    if (file.exists(paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))) {
      models <- read.csv(paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"), header = T)
      # Find best model
      modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), full.names = T)
      for (i in 1:length(modelFiles)) {
        modelFile <- grep(models[i, "GCM"], modelFiles, value = T)
        # Read relevant predictor file
        buffer <- vect(paste0(outPath, sp, "/Buffer/buffer.shp"))
        predFiles <- list.files(predPath, full.names = T)
        predFiles <- grep(strsplit(strsplit(modelFile, "ase_")[[1]][2], "_current")[[1]][1], predFiles, value = T)
        predFilesCurrent <- predFiles[grepl("current", predFiles)]
        preds <- rast(predFilesCurrent)
        # Extract environment
        backgroundEnv <- terra::extract(preds, background_filtered, ID = F)
        presencesEnv <- terra::extract(preds, presences, ID = F)
        names(backgroundEnv) <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC",
                                  "HAD_mean", "HAD_median", "HAD_q0.75",
                                  "HAS_mean", "HAS_median", "HAS_q0.75",
                                  "HMD_mean", "HMD_median", "HMD_q0.75",
                                  "HMS_mean", "HMS_median", "HMS_q0.75")
        names(presencesEnv) <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC",
                                 "HAD_mean", "HAD_median", "HAD_q0.75",
                                 "HAS_mean", "HAS_median", "HAS_q0.75",
                                 "HMD_mean", "HMD_median", "HMD_q0.75",
                                 "HMS_mean", "HMS_median", "HMS_q0.75")
        # Predict suitability
        load(modelFile)
        suitValsPres <- predict(model@models$fc.LQP_rm.1, newdata = presencesEnv, clamp = F, type = "cloglog")
        suitValsPres <- as.data.frame(suitValsPres)
        suitValsPres$Presence <- 1
        suitValsBackground <- predict(model@models$fc.LQP_rm.1, newdata = backgroundEnv, clamp = F, type = "cloglog")
        suitValsBackground <- as.data.frame(suitValsBackground)
        suitValsBackground$Presence <- 0
        suitVals <- rbind(suitValsPres, suitValsBackground)
        names(suitVals)[1] <- "Suitability"
        # Calculate AUC
        roc <- roc(suitVals$Presence, suitVals$Suitability, quiet = T)
        models[i, paste0("AUC_", distances[j])] <- as.numeric(auc(roc))
      }
      write.csv(models, paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"), row.names = F)
    }
  }
}

# Create suitability and binary maps for best models and calculate range size statistics
for (currentSpecies in 1:length(spList)) {
  # Get species name
  sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
  sp <- sp[length(sp)]
  sp <- strsplit(sp, ".", fixed = T)[[1]][1]
  print(sp)
  
  # Delete previous runs
  unlink(paste0(outPath, sp, "/BinaryMapsBestModel/"), recursive = T)
  unlink(paste0(outPath, sp, "/BinaryMapsNoDispBestModel/"), recursive = T)
  unlink(paste0(outPath, sp, "/RangeSizesBestModel/"), recursive = T)
  
  # Get best model
  if (file.exists(paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"))) {
    models <- read.csv(paste0(outPath, sp, "/ModelSetComparison/modelSetComparison.csv"), header = T)
    models <- subset(models, AverageValidationAUC > aucThreshold)
    models <- subset(models, AUC_50000 > aucDistanceThreshold)
    if (nrow(models) == 0) {
      print("No models above AUC threshold.")
      next
    }
    else {
      # Find best model
      modelFiles <- list.files(paste0(outPath, sp, "/ENM/"), full.names = T)
      modelFiles <- grep(models[1, "GCM"], modelFiles, value = T)
      # Crop predictor files using buffer
      print(paste0("Cropping predictor files for species: ", sp))
      buffer <- vect(paste0(outPath, sp, "/Buffer/buffer.shp"))
      predFiles <- list.files(predPath, full.names = T)
      predFiles <- grep(strsplit(strsplit(modelFiles, "ase_")[[1]][2], "_current")[[1]][1], predFiles, value = T)
      predFilesCurrent <- predFiles[grepl("current", predFiles)]
      predFilesFuture <- predFiles[!grepl("current", predFiles)]
      dir.create(paste0(outPath, sp, "/CroppedPredictors/"))
      dir.create(paste0(outPath, sp, "/CroppedPredictors/Current/"))
      dir.create(paste0(outPath, sp, "/CroppedPredictors/Future/"))
      for (i in 1:length(predFilesCurrent)) {
        preds <- predFilesCurrent[i]
        cropPredictors(preds, buffer, paste0(outPath, sp, "/CroppedPredictors/Current/"))
      }
      for (i in 1:length(predFilesFuture)) {
        preds <- predFilesFuture[i]
        cropPredictors(preds, buffer, paste0(outPath, sp, "/CroppedPredictors/Future/"))
      }
      # Predict suitability
      dir.create(paste0(outPath, sp, "/SuitabilityMapsBestModel/"))
      dir.create(paste0(outPath, sp, "/SuitabilityMapsBestModel/Current/"))
      dir.create(paste0(outPath, sp, "/SuitabilityMapsBestModel/Future/"))
      croppedPreds <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current"), full.names = T)
      for (i in 1:length(croppedPreds)) {
        predictSuitability(modelFiles, croppedPreds[i], outDir = paste0(outPath, sp, "/SuitabilityMapsBestModel/Current/"))
      }
      croppedPreds <- list.files(paste0(outPath, sp, "/CroppedPredictors/Future"), full.names = T)
      for (i in 1:length(croppedPreds)) {
        predictSuitability(modelFiles, croppedPreds[i], outDir = paste0(outPath, sp, "/SuitabilityMapsBestModel/Future/"))
      }
      unlink(paste0(outPath, sp, "/CroppedPredictors/"), recursive = T)
      # Get threshold
      thresholds <- read.csv(paste0(outPath, sp, "/Thresholds/thresholds.csv"), header = T)
      threshold <- thresholds[grepl(models[1, "GCM"], thresholds$Model), "Threshold"]
      # Threshold suitability maps
      dir.create(paste0(outPath, sp, "/BinaryMapsBestModel/"))
      dir.create(paste0(outPath, sp, "/BinaryMapsBestModel/Current/"))
      dir.create(paste0(outPath, sp, "/BinaryMapsBestModel/Future/"))
      suitMaps <- list.files(paste0(outPath, sp, "/SuitabilityMapsBestModel/"), pattern = "tif$", full.names = T, include.dirs = F, recursive = T)
      for (i in 1:length(suitMaps)) {
        tmpMap <- rast(suitMaps[i])
        binaryMap <- tmpMap > threshold
        binaryMap <- sum(binaryMap)
        outName <- gsub("SuitabilityMapsBestModel", "BinaryMapsBestModel", suitMaps[i])
        writeRaster(binaryMap, outName, overwrite = T)
      }
      unlink(paste0(outPath, sp, "/SuitabilityMapsBestModel/"), recursive = T)
      # Calculate range sizes
      dir.create(paste0(outPath, sp, "/RangeSizesBestModel/"))
      rangeSizes <- data.frame(Map = character(),
                               RangeSize = numeric(),
                               RangeSizeChange = numeric())
      binaryMaps <- list.files(paste0(outPath, sp, "/BinaryMapsBestModel/"), pattern = "tif$", full.names = T, include.dirs = F, recursive = T)
      for (i in 1:length(binaryMaps)) {
        tmpMap <- rast(binaryMaps[i])
        size <- freq(tmpMap, digits = NA, value = 1, bylayer = F)[1, 2]
        rangeSizes[i, "Map"] <- binaryMaps[i]
        rangeSizes[i, "GCM"] <- strsplit(strsplit(binaryMaps[i], "ase_")[[1]][2], "_", fixed = T)[[1]][1]
        rangeSizes[i, "Time"] <- strsplit(strsplit(strsplit(binaryMaps[i], "ase_")[[1]][2], "_", fixed = T)[[1]][2], ".", fixed = T)[[1]][1]
        rangeSizes[i, "RangeSize"] <- size
        rangeSizes[i, "GeneratingModel"] <- strsplit(strsplit(modelFiles, "ase_")[[1]][2], "_current")[[1]][1]
      }
      GCMs <- unique(rangeSizes$GCM)
      for (i in 1:length(GCMs)) {
        currentSize <- subset(rangeSizes, GCM == GCMs[i] & Time == "current")$RangeSize
        for (j in 1:nrow(rangeSizes)) {
          if (rangeSizes[j, "GCM"] == GCMs[i] & rangeSizes[j, "Time"] != "current") {
            rangeSizes[j, "RangeSizeChange"] <- rangeSizes[j, "RangeSize"] / currentSize
          }
        }
      }
      write.csv(rangeSizes, paste0(outPath, sp, "/RangeSizesBestModel/rangeSizes.csv"), row.names = F)
      # Read zero-dispersal buffer
      buffer0 <- vect(paste0(outPath, sp, "/BufferNull/buffer.shp"))
      # Crop binary maps
      dir.create(paste0(outPath, sp, "/BinaryMapsNoDispBestModel/"))
      dir.create(paste0(outPath, sp, "/BinaryMapsNoDispBestModel/Current/"))
      dir.create(paste0(outPath, sp, "/BinaryMapsNoDispBestModel/Future/"))
      binaryMaps <- list.files(paste0(outPath, sp, "/BinaryMapsBestModel/"), pattern = "tif$", full.names = T, include.dirs = F, recursive = T)
      for (i in 1:length(binaryMaps)) {
        tmpMap <- rast(binaryMaps[i])
        croppedMap <- crop(tmpMap, buffer0)
        outName <- gsub("BinaryMapsBestModel", "BinaryMapsNoDispBestModel", binaryMaps[i])
        writeRaster(croppedMap, outName, overwrite = T)
      }
      # Calculate range sizes without dispersal
      rangeSizes <- data.frame(Map = character(),
                               RangeSize = numeric(),
                               RangeSizeChange = numeric())
      binaryMaps <- list.files(paste0(outPath, sp, "/BinaryMapsNoDispBestModel/"), pattern = "tif$", full.names = T, include.dirs = F, recursive = T)
      for (i in 1:length(binaryMaps)) {
        tmpMap <- rast(binaryMaps[i])
        size <- freq(tmpMap, digits = NA, value = 1, bylayer = F)[1, 2]
        rangeSizes[i, "Map"] <- binaryMaps[i]
        rangeSizes[i, "GCM"] <- strsplit(strsplit(binaryMaps[i], "ase_")[[1]][2], "_", fixed = T)[[1]][1]
        rangeSizes[i, "Time"] <- strsplit(strsplit(strsplit(binaryMaps[i], "ase_")[[1]][2], "_", fixed = T)[[1]][2], ".", fixed = T)[[1]][1]
        rangeSizes[i, "RangeSize"] <- size
        rangeSizes[i, "GeneratingModel"] <- strsplit(strsplit(modelFiles, "ase_")[[1]][2], "_current")[[1]][1]
      }
      GCMs <- unique(rangeSizes$GCM)
      for (i in 1:length(GCMs)) {
        currentSize <- subset(rangeSizes, GCM == GCMs[i] & Time == "current")$RangeSize
        for (j in 1:nrow(rangeSizes)) {
          if (rangeSizes[j, "GCM"] == GCMs[i] & rangeSizes[j, "Time"] != "current") {
            rangeSizes[j, "RangeSizeChange"] <- rangeSizes[j, "RangeSize"] / currentSize
          }
        }
      }
      write.csv(rangeSizes, paste0(outPath, sp, "/RangeSizesBestModel/rangeSizesNoDisp.csv"), row.names = F)
    }
  }
  else {
    print("No models found.")
    next
  }
}