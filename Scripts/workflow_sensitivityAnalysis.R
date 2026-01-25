library(terra)

# Set directory structure and load workflow functions
predPath <- "D:/Research/DroughtForecasts/Data/Predictors/"
spPath <- "D:/Research/DroughtForecasts/Data/Occurrences/BySpecies/Thinned/"
scriptPath <- "D:/Research/DroughtForecasts/Scripts/functions.R"
outPath <- "D:/Research/DroughtForecasts/OutputsSensitivityAnalysis/"
# Predictor rasters to fit models - all must have the same projection
predFilesCurrent <- list.files(predPath, "current.tif$", full.names = T)
# Predictor rasters for future projection - all must have the same projection - set to NULL to disable
predFilesFuture <- list.files(predPath, "SSP", full.names = T)
# Predictor raster with all possible variables, to ensure minimum sample size requirements are met
mostPreds <- rast("D:/Research/DroughtForecasts/Data/Predictors/ase_GFDL-ESM4_current.tif")

# Create directories if needed
if (!dir.exists(outPath)) dir.create(outPath)

# Set seed number
seedNo <- 2025

# Set number of cores for parallel processing - for predictor cropping, model fitting, and mapping 
numCores <- 32 

# Set terra options
terraOptions(memfrac = 0.9)

# Set modeling parameters
crs <- "+proj=longlat" # Projection for occurrence data
bgNo <- 20000 # Number of background points to be sampled
minNo <- 10 # Minimum number of occurrences for modeling
bufferSize <- c(100000, 500000) # Buffer size (m)
maxClusters <- 100 # Maximum number of clusters for spatial stratification
minClusters <- 10 # Maximum number of clusters for spatial stratification
finalFolds <- 5 # Final number of clusters after cluster size balancing
spatStratBootstraps <- 100 # Number of bootstraps for clustering
corrThreshold <- 0.7 # Threshold for removing correlated predictor variables
tune.args <- "lqp" # Set feature complexity
responseCurves <- F # Create response curves
varImportance <- 5 # Calculate variable importance with this number of permutations - set to 0 to disable
modelStats <- T # Calculate model statistics (AUC, BIC, Boyce index)
calcThresholds <- T # Calculate common thresholds (max. TSS, equal sens. and spec., 10% omission rate) for suitability maps
statsDistances <- NULL # Vector of distances (m) used for AUC and Boyce index subsetting - set to NULL to disable

# Set mapping parameters
mapping <- T # Create suitability and binary maps - calcThresholds needs to be set to T
aucThreshold <- 0 # Threshold for discarding models based on AUC - set to 0 to disable
boyceThreshold <- NA # Threshold for discarding models based on Boyce index - set to NA to disable
bgThreshold <- F # Discard models where the number of sampled background points is less than bgNo
clampSD <- c(0, 1) # Number of standard deviations for clamping future predictions (for both variables and features) - set to 0 only to disable clamping
cropSizes <- c(0, 100000, 500000) # Maximum dispersal sizes (m) - set to NULL to disable cropping

# Set cleaning options
rmCroppedPreds <- T # Removes cropped predictor files after modeling and mapping
rmSuitMaps <- F # Remove suitability maps after thresholding - plots are always saved

# Load modeling functions
source(scriptPath)

# List species to be modeled
spList <- list.files(spPath, "shp$", full.names = T)
spList <- spList[!grepl("envT_extT", spList, fixed = T)]

# Choose random species
species <- sample(unique(sapply(strsplit(basename(spList), "_"), `[`, 1)), 100)
write.csv(species, paste0(outPath, "speciesSample.csv"), row.names = F)
spList <- spList[grepl(paste(species, collapse = "|"), spList)]

# For debugging warnings
options(warn = 2)

# Fit models
for (currentSpecies in 1:length(spList)) {
  # Get species name
  sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
  sp <- sp[length(sp)]
  sp <- strsplit(sp, ".", fixed = T)[[1]][1]
  
  # Print species name
  print(paste0("Starting modeling workflow for species: ", sp))
  
  # Create output directories
  if (!dir.exists(paste0(outPath, sp))) {
    dir.create(paste0(outPath, sp))
    dir.create(paste0(outPath, sp, "/Background/"))
    dir.create(paste0(outPath, sp, "/Buffer/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/UncorrelatedPredictors/"))
    dir.create(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"))
    dir.create(paste0(outPath, sp, "/ENM/"))
    dir.create(paste0(outPath, sp, "/Presences/"))
    dir.create(paste0(outPath, sp, "/StatusReport/"))
    dir.create(paste0(outPath, sp, "/ModelSetComparison/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Current/"))
    dir.create(paste0(outPath, sp, "/CroppedPredictors/Future/"))
    dir.create(paste0(outPath, sp, "/SuitabilityMaps/Future/"))
    dir.create(paste0(outPath, sp, "/MarginalEffects/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Current/"))
    dir.create(paste0(outPath, sp, "/BinaryMaps/Future/"))
    dir.create(paste0(outPath, sp, "/RangeSizes/"))
    dir.create(paste0(outPath, sp, "/SpatialStratification/"))
  }
  
  # Match occurrence data projection to environmental data
  spOcc <- vect(spList[currentSpecies])
  spOcc <- project(spOcc, crs(rast(predFilesCurrent[1])))
  occProj <- crds(spOcc, df = T)
  names(occProj) <- c("lon", "lat")
  occProj$sp <- sp
  occProj$ID <- spOcc$AccNo
  
  # Prepare data frame to track status
  statusReport <- data.frame(Stage = character(), Status = character())
  
  # Check sample size
  print(paste0("Checking sample size requirements for species: ", sp))
  sampleStatus <- logical(length(predFilesCurrent))
  for (i in 1:length(predFilesCurrent)) {
    preds <- predFilesCurrent[i]
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
  # NOTE: this is proper behavior, as those points may include true presences.
  # NOTE: such points may not be caught by occTest, as occTest uses an unprojected environmental raster for extractions, 
  # NOTE: whereas in this workflow we use a projected environmental raster for extractions.
  print(paste0("Construct buffered convex hulls for species: ", sp))
  statusReport[3, 1] <- "Buffer constructed"
  statusReport[3, 2] <- 0
  for (i in 1:length(bufferSize)) {
    statusReport[3, 2] <- as.integer(statusReport[3, 2]) + buildBuffer(vect(spList[currentSpecies], crs = crs), bufferSize[i], predFilesCurrent[1], paste0(outPath, sp, "/Buffer/buffer_", format(bufferSize[i], scientific = F), ".shp"), plot = F)  
  }
  if (as.integer(statusReport[3, 2]) == length(bufferSize)) statusReport[3, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Crop predictor files using buffer
  print(paste0("Cropping predictor files for species: ", sp))
  buffers <- list.files(paste0(outPath, sp, "/Buffer/"), "shp$")
  success <- FALSE
  numCoresTmp <- numCores
  repeat {
    message(paste("Trying with", numCoresTmp, "cores."))
    if (numCoresTmp > 0) {
      cl <- makeCluster(numCoresTmp)
      registerDoParallel(cl)
    }
    tryCatch({
      foreach(i = 1:length(predFilesCurrent), .packages = c("terra")) %dopar% {
        preds <- predFilesCurrent[i]
        for (j in 1:length(buffers)) {
          cropPredictors(preds, vect(paste0(outPath, sp, "/Buffer/", buffers[j])), paste0(outPath, sp, "/CroppedPredictors/Current/"), buffers[j])
          gc()
        }
      }
      success <- TRUE
    }, error = function(e) {
      message("Error encountered: ", conditionMessage(e))
      if (grepl("std::bad_alloc", conditionMessage(e), fixed = TRUE) | grepl("paging file", conditionMessage(e), fixed = TRUE) | grepl("cannot allocate vector", conditionMessage(e), fixed = TRUE)) {
        if (numCoresTmp > 1) {
          numCoresTmp <<- max(1, floor(numCoresTmp / 2))
          message("Memory allocation failed. Retrying with ", numCoresTmp, " cores.")
        } else {
          stop("Failed even with 1 core — aborting.")
        }
      } else {
        stop(e)
      }
    }, finally = {
      if (exists("cl") && inherits(cl, "cluster")) stopCluster(cl)
    })
    if (success) {
      message("Processing completed successfully.")
      break
    }
  }
  statusReport[4, 1] <- "Cropped predictor files"
  statusReport[4, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Create background samples
  # NOTE: if fewer than the requested number of background points can be sampled (in cases of restricted ranges),
  # NOTE: then the function will proceed with a warning. The number of background points is always properly recorded.
  set.seed(seedNo)
  print(paste0("Sampling background points for species: ", sp))
  croppedPreds <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  statusReport[5, 1] <- "Minimum number of background points"
  for (i in 1:length(bufferSize)) {
    tmpCroppedPreds <- grepv(format(bufferSize[i], scientific = F), croppedPreds)[1]
    bgPoints <- sampleBackground(tmpCroppedPreds, bgNo)
    write.csv(bgPoints, paste0(outPath, sp, "/Background/background_buffer", format(bufferSize[i], scientific = F), ".csv"), row.names = F)
    bgSp <- bgPoints
    names(bgSp) <- c("lon", "lat")
    writeVector(vect(bgSp, crs = crs(spOcc)), paste0(outPath, sp, "/Background/background_buffer", format(bufferSize[i], scientific = F), ".shp"), overwrite = T)
    statusReport[5, 2] <- min(nrow(bgPoints), statusReport[5, 2], na.rm = T)
  }
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Remove correlated predictors
  print(paste0("Removing variables with correlations above threshold: ", corrThreshold))
  croppedPredFiles <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = T)
  croppedPredFilesShort <- list.files(paste0(outPath, sp, "/CroppedPredictors/Current/"), "tif$", full.names = F)
  predCountSumm <- data.frame(PredictorStack = character(),
                              RawCount = integer(),
                              FinalCount = integer(),
                              NoExcluded = integer())
  for (i in 1:length(croppedPredFiles)) {
    preds <- croppedPredFiles[i]
    predCountSumm[i, "PredictorStack"] <- croppedPredFilesShort[i]
    predCount <- nlyr(rast(croppedPredFiles[i]))
    predCountSumm[i, "RawCount"] <- predCount
    tossed <- removeCorrelatedPredictors(rast(preds), corrThreshold, paste0(outPath, sp, "/UncorrelatedPredictors/Current/", croppedPredFilesShort[i]), seedNo = seedNo)
    predCountSumm[i, "NoExcluded"] <- length(tossed)
    predCount <- predCount - length(tossed)
    predCountSumm[i, "FinalCount"] <- predCount
    write.csv(tossed, paste0(outPath, sp, "/UncorrelatedPredictors/Current/", croppedPredFilesShort[i], ".csv"), row.names = F)
  }
  write.csv(predCountSumm, paste0(outPath, sp, "/UncorrelatedPredictors/corrFilterSummary.csv"), row.names = F)
  statusReport[6, 1] <- "Removed correlated predictors"
  statusReport[6, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Spatially stratify occurrences
  print(paste0("Spatially stratifying occurrences: ", sp))
  statusReport[7, 1] <- "Number of strata (optimal, chosen)"
  strata <- spatialStratify(finalPres, seedNo, minClusters, maxClusters, spatStratBootstraps, finalFolds)
  statusReport[7, 2] <- paste0(strata$optimal_k, ", ", strata$chosen_k)
  save(strata, file = paste0(outPath, sp, "/SpatialStratification/spatialStratification.rda"))
  folds <- strata$focal.pres$folds
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Fit models, evaluate performance, and plot response curves
  print(paste0("Fitting models for species: ", sp))
  uncorrPredFiles <- list.files(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"), "tif$", full.names = T)
  uncorrPredFilesShort <- list.files(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"), "tif$", full.names = F)
  bgFiles <- list.files(paste0(outPath, sp, "/Background/"), "shp$", full.names = T)
  modelPerformance <- data.frame(model = character(),
                                 noPres = integer(),
                                 noBg = integer(),
                                 bgSize = character(),
                                 lambdaRule = character(),
                                 nvarraw = integer(),
                                 ncoefraw = integer(),
                                 varraw_drought = logical(),
                                 nvarfit = integer(),
                                 ncoeffit = integer(),
                                 varfit_drought = logical(),
                                 EPV = numeric(),
                                 AUC = numeric(),
                                 boyce = numeric(),
                                 BIC = numeric(),
                                 deltaBIC = numeric(),
                                 maxTSS = numeric(),
                                 eqsensspec = numeric(),
                                 OR10 = numeric(),
                                 sp = character(),
                                 varnamesfit = character(),
                                 varImportance = numeric())
  for (i in 1:length(uncorrPredFiles)) {
    print(uncorrPredFilesShort[i])
    modelPerformance[i, "model"] <- uncorrPredFilesShort[i]
    bgTmp <- strsplit(strsplit(uncorrPredFiles[i], "_buffer", fixed = T)[[1]][2], ".tif", fixed = T)[[1]][1]
    modelFit <- fitModel(pres = vect(paste0(outPath, sp, "/Presences/finalPres.shp")),
                         preds = rast(uncorrPredFiles[i]),
                         bg = vect(grepv(bgTmp, bgFiles)),
                         sp = sp,
                         tune.args = tune.args,
                         numCores = numCores,
                         folds = folds,
                         varImportance = varImportance,
                         modelStats = modelStats,
                         calcThresholds = calcThresholds,
                         statsDistances = statsDistances)
    if (inherits(modelFit, "try-error")) stop(paste0("Model fitting failed for taxon ", sp, " and predictor file ", uncorrPredFilesShort[i], " (index ", i, ")."))
    modelFit$sp <- sp
    modelFit$predFile <- uncorrPredFilesShort[i]
    modelFit$tune.args <- tune.args
    modelFit$folds <- folds
    modelFit$bgSize <- paste0(as.integer(bgTmp) / 1000, " km")
    modelFit$noPres <- length(vect(paste0(outPath, sp, "/Presences/finalPres.shp")))
    modelFit$noBg <- length(vect(grepv(bgTmp, bgFiles)))
    modelPerformance[i, "sp"] <- modelFit$sp
    modelPerformance[i, "noPres"] <- modelFit$noPres
    modelPerformance[i, "noBg"] <- modelFit$noBg
    modelPerformance[i, "bgSize"] <- modelFit$bgSize
    modelPerformance[i, "lambdaRule"] <- modelFit$lambdaRule
    modelPerformance[i, "nvarraw"] <- modelFit$nvarraw
    modelPerformance[i, "ncoefraw"] <- modelFit$ncoefraw
    modelPerformance[i, "varraw_drought"] <- modelFit$varraw_drought
    modelPerformance[i, "nvarfit"] <- modelFit$nvarfit
    modelPerformance[i, "varnamesfit"] <- paste0(modelFit$varnamesfit, collapse = " - ")
    modelPerformance[i, "varImportance"] <- ifelse(varImportance > 0, paste0(modelFit$varImportance, collapse = " - "), NA)
    modelPerformance[i, "ncoeffit"] <- modelFit$ncoeffit
    modelPerformance[i, "varfit_drought"] <- modelFit$varfit_drought
    modelPerformance[i, "EPV"] <- modelFit$EPV
    if (modelStats) {
      modelPerformance[i, "AUC"] <- modelFit$AUC
      modelPerformance[i, "boyce"] <- modelFit$boyce
      modelPerformance[i, "BIC"] <- ifelse(length(modelFit$cv_bic) == 0, NA, modelFit$cv_bic)
    } else {
      modelPerformance[i, "AUC"] <- NA
      modelPerformance[i, "boyce"] <- NA
      modelPerformance[i, "BIC"] <- NA
    }
    if (calcThresholds) {
      modelPerformance[i, "maxTSS"] <- modelFit$maxTSS
      modelPerformance[i, "eqsensspec"] <- modelFit$eqsensspec
      modelPerformance[i, "OR10"] <- modelFit$OR10
    } else {
      modelPerformance[i, "maxTSS"] <- NA
      modelPerformance[i, "eqsensspec"] <- NA
      modelPerformance[i, "OR10"] <- NA
    }
    # Plot error curve
    png(paste0(outPath, sp, "/ENM/", gsub(".tif", ".png", uncorrPredFilesShort[i], fixed = T)))
    plot(modelFit)
    if (!is.na(modelFit$lambda.value)) abline(v = -log(modelFit$lambda.value), col = "red", lwd = 2)
    dev.off()
    # Plot marginal effect curves
    # Note: response curve limits may be nonsensical (e.g. negative precipitation)
    if (responseCurves & length(modelFit$varnamesfit) > 0) {
      png(paste0(outPath, sp, "/MarginalEffects/", gsub(".tif", ".png", uncorrPredFilesShort[i], fixed = T)), width = 2000, height = 1500, res = 200)
      par(mfrow = c(ceiling(length(modelFit$varnamesfit) / 3), 3), mar = c(4, 4, 2, 1))
      for (j in 1:length(modelFit$varnamesfit)) {
        predFrame <- matrix(data = NA, nrow = 1000, ncol = length(modelFit$varnamesfit))
        colnames(predFrame) <- modelFit$varnamesfit
        seqMin <- modelFit$varmin[modelFit$varnamesfit[j]] - 2 * modelFit$varsd[modelFit$varnamesfit[j]]
        seqMax <- modelFit$varmax[modelFit$varnamesfit[j]] + 2 * modelFit$varsd[modelFit$varnamesfit[j]]
        seq <- seq(from = seqMin, to = seqMax, length.out = 1000)
        predFrame[, j] <- seq
        for (k in 1:length(modelFit$varnamesfit)) {
          if (k == j) next
          predFrame[, k] <- as.numeric(modelFit$samplemeans[modelFit$varnamesfit[k]])
        }
        predFrame <- as.data.frame(predFrame)
        predictions <- predict.cv.maxnet(object = modelFit, newdata = predFrame, type = "cloglog", clampVars = F, clampFeatures = F)
        predictionsClamped <- predict.cv.maxnet(object = modelFit, newdata = predFrame, type = "cloglog", clampVars = T, clampFeatures = T, clampRule = 1)
        plot(predictions ~ seq, type = "l", main = modelFit$varnamesfit[j], xlab = "Predictor", ylab = "Suitability")
        lines(x = seq, y = predictionsClamped, col = "red")
        abline(v = as.numeric(modelFit$samplemeans[modelFit$varnamesfit[j]]), col = "black", lwd = 2, lty = 2)
        abline(v = as.numeric(modelFit$samplemedians[modelFit$varnamesfit[j]]), col = "black", lwd = 2, lty = 3)
        abline(v = as.numeric(modelFit$varmin[modelFit$varnamesfit[j]]), col = "black", lwd = 1, lty = 2)
        abline(v = as.numeric(modelFit$varmax[modelFit$varnamesfit[j]]), col = "black", lwd = 1, lty = 2)
        abline(v = as.numeric(modelFit$varminobs[modelFit$varnamesfit[j]]), col = "green", lwd = 1, lty = 2)
        abline(v = as.numeric(modelFit$varmaxobs[modelFit$varnamesfit[j]]), col = "green", lwd = 1, lty = 2)
        abline(v = modelFit$varmin[modelFit$varnamesfit[j]] - 1 * modelFit$varsd[modelFit$varnamesfit[j]], col = "red", lwd = 1, lty = 2)
        abline(v = modelFit$varmax[modelFit$varnamesfit[j]] + 1 * modelFit$varsd[modelFit$varnamesfit[j]], col = "red", lwd = 1, lty = 2)
      }
      dev.off()
      par(mfrow = c(1, 1))
    }
    save(modelFit, file = paste0(outPath, sp, "/ENM/", gsub(".tif", ".rda", uncorrPredFilesShort[i], fixed = T)))
  }
  if (modelStats) {
    if (all(is.na(modelPerformance$BIC))) {
      modelPerformance$deltaBIC <- NA
    } else {
      modelPerformance$deltaBIC <- modelPerformance$BIC - min(modelPerformance$BIC, na.rm = TRUE)
    }
  }
  write.csv(modelPerformance, paste0(outPath, sp, "/ModelSetComparison/modelComparison.csv"), row.names = F)
  statusReport[8, 1] <- "Built ENMs"
  print(paste0("ENMs built for species: ", sp))
  statusReport[8, 2] <- "Passed"
  write.csv(statusReport, paste0(outPath, sp, "/StatusReport/statusReport.csv"), row.names = F)
  
  # Post-modeling cleaning
  if (rmCroppedPreds & !mapping) unlink(paste0(outPath, sp, "/CroppedPredictors/"), recursive = T)
  if (rmCroppedPreds) unlink(paste0(outPath, sp, "/UncorrelatedPredictors/Current/"), recursive = T)
  # Remove temporary files
  tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
  gc()
}

# Create suitability and binary maps
if (mapping & calcThresholds) {
  for (currentSpecies in 1:length(spList)) {
    # Get species name
    sp <- strsplit(spList[currentSpecies], "/", fixed = T)[[1]]
    sp <- sp[length(sp)]
    sp <- strsplit(sp, ".", fixed = T)[[1]][1]
    
    # Print species name
    print(paste0("Starting mapping workflow for species: ", sp))
    
    # Subset models
    print(paste0("Subsetting suitable models for species: ", sp))
    if (!file.exists(paste0(outPath, sp, "/ModelSetComparison/modelComparison.csv"))) {
      print(paste0("Mapping failed for taxon ", sp, ". No models found."))
      next
    } else {
      modelPerf <- read.csv(paste0(outPath, sp, "/ModelSetComparison/modelComparison.csv"), header = T)
      modelPerf <- subset(modelPerf, ncoeffit > 0)
      if (aucThreshold > 0) modelPerf <- subset(modelPerf, AUC > aucThreshold)
      if (!is.na(boyceThreshold)) modelPerf <- subset(modelPerf, boyce > boyceThreshold)
      if (bgThreshold) modelPerf <- subset(modelPerf, noBg == bgNo)
      if (nrow(modelPerf) == 0) {
        print(paste0("Mapping failed for taxon ", sp, ". No suitable models."))
        next
      }  
    }
    
    # Crop future predictors
    print(paste0("Cropping future predictor files for species: ", sp))
    if (!is.null(predFilesFuture)) {
      buffers <- list.files(paste0(outPath, sp, "/Buffer/"), "shp$")
      buffersTmp <- buffers[grep(
        paste(format(as.integer(unlist(strsplit(unique(modelPerf$bgSize), " km"))) * 1000, scientific = FALSE),
              collapse = "|"),
        buffers
      )]
      success <- FALSE
      numCoresTmp <- numCores
      repeat {
        message(paste("Trying with", numCoresTmp, "cores."))
        if (numCoresTmp > 0) {
          cl <- makeCluster(numCoresTmp)
          registerDoParallel(cl)
        }
        tryCatch({
          foreach(i = 1:length(predFilesFuture), .packages = c("terra")) %dopar% {
            preds <- predFilesFuture[i]
            for (j in 1:length(buffers)) {
              cropPredictors(
                preds,
                vect(paste0(outPath, sp, "/Buffer/", buffers[j])),
                paste0(outPath, sp, "/CroppedPredictors/Future/"),
                buffers[j]
              )
              gc()
            }
          }
          success <- TRUE
        }, error = function(e) {
          message("Error encountered: ", conditionMessage(e))
          if (grepl("std::bad_alloc", conditionMessage(e), fixed = TRUE) | grepl("paging file", conditionMessage(e), fixed = TRUE) | grepl("cannot allocate vector", conditionMessage(e), fixed = TRUE)
              | grepl("error writing to connection", conditionMessage(e), fixed = TRUE)) {
            if (numCoresTmp > 1) {
              numCoresTmp <<- max(1, floor(numCoresTmp / 2))
              message("Memory allocation failed. Retrying with ", numCoresTmp, " cores.")
            } else {
              stop("Failed even with 1 core — aborting.")
            }
          } else {
            stop(e)
          }
        }, finally = {
          if (exists("cl") && inherits(cl, "cluster")) stopCluster(cl)
        })
        if (success) {
          message("Processing completed successfully.")
          break
        }
      }
    }
    
    # Create maps
    rangeSizes <- data.frame()
    for (i in 1:nrow(modelPerf)) {
      # Load model
      load(gsub(".tif", ".rda", paste0(outPath, sp, "/ENM/", modelPerf[i, "model"]), fixed = T))
      print(paste0("Loaded model: ", gsub(".tif", ".rda", paste0(outPath, sp, "/ENM/", modelPerf[i, "model"]), fixed = T)))
      
      # Associate predictor maps
      print(paste0("Associating predictor maps for model: ", modelPerf[i, "model"]))
      predMaps <- list.files(paste0(outPath, sp, "/CroppedPredictors/"), "tif$", recursive = T)
      varSet <- paste0(strsplit(modelPerf[i, "model"], "_", fixed = T)[[1]][1], "_")
      bufferSet <- paste0(strsplit(modelPerf[i, "model"], "buffer", fixed = T)[[1]][2])
      predMaps <- grepv(varSet, predMaps)
      predMaps <- grepv(bufferSet, predMaps)
      if (varSet == "ase_") {
        GCMSet <- strsplit(modelPerf[i, "model"], "_", fixed = T)[[1]][2]
        predMaps <- grepv(GCMSet, predMaps)
      }
      
      # Predict suitability
      print(paste0("Creating suitability maps for model: ", modelPerf[i, "model"]))
      numCoresTmp <- numCores
      success <- FALSE
      repeat {
        message(paste("Trying with", numCoresTmp, "cores..."))
        if (numCoresTmp > 0) {
          cl <- makeCluster(numCoresTmp)
          registerDoParallel(cl)
        }
        tryCatch({
          results <- foreach(j = 1:length(predMaps), .packages = c("terra"), .combine = rbind) %dopar% {
            outList <- list()
            # Read in predictor values
            predVals <- as.data.frame(rast(paste0(outPath, sp, "/CroppedPredictors/", predMaps[j])), na.rm = FALSE)
            # Create template
            suitMap <- rast(paste0(outPath, sp, "/CroppedPredictors/", predMaps[j]))[[1]]
            names(suitMap) <- "suit"
            values(suitMap) <- NA
            # Make predictions
            for (k in 1:length(clampSD)) {
              if (clampSD[k] == 0) { 
                suitVals <- predict.cv.maxnet(
                  object = modelFit, newdata = predVals, type = "cloglog",
                  clampVars = FALSE, clampFeatures = FALSE
                )
              } else {
                suitVals <- predict.cv.maxnet(
                  object = modelFit, newdata = predVals, type = "cloglog",
                  clampVars = TRUE, clampFeatures = TRUE, clampRule = clampSD[k]
                )
              }
              values(suitMap) <- suitVals
              # Write suitability map
              outNameSuit <- paste0(
                outPath, sp, "/SuitabilityMaps/",
                sub("\\.tif$", paste0("_clamp", clampSD[k], ".tif"), predMaps[j])
              )
              writeRaster(suitMap * 10000, outNameSuit, overwrite = TRUE, datatype = "INT2U")
              png(gsub(".tif", ".png", outNameSuit, fixed = T))
              plot(suitMap, range = c(0, 1))
              plot(vect(paste0(outPath, sp, "/Presences/finalPres.shp")), add = TRUE)
              dev.off()
              # Threshold maps
              for (l in c("maxTSS", "eqsensspec", "OR10")) {
                tmpThreshold <- modelFit[l]
                binaryMap <- suitMap > as.numeric(tmpThreshold)
                binaryMap <- sum(binaryMap)
                if (is.null(cropSizes)) {
                  outNameBinary <- paste0(
                    outPath, sp, "/BinaryMaps/",
                    sub("\\.tif$", paste0("_clamp", clampSD[k], "_cropNA.tif"), predMaps[j])
                  )
                  outNameBinary <- sub("\\.tif$", paste0("_th-", l , ".tif"), outNameBinary)
                  writeRaster(binaryMap, outNameBinary, overwrite = TRUE, datatype = "INT1U")
                  png(gsub(".tif", ".png", outNameBinary, fixed = T))
                  plot(binaryMap, range = c(0, 1))
                  plot(vect(paste0(outPath, sp, "/Presences/finalPres.shp")), add = TRUE)
                  dev.off()
                  rangeSize <- freq(binaryMap, value = 1, bylayer = FALSE)[1, 2]
                  outList[[length(outList) + 1]] <- data.frame(basename(outNameBinary), rangeSize)
                } else {
                  cropSizesTmp <- cropSizes[
                    cropSizes <= as.integer(unlist(strsplit(modelFit$bgSize, " km"))) * 1000
                  ]
                  for (cropDist in cropSizesTmp) {
                    if (cropDist != as.integer(unlist(strsplit(modelFit$bgSize, " km"))) * 1000) {
                      chull <- convHull(vect(spList[currentSpecies]))
                      tmpBuffer <- if (cropDist > 0) buffer(chull, cropDist) else chull
                      tmpBuffer <- project(tmpBuffer, crs(rast(paste0(outPath, sp, "/CroppedPredictors/", predMaps[j]))))
                      croppedMap <- crop(binaryMap, tmpBuffer, mask = TRUE)
                    } else {
                      croppedMap <- binaryMap
                    }
                    outNameBinary <- paste0(
                      outPath, sp, "/BinaryMaps/",
                      sub("\\.tif$", paste0("_clamp", clampSD[k], "_crop", as.character(format(cropDist, scientific = FALSE)), ".tif"), predMaps[j])
                    )
                    outNameBinary <- sub("\\.tif$", paste0("_th-", l , ".tif"), outNameBinary)
                    
                    writeRaster(croppedMap, outNameBinary, overwrite = TRUE, datatype = "INT1U")
                    png(gsub(".tif", ".png", outNameBinary, fixed = T))
                    plot(croppedMap, range = c(0, 1))
                    plot(vect(paste0(outPath, sp, "/Presences/finalPres.shp")), add = TRUE)
                    dev.off()
                    rangeSize <- freq(croppedMap, value = 1, bylayer = FALSE)[1, 2]
                    outList[[length(outList) + 1]] <- data.frame(basename(outNameBinary), rangeSize)
                  }
                }
              }
            }
            do.call(rbind, outList)
          }
          success <- TRUE
        }, error = function(e) {
          message("Error encountered: ", conditionMessage(e))
          if (grepl("std::bad_alloc", conditionMessage(e), fixed = TRUE) | grepl("paging file", conditionMessage(e), fixed = TRUE) | grepl("cannot allocate vector", conditionMessage(e), fixed = TRUE) |
              grepl("error writing to connection", conditionMessage(e), fixed = TRUE)) {
            if (numCoresTmp > 1) {
              numCoresTmp <<- max(1, floor(numCoresTmp / 2))
              message("Memory allocation failed. Retrying with ", numCoresTmp, " cores.")
            } else {
              stop("Failed even with 1 core — aborting.")
            }
          } else {
            stop(e)
          }
        }, finally = {
          if (exists("cl") && inherits(cl, "cluster")) stopCluster(cl)
        })
        if (success) {
          message("Processing completed successfully.")
          break
        }
      }
      resultsAll <- cbind(results, modelPerf[rep(i, nrow(results)), ])
      rangeSizes <- rbind(rangeSizes, resultsAll)
    }
    # Export range size data
    write.csv(rangeSizes, paste0(outPath, sp, "/RangeSizes/rangeSizes.csv"), row.names = F)
    # Post-mapping cleaning
    if (rmSuitMaps) invisible(file.remove(list.files(paste0(outPath, sp, "/SuitabilityMaps/"), full.names = TRUE, pattern = "tif$", recursive = T)))
    if (rmCroppedPreds) unlink(paste0(outPath, sp, "/CroppedPredictors/"), recursive = T)
    # Remove temporary files
    tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    gc()
  }
}