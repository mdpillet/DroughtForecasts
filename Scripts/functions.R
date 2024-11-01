library(terra)
library(dplyr)
library(paran)
library(ggplot2)

# This function checks if a minimum of n presence points are available for model fitting
# 1. At least n records must be present
# 2. At least n unique grid cells must be present
# 3. No NAs may be present in the environmental layers for those grid cells
# 4. At least n unique combinations of predictors must be present

checkSampleSize <- function(occ, preds, n) {
  if (nrow(occ) >= n) {
    preds <- rast(preds)
    spOcc <- vect(occ, crs = crs(preds))
    extr <- terra::extract(preds, spOcc, cells = T, ID = F)
    extr <- extr[complete.cases(extr),]
    if (length(unique(extr$cell)) >= n) {
      usable <- nrow(distinct(extr[, 1:nlyr(preds)]))
      if (usable >= n) return(TRUE)
    }
  }
  return(FALSE)
}

# This function extracts all cells with presences that have unique predictor values
# NOTE: mostPreds should represent an environmental layers file with a complete set of predictors as to build equivalent models when comparing variable selection methods
# NOTE: Spatial thinning is implicitly performed by only retaining a single point per grid cell
getPresenceRecords <- function(occ, mostPreds) {
  preds <- rast(mostPreds)
  spOcc <- vect(occ, crs = crs(preds))
  extr <- terra::extract(preds, spOcc, cells = F, xy = T, ID = F)
  extr <- extr[complete.cases(extr),]
  keepRows <- unique.data.frame(extr[, 1:nlyr(preds)])
  extr <- extr[row.names(keepRows), c("x", "y")]
  return(extr)
}

# This function builds a buffered convex hull around a set of occurrences
# NOTE: occurrence data should be unprojected
buildBuffer <- function(occ, distance, preds, outName, plot = F) {
  chull <- convHull(occ)
  buffer <- buffer(chull, distance)
  preds <- rast(preds)
  buffer <- project(buffer, crs(preds))
  writeVector(buffer, outName, overwrite = T)
  if (plot) {
    plot(buffer)
    points(project(occ, crs(preds)))
  }
  return(TRUE)
}

# This function crops a predictor file based on a buffer
cropPredictors <- function(preds, buffer, outDir) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  layerTemplate <- strsplit(outName, "_", fixed = T)[[1]][1]
  if (layerTemplate == "a") layerNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9))
  if (layerTemplate == "as") layerNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC")
  if (layerTemplate == "ae") layerNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), 
                                             "HAD_mean", "HAD_median", "HAD_q0.75",
                                             "HAS_mean", "HAS_median", "HAS_q0.75",
                                             "HMD_mean", "HMD_median", "HMD_q0.75",
                                             "HMS_mean", "HMS_median", "HMS_q0.75")
  if (layerTemplate == "ase") layerNames <- c("BIO1", paste0("BIO", 10:19), paste0("BIO", 2:9), "CEC", "CFVO", "clay", "pH", "sand", "silt", "SOC",
                                              "HAD_mean", "HAD_median", "HAD_q0.75",
                                              "HAS_mean", "HAS_median", "HAS_q0.75",
                                              "HMD_mean", "HMD_median", "HMD_q0.75",
                                              "HMS_mean", "HMS_median", "HMS_q0.75")
  preds <- rast(preds)
  croppedPreds <- crop(preds, buffer, mask = T)
  names(croppedPreds) <- layerNames
  writeRaster(croppedPreds, paste0(outDir, outName), overwrite = T)
  return(TRUE)
}

# This function samples background points from an environmental layers file
sampleBackground <- function(preds, pres, no) {
  preds <- raster::stack(preds)
  bgSample <- dismo::randomPoints(preds, no, pres, lonlatCorrection = F)
  return(bgSample)
}

# This function fits a PCA model on a set of background points. Parallel analysis is performed.
# NOTE: zero-variance variables are removed
fitPCA <- function(preds, bgPoints, outDir, seedNo) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  preds <- rast(preds)
  extr <- terra::extract(preds, bgPoints, cells = F, xy = F, ID = F)
  extr <- extr[, which(apply(extr, 2, var) != 0)]
  model <- prcomp(extr, center = T, scale. = T)
  parHorn <- paran(extr, centile = 95, seed = seedNo, iterations = 1000, quietly = T, status = F)$Retained
  save(model, parHorn, file = paste0(outDir, outName, ".rda"))
  return(parHorn)
}

# This function transforms predictors using a PCA model
# NOTE: a minimum of 2 PCs is always retained
transformPCA <- function(preds, modelDir, outDir) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  preds <- rast(preds)
  predsValues <- values(preds)
  load(paste0(modelDir, outName, ".rda"))
  trans <- predict(model, predsValues)
  predsTrans <- preds[[1:ncol(trans)]]
  values(predsTrans) <- trans
  predsTrans <- predsTrans[[1:max(parHorn, 2)]]
  names(predsTrans) <- paste0("PC", 1:max(parHorn, 2))
  writeRaster(predsTrans, paste0(outDir, outName, ".tif"), overwrite = T)
  return(TRUE)
}

# This function fits Maxent models
fitModel <- function(pres, preds, bg, sp, outDir, tune.args, numCores, plot) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  model <- ENMevaluate(occs = pres,
                       envs = preds,
                       bg = bg,
                       tune.args = tune.args,
                       algorithm = "maxnet",
                       partitions = "block",
                       doClamp = F,
                       taxon.name = sp,
                       overlap = F,
                       parallel = T,
                       numCores = numCores,
                       quiet = T)
  save(model, file = paste0(outDir, outName, ".rda"))
  if (plot) {
    modelStats <- evalplot.stats(model, stats = c("or.10p", "auc.val"), x.var = "rm", color = "fc", dodge = 0.5)
    ggsave(paste0(outDir, outName, ".png"), modelStats)
  }
  return(TRUE)
}

# This function performs model selection by selecting models with the lowest average test omission rate, breaking ties with the highest average validation AUC,
# then lowest AICc, then by highest regularization multiplier, then by preferring H > LQ > LQH.
# Null models are (optionally) created.
selectModel <- function(modelSet, outDirModel, numCores, nullModels, outDirNulls) {
  # Select best model
  outName <- strsplit(modelSet, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  load(modelSet)
  res <- eval.results(model)
  opt.seq <- res %>% 
    filter(or.10p.avg == min(or.10p.avg)) %>% 
    filter(auc.val.avg == max(auc.val.avg)) %>%
    filter(AICc == min(AICc)) %>%
    filter(as.numeric(rm) == max(as.numeric(rm))) %>%
    filter(nchar(as.character(fc)) == min(nchar(as.character(fc))))
  bestModel <- eval.models(model)[[opt.seq$tune.args]]
  save(bestModel, opt.seq, file = paste0(outDirModel, outName, ".rda"))
  # Save metadata
  rmm <- eval.rmm(model)
  rmm$model$selectionRules <- "Lowest 10 percentile omission rate, break ties with average validation AUC and then AICc and then regularization and then feature complexity."
  rmm$model$finalModelSettings <- paste0(as.character(opt.seq$fc), as.character(opt.seq$rm))
  rangeModelMetadata::rmmToCSV(rmm, paste0(outDirModel, outName, ".csv"))
  # Create null models
  if (nullModels) {
    mod.null <- ENMnulls(model, mod.settings = list(fc = as.character(opt.seq$fc), rm = as.numeric(opt.seq$rm)), no.iter = 100,
                         eval.stats = c("auc.val", "or.10p"),
                         parallel = F,
                         numCores = numCores,
                         quiet = T)
    summ.table <- null.emp.results(mod.null)
    save(mod.null, summ.table, file = paste0(outDirNulls, outName, ".rda"))
    nullPlots <- evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")
    ggsave(paste0(outDirNulls, outName, ".png"), nullPlots)
  }
  return(TRUE)
}

# This function compares model sets based on AICc.
compareModelSets <- function(modelSets, nullModels, outDir) {
  nrow <- length(modelSets)
  comp <- data.frame(VariableSet = character(),
                     GCM = character(),
                     FeatureComplexity = character(),
                     RegularizationMultiplier = numeric(),
                     AICc = numeric(),
                     DeltaAICc = numeric(),
                     AverageValidationAUC = numeric(),
                     SignificanceAUC = numeric(),
                     AverageOmissionRate = numeric(),
                     SignificanceOmissionRate = numeric())
  for (i in 1:nrow) {
    tmp <- modelSets[i]
    outName <- strsplit(tmp, "/", fixed = T)
    outName <- outName[[1]][length(outName[[1]])]
    outName <- strsplit(outName, ".", fixed = T)
    outName <- outName[[1]][1]
    comp[i, "VariableSet"] <- toupper(strsplit(outName, "_", fixed = T)[[1]][1])
    if (strsplit(outName, "_", fixed = T)[[1]][2] == "current") {
      comp[i, "GCM"] <- NA
    } else {
      comp[i, "GCM"] <- strsplit(outName, "_", fixed = T)[[1]][2]
    }
    load(tmp)
    comp[i, "FeatureComplexity"] <- as.character(opt.seq$fc)
    comp[i, "RegularizationMultiplier"] <- as.numeric(levels(opt.seq$rm))[opt.seq$rm]
    comp[i, "AICc"] <- opt.seq$AICc
    comp[i, "AverageValidationAUC"] <- opt.seq$auc.val.avg
    comp[i, "AverageOmissionRate"] <- opt.seq$or.10p.avg
    if (length(nullModels) == nrow) {
      load(nullModelFiles[i])
      comp[i, "SignificanceAUC"] <- summ.table[summ.table$statistic == "pvalue", "auc.val"]
      comp[i, "SignificanceOmissionRate"] <- summ.table[summ.table$statistic == "pvalue", "or.10p"]
    }
  }
  comp$DeltaAICc <- comp$AICc - min(comp$AICc)
  comp <- arrange(comp, DeltaAICc)
  write.csv(comp, paste0(outDir, "modelSetComparison.csv"), row.names = F)
}

# This function creates a suitability raster
predictSuitability <- function(model, preds, outDir) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  load(model)
  preds <- raster::stack(preds)
  suit <- enm.maxnet@predict(bestModel, envs = preds, other.settings = list(doClamp = F, pred.type = "cloglog"))
  suit <- rast(suit)
  writeRaster(suit, paste0(outDir, outName, ".tif"), overwrite = T)
  png(filename = paste0(outDir, outName, ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(suit, range = c(0, 1))
  dev.off()
}

# This function transforms predictors using a PCA model for future environments
# NOTE: a minimum of 2 PCs is always retained
transferPCA <- function(preds, PCAmodel, outDir) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  preds <- rast(preds)
  predsValues <- values(preds)
  load(PCAmodel)
  trans <- predict(model, predsValues)
  predsTrans <- preds[[1:ncol(trans)]]
  values(predsTrans) <- trans
  predsTrans <- predsTrans[[1:max(parHorn, 2)]]
  names(predsTrans) <- paste0("PC", 1:max(parHorn, 2))
  writeRaster(predsTrans, paste0(outDir, outName, ".tif"), overwrite = T)
  return(TRUE)
}

# This function extracts summaries from predictors for plotting
extractPredSummary <- function(predsCurrent, predsFutureDir, occ, outDir) {
  outName <- strsplit(predsCurrent, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  tmpRast <- rast(predsCurrent)
  valRast <- values(tmpRast)
  summStats <- matrix(nrow = 5, ncol = nlyr(tmpRast))
  colnames(summStats) <- names(tmpRast)
  summStats[1,] <- apply(valRast, 2, min, na.rm = T)
  summStats[3,] <- apply(valRast, 2, mean, na.rm = T)
  presExtr <- terra::extract(tmpRast, occ, ID = F)
  summStats[2,] <- apply(presExtr, 2, mean, na.rm = T)
  predFilesFuture <- list.files(predsFutureDir, "tif$", full.names = T)
  tmpVar <- strsplit(outName, "_", fixed = T)[[1]][1]
  tmpGCM <- strsplit(outName, "_", fixed = T)[[1]][2]
  if (tmpVar %in% c("a", "as", "ac", "asc")) {
    tmpPreds <- grep(paste0(tmpVar, "_"), predFilesFuture, fixed = T, value = T)
  } else {
    tmpPreds <- grep(paste0(tmpVar, "_", tmpGCM), predFilesFuture, fixed = T, value = T)
  }
  for (i in 1:length(tmpPreds)) {
    tmpFut <- rast(tmpPreds[i])
    apply(values(tmpFut), 2, min, na.rm = T)
    summStats[4,] <- pmin(apply(values(tmpFut), 2, min, na.rm = T), summStats[4,], na.rm = T)
    summStats[5,] <- pmin(apply(values(tmpFut), 2, max, na.rm = T), summStats[5,], na.rm = T)
  }
  rownames(summStats) <- c("MinPresent", "MaxPresent", "MeanPresent", "MinFuture", "MaxFuture")
  summStats <- as.data.frame(summStats)
  write.csv(summStats, paste0(outDir, outName, ".csv"))
}

# This function creates marginal effect plots
plotMarginal <- function(summFile, pcaFile, modelFile, outDir) {
  outName <- strsplit(summFile, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  summStats <- read.csv(summFile, header = T, row.names = 1)
  load(pcaFile)
  load(modelFile)
  png(filename = paste0(outDir, outName, ".png"), width = 20, height = 20, units = "in", res = 320)
  par(mfrow = c(ceiling(ncol(summStats) / 4), 4))
  for (i in 1:ncol(summStats)) {
    tmpSeq <- matrix(nrow = 1000, ncol = ncol(summStats))
    for (j in 1:ncol(tmpSeq)) {
      if (j == i) {
        tmpSeq[, j] <- seq(min(summStats[1, j], summStats[4, j]), max(summStats[2, j], summStats[5, j]), length.out = 1000)
      } else {
        tmpSeq[, j] <- summStats[3, j]
      }
    }
    colnames(tmpSeq) <- colnames(summStats)
    tmpSeq <- as.data.frame(tmpSeq)
    tmpSeqTrans <- predict(model, tmpSeq)
    tmpSeqTrans <- tmpSeqTrans[, 1:max(parHorn, 2)]
    suit <- enm.maxnet@predict(bestModel, envs = tmpSeqTrans, other.settings = list(doClamp = F, pred.type = "cloglog"))
    plot(tmpSeq[, i], suit, type = "l", xlab = colnames(tmpSeq)[i], ylab = "Suitability")
  }
  dev.off()
}

# This function calculates thresholds by controlling the omission rate at a certain percentile.
getThresholds <- function(suitMaps, pres, centile, outDir) {
  thresholdTable <- data.frame(Model = character(),
                               Threshold = numeric())
  for (i in 1:length(suitMaps)) {
    suits <- terra::extract(rast(suitMaps[i]), pres, ID = F)
    p10 <- quantile(suits$layer, centile, na.rm = T, names = F)
    modelName <- rev(strsplit(suitMaps[i], "/")[[1]])[1]
    modelName <- strsplit(modelName, ".", fixed = T)[[1]][1]
    modelName <- strsplit(modelName, "_current", fixed = T)[[1]][1]
    thresholdTable[i, "Model"] <- modelName
    thresholdTable[i, "Threshold"] <- p10
  }
  write.csv(thresholdTable, paste0(outDir, "thresholds.csv"), row.names = F)
}

# This function creates binary maps from suitability maps and calculates projected range size
threshold <- function(suitMaps, thresholds, outDir) {
  for (i in 1:nrow(thresholds)) {
    subStr <- paste0("/", thresholds[i, "Model"], "_")
    maps <- grep(subStr, suitMaps, fixed = T, value = T)
    threshold <- thresholds[i, "Threshold"]
    rangeSizes <- data.frame(Map = character(),
                             RangeSize = numeric(),
                             RangeSizeChange = character())
    for (j in 1:length(maps)) {
      tmpMap <- rast(maps[j])
      binaryMap <- tmpMap > threshold
      binaryMap <- sum(binaryMap)
      outName <- gsub("SuitabilityMaps", "BinaryMaps", maps[j])
      writeRaster(binaryMap, outName, overwrite = T)
      plotName <- gsub(".tif", ".png", outName, fixed = T)
      png(filename = plotName, width = 10, height = 10, units = "in", res = 320)
      plot(binaryMap, range = c(0, 1))
      dev.off()
      rangeSizes[j, "Map"] <- strsplit(rev(strsplit(maps[j], "/")[[1]])[1], ".", fixed = T)[[1]][1]
      rangeSizes[j, "RangeSize"] <- rangeSize <- freq(binaryMap, digits = NA, value = 1, bylayer = F)[1, 2]
    }
    rangeSizes[, "RangeSizeChange"] <- paste0(round((rangeSizes[, "RangeSize"] - rangeSizes[1, "RangeSize"]) / rangeSizes[1, "RangeSize"] * 100, 1), "%")
    write.csv(rangeSizes, paste0(outDir, rangeSizes[1, 1], ".csv"), row.names = F)
  }
}

# This function calculates AICc weights (cf. Wagenmakers & Farrell, 2004) for a set of suitability maps
calculateWeights <- function(suitMaps, modelQuality) {
  AIC <- read.csv(modelQuality, header = T)
  weights <- data.frame(VariableSet = character(length(suitMaps)),
                        GCM = character(length(suitMaps)),
                        AICc = numeric(length(suitMaps)))
  for (i in 1:length(suitMaps)) {
    modelName <- unlist(strsplit(suitMaps[i], "/"))
    modelName <- modelName[length(modelName)]
    modelName <- unlist(strsplit(modelName, "_"))
    weights[i, "VariableSet"] <- toupper(modelName[1])
    if (modelName[1] %in% c("a", "as", "ac", "asc")) {
      weights[i, "GCM"] <- NA
      weights[i, "AICc"] <- subset(AIC, VariableSet == weights[i, "VariableSet"])$AICc
    }
    else {
      weights[i, "GCM"] <- toupper(modelName[2])
      weights[i, "AICc"] <- subset(AIC, VariableSet == weights[i, "VariableSet"] & GCM == weights[i, "GCM"])$AICc
    }
  }
  weights$DeltaAICc <- weights$AICc - min(weights$AICc)
  weights$weightNum <- exp(-0.5*weights$DeltaAICc)
  weights$weightDenom <- sum(weights$weightNum)
  weights$weight <- weights$weightNum / weights$weightDenom
  return(weights)
}
