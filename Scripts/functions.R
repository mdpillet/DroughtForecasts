library(terra)
library(dplyr)
library(paran)
library(ggplot2)
library(flexclust)
library(glmnet)

# This function predicts with a cross-validated Maxent model
# (Code for fitting cross-validated Mexent models is implemented in the custom fork of the ENmeval package)
predict.CVmaxnet <- function(
    object, newdata, clamp = TRUE,
    type = c("link", "exponential", "cloglog", "logistic"),
    offset = NULL, ...) {
  # |GEPB|: First lines in predict.maxnet()
  # na_action <- options("na.action")[[1]]
  # on.exit(options(na.action = na_action))
  # options(na.action = "na.pass")
  
  newdata <- data.frame(newdata) # sparse matrices cause issues
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }
  
  # maxnet has a bug where some hinges are duplicated, so combining them
  #  cleaning up the formatting of these types for use with formula
  b.names <- names(object$betas)
  to.combine <- b.names[which(duplicated(b.names))]
  # CHECK GEPB betas
  if (length(to.combine) > 0) {
    for (tc in to.combine) {
      tmp <- which(b.names %in% tc)
      object$betas[tmp] <- sum(object$betas[tmp])
      object$betas <- object$betas[-tmp[-1]]
      b.names <- names(object$betas)
    }
  }
  if (!length(unique(names(object$betas))) == length((names(object$betas)))) {
    print('some duplicated feature names still')
  }
  
  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\\2)",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  # changed f to form because f is a command when using browser()
  form <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  
  mm <- model.matrix(form, data.frame(newdata))
  if (clamp) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  }
  # CM added offset; CM already logs this on input, so this isn't general
  # CM ditched alpha because it required a wasted projection step in maxnet, and
  # its just a normalizing constant that you can take care of later, i think.
  # maybe its an issue for the rescaled ones, but don't care...
  t1 <- proc.time()
  link <- (mm %*% object$betas) + ifelse(!is.null(offset), offset, 0) + object$alpha
  t2 <-  proc.time() - t1
  type <- match.arg(type)
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    out <- exp(link)
    return(out)
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1 / (1 + exp(-object$entropy - link)))
  }
}

categoricalval <- function(x, category) {
  ifelse(x == category, 1, 0)
}

thresholdval <- function(x, knot) {
  ifelse(x >= knot, 1, 0)
}

hingeval <- function(x, min, max) {
  pmin(1, pmax(0, (x - min) / (max - min)))
}

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
  if (distance > 0) { 
    buffer <- buffer(chull, distance)
  } else {
    buffer <- chull
  }
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
sampleBackground <- function(preds, no) {
  preds <- raster::stack(preds)
  bgSample <- dismo::randomPoints(preds, no, lonlatCorrection = F)
  return(bgSample)
}

# This function removes correlated predictors
removeCorrelatedPredictors <- function(env, corrThreshold, outName, seedNo) {
  set.seed(seedNo)
  c1 <- cor(values(env), use = "complete.obs")
  # Toss variables that were all NULL, NA, NaN, or 0 SD
  tossed <- NULL
  bad <- which(apply(c1, 1, function(x) {
    all(is.na(x)) | all(is.nan(x)) | all(is.null(x)) | sum(x, na.rm = T) == 1
  }))
  if (length(bad) > 0) {
    tossed <- c(tossed, names(bad))
    c1.index <- which(colnames(c1) %in% tossed)
    c1 <- c1[-c1.index, -c1.index]
  }
  # Toss variables that are too correlated with other variables
  too.cor <- apply(abs(c1) > corrThreshold, 1, sum) - 1
  
  while (any(too.cor >= 1)) {
    most.cor <- too.cor[too.cor == max(too.cor)]
    toss <- sample(1:length(most.cor), 1)
    tossed <- c(tossed, names(most.cor[toss]))
    c1.index <- which(colnames(c1) %in% tossed)
    c1 <- c1[-c1.index, -c1.index]
    too.cor <- apply(abs(c1) > corrThreshold, 1, sum) - 1
  }
  
  if (!is.null(tossed)) env <- env[[-which(names(env) %in% tossed)]]
  writeRaster(env, outName)
  return(tossed)
}

# This function spatially stratifies occurrence data
spatialStratify <- function(sp.pts, seedNo) {
  set.seed(seedNo)
  out <- try({
    n <- nrow(sp.pts)
    focal.pres <- sp.pts
    combine.subclusters <- function(n.sub, folds, folds2) {
      # Combine subclusters to get desired number of folds
      combine.folds <- lapply(1:50, function(z) { matrix(sample(1:nrow(folds.tmp$centers), n.sub, replace = FALSE), ncol = n.sub / 5)})
      fold.candidates <- lapply(1:length(combine.folds), function(jj) {
        for (ii in 1:nrow(combine.folds[[jj]])) { folds[folds2 %in% c(combine.folds[[jj]][ii,])] = ii }
        folds
      })
      fold.entropy <- sapply(1:length(fold.candidates), function(ii) {
        t <- table(fold.candidates[[ii]]) / length(fold.candidates[[ii]])
        -1*sum(t*log(t))
      })
      best <- which.max(fold.entropy)
      fold.candidates[[best]]
    }
    
    if (n > 5 & n <= 15) {
      folds.tmp <- kmeans(sp.pts, 5)
      folds.tmp1 <- flexclust::as.kcca(folds.tmp, data = sp.pts)
      focal.pres$folds <- folds.tmp$cluster
    }
    if (n > 15 & n <= 30) {
      n.sub <- 10
      folds.tmp <- kmeans(sp.pts, n.sub)
      folds.tmp1 <- flexclust::as.kcca(folds.tmp, data = sp.pts)
      focal.pres$folds <- combine.subclusters(n.sub, folds.tmp$clust, folds.tmp$clust)
    }
    if (n > 30 & n <= 45) {
      n.sub <- 15
      folds.tmp <- kmeans(sp.pts, n.sub)
      folds.tmp1 <- flexclust::as.kcca(folds.tmp, data = sp.pts)
      focal.pres$folds <- combine.subclusters(n.sub, folds.tmp$clust, folds.tmp$clust)
    }
    if (n > 45 & n <= 60) {
      n.sub <- 20
      folds.tmp <- kmeans(sp.pts, n.sub)
      folds.tmp1 <- flexclust::as.kcca(folds.tmp, data = sp.pts)
      focal.pres$folds <- combine.subclusters(n.sub, folds.tmp$clust, folds.tmp$clust)
    }
    if (n > 60) {
      n.sub <- 25
      folds.tmp <- kmeans(sp.pts, n.sub)
      folds.tmp1 <- flexclust::as.kcca(folds.tmp, data = sp.pts)
      focal.pres$folds <- combine.subclusters(n.sub, folds.tmp$clust, folds.tmp$clust)
    }
    
    sp.pts <- focal.pres
    return(sp.pts)
  })
  return(out)
}

# This function fits Maxent models
fitModel <- function(pres, preds, bg, sp, outDir, tune.args, numCores, plot, folds) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  partitions <- list(occs.grp = folds,
                     bg.grp = rep(0, nrow(bg)))
  colnames(bg) <- colnames(pres)
  model <- try(ENMevaluate(occs = pres,
                       envs = preds,
                       bg = bg,
                       tune.args = tune.args,
                       algorithm = "CVmaxnet",
                       partitions = "user",
                       user.grp = partitions,
                       doClamp = F,
                       taxon.name = sp,
                       overlap = F,
                       parallel = F,
                       numCores = numCores,
                       quiet = T))
  if (class(model) == "try-error") {
      model <- NULL
    # other.args <- list(regfun = maxnet.default.regularization)
    # model <- ENMevaluate(occs = pres,
    #                      envs = preds,
    #                      bg = bg,
    #                      tune.args = tune.args,
    #                      algorithm = "maxnet",
    #                      partitions = "user",
    #                      user.grp = partitions,
    #                      doClamp = F,
    #                      taxon.name = sp,
    #                      overlap = F,
    #                      parallel = F,
    #                      numCores = numCores,
    #                      quiet = T,
    #                      other.settings = list(other.args = other.args))
  }
  if (!is.null(model)) save(model, file = paste0(outDir, outName, ".rda"))
  if (plot) {
    modelStats <- evalplot.stats(model, stats = c("or.10p", "auc.val"), x.var = "rm", color = "fc", dodge = 0.5)
    ggsave(paste0(outDir, outName, ".png"), modelStats)
  }
  return(model@algorithm)
}

# This function calculates BIC of a cv.glmnet model, following Renner et al. (2021).
cv_bic <- function(fit, finalPres){
  whlm <- which(fit$lambda == fit$lambda.value)
  with(fit$glmnet.fit,
       {
         tLL <- nulldev - nulldev * (1 - dev.ratio)[whlm]
         k <- df[whlm]
         return(log(finalPres) * k - tLL)
       })
}

# This function calculates BIC of a glmnet model, following Renner et al. (2021).
bic <- function(fit, finalPres){
  with(fit,
       {
         tLL <- nulldev - nulldev * (1 - dev.ratio)[200]
         k <- df[200]
         return(log(finalPres) * k - tLL)
       })
}

# This function compares model sets based on BIC.
compareModelSets <- function(modelSets, nullModels, outDir) {
  nrow <- length(modelSets)
  comp <- data.frame(VariableSet = character(),
                     GCM = character(),
                     Algorithm = character(),
                     LambdaRule = character(),
                     BIC = numeric(),
                     DeltaBIC = numeric(),
                     AverageValidationAUC = numeric(),
                     SignificanceAUC = numeric(),
                     AverageOmissionRate = numeric(),
                     SignificanceOmissionRate = numeric(),
                     Variables = character())
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
    comp[i, "AverageValidationAUC"] <- eval.results(model)$auc.val.avg
    comp[i, "AverageOmissionRate"] <- eval.results(model)$or.10p.avg
    comp[i, "Algorithm"] <- model@algorithm
    if (model@algorithm == "CVmaxnet") {
      comp[i, "LambdaRule"] <- model@models$fc.LQP_rm.1$lambdaRule
      comp[i, "BIC"] <- cv_bic(model@models$fc.LQP_rm.1, nrow(model@occs))
    }
    else {
      comp[i, "LambdaRule"] <- NA
      comp[i, "BIC"] <- bic(model@models$fc.LQP_rm.1, nrow(model@occs))
    }
    comp[i, "Variables"] <- paste0(names(model@models$fc.LQP_rm.1$betas), collapse = ",")  
    if (length(nullModels) == nrow) {
      load(nullModelFiles[i])
      comp[i, "SignificanceAUC"] <- summ.table[summ.table$statistic == "pvalue", "auc.val"]
      comp[i, "SignificanceOmissionRate"] <- summ.table[summ.table$statistic == "pvalue", "or.10p"]
    }
  }
  comp$DeltaBIC <- comp$BIC - min(comp$BIC)
  comp <- arrange(comp, DeltaBIC)
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
  options(na.action = "na.pass")
  suitVals <- predict.CVmaxnet(model@models$fc.LQP_rm.1, newdata = values(preds), clamp = F, type = "cloglog")
  options(na.action = "na.omit")
  suit <- preds[[1]]
  values(suit) <- suitVals
  writeRaster(suit, paste0(outDir, outName, ".tif"), overwrite = T)
  suitMap <- rast(paste0(outDir, outName, ".tif"))
  png(filename = paste0(outDir, outName, ".png"), width = 10, height = 10, units = "in", res = 320)
  plot(suitMap, range = c(0, 1))
  dev.off()
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
  summStats[2,] <- apply(valRast, 2, max, na.rm = T)
  presExtr <- terra::extract(tmpRast, occ, ID = F)
  summStats[3,] <- apply(presExtr, 2, mean, na.rm = T)
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
    summStats[4,] <- pmin(apply(values(tmpFut), 2, min, na.rm = T), summStats[4,], na.rm = T)
    summStats[5,] <- pmin(apply(values(tmpFut), 2, max, na.rm = T), summStats[5,], na.rm = T)
  }
  rownames(summStats) <- c("MinPresent", "MaxPresent", "MeanPresent", "MinFuture", "MaxFuture")
  summStats <- as.data.frame(summStats)
  write.csv(summStats, paste0(outDir, outName, ".csv"))
}

# This function creates marginal effect plots
plotMarginal <- function(summFile, modelFile, outDir) {
  outName <- strsplit(summFile, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  outName <- strsplit(outName, ".", fixed = T)
  outName <- outName[[1]][1]
  summStats <- read.csv(summFile, header = T, row.names = 1)
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
    suit <- predict(model@models$fc.LQP_rm.1, newdata = tmpSeq, clamp = F, type = "cloglog")
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
    p10 <- quantile(suits[[1]], centile, na.rm = T, names = F)
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