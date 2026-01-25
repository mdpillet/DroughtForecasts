library(terra)
library(dplyr)
library(ggplot2)
library(glmnet)
library(usdm)
library(cluster)
library(dismo)
library(doParallel)
library(pROC)
library(ecospat)

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
  preds <- mostPreds
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
  return(1)
}

# This function crops a predictor file based on a buffer
cropPredictors <- function(preds, buffer, outDir, bufferName) {
  outName <- strsplit(preds, "/", fixed = T)
  outName <- outName[[1]][length(outName[[1]])]
  tmp <- gsub(".shp", "", gsub("_", "", bufferName, fixed = T), fixed = T)
  tmp <- paste0("_", tmp)
  outName <- gsub(".tif", paste0(tmp, ".tif"), outName, fixed = T)
  preds <- rast(preds)
  croppedPreds <- crop(preds, buffer, mask = T)
  writeRaster(croppedPreds, paste0(outDir, outName), overwrite = T)
  return(TRUE)
}

# This function samples background points from an environmental layers file
sampleBackground <- function(preds, no) {
  preds <- raster::stack(preds)
  bgSample <- dismo::randomPoints(preds, no, lonlatCorrection = F, warn = 0)
  return(bgSample)
}

# This function removes correlated predictors
removeCorrelatedPredictors <- function(env, corrThreshold, outName, seedNo) {
  set.seed(seedNo)
  excluded <- withCallingHandlers(
    usdm::vifcor(env, th = corrThreshold)@excluded,
    warning = function(w) {
      if (grepl("essentially perfect fit", w$message) ||
          grepl("standard deviation is zero", w$message) ||
          grepl("fewer values returned", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!is.null(excluded)) env <- env[[-which(names(env) %in% excluded)]]
  writeRaster(env, outName, overwrite = T)
  return(excluded)
}

# This function spatially stratifies occurrence data
spatialStratify <- function(sp.pts, seedNo, minClusters, maxClusters, bootstraps, finalFolds) {
  # Decide on initial cluster size based on gap statistics
  set.seed(seedNo)
  gap <- clusGap(x = sp.pts, K.max = min(nrow(sp.pts) - 1, maxClusters), B = bootstraps, verbose = F,
                 FUNcluster = function(x, k) {
                   withCallingHandlers(kmeans(x, centers = k, iter.max = 1000, nstart = 10),
                                       warning = function(w) {
                                         if (grepl("did not converge", w$message) | grepl("exceeded maximum", w$message)) invokeRestart("muffleWarning")
                                       })})
  optimal_k <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = "firstSEmax")
  chosen_k <- min(max(minClusters, optimal_k), nrow(sp.pts) - 1)
  folds <- withCallingHandlers(
    kmeans(sp.pts, centers = chosen_k, nstart = 10, iter.max = 1000),
    warning = function(w) {
      if (grepl("did not converge", w$message) | grepl("exceeded maximum", w$message)) invokeRestart("muffleWarning")
    }
  )
  # Combine initial clusters to get finalFolds clusters while balancing cluster sizes
  focal.pres <- sp.pts
  combine.subclusters <- function(n.sub, folds, folds2) {
    combine.folds <- lapply(1:bootstraps, function(z) {
      vals <- sample(1:nrow(folds.tmp$centers), n.sub, replace = FALSE)
      assigns <- matrix(NA, nrow = finalFolds, ncol = ceiling(n.sub / finalFolds))
      fullCols <- floor(n.sub / finalFolds) * finalFolds
      assigns[1:fullCols] <- vals[1:fullCols]
      if (sum(is.na(assigns)) > 0) {
        remaining <- vals[(fullCols + 1):n.sub]
        assigns[sample(1:finalFolds, length(remaining), replace = FALSE), ncol(assigns)] <- remaining  
      }
      return(assigns)
    })
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
  folds.tmp <- folds
  focal.pres$folds <- combine.subclusters(chosen_k, folds.tmp$clust, folds.tmp$clust)
  return(list(focal.pres = focal.pres,
              gap = gap,
              optimal_k = optimal_k,
              chosen_k = chosen_k))
}

# Helper functions for construction of model matrix
categoricalval <- function(x, category) {
  ifelse(x == category, 1, 0)
}

thresholdval <- function(x, knot) {
  ifelse(x >= knot, 1, 0)
}

hingeval <- function(x, min, max) {
  pmin(1, pmax(0, (x - min) / (max - min)))
}

# Helper function for model prediction
predict.cv.maxnet <- function(object,
                              newdata,
                              clampVars = F, 
                              clampFeatures = F, 
                              clampRule = 1,
                              type = c("link", "exponential", "cloglog", "logistic"),
                              offset = NULL,
                              ...) {
  # Ensure NA values are returned
  na_action <- options("na.action")[[1]]
  on.exit(options(na.action = na_action))
  options(na.action = "na.pass")
  # Format as data frame
  newdata <- data.frame(newdata)
  # Clamp predictor variables
  if (clampVars) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v] - clampRule * object$varsd[v]), object$varmax[v] + clampRule * object$varsd[v])
    }
  }
  # Construct model matrix
  b.names <- names(object$betas)
  to.combine <- b.names[which(duplicated(b.names))]
  if (length(to.combine) > 0) {
    for (tc in to.combine) {
      tmp <- which(b.names %in% tc)
      object$betas[tmp] <- sum(object$betas[tmp])
      object$betas <- object$betas[-tmp[-1]]
      b.names <- names(object$betas)
    }
  }
  if (!length(unique(names(object$betas))) == length((names(object$betas)))) print("Duplicated feature names exist.")
  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)", names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\\2)", terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)", terms)
  form <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  mm <- model.matrix(form, newdata)
  # Clamp features
  if (clampFeatures) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)] - clampRule * object$featuresd[names(object$betas)]),
                 object$featuremaxs[names(object$betas)] + clampRule * object$featuresd[names(object$betas)]))
  } 
  # Return predictions
  link <- (mm %*% object$betas) + object$alpha + ifelse(!is.null(offset), offset, 0)
  type <- match.arg(type)
  if (type == "link") predOut <- link
  if (type == "exponential") predOut <- exp(link)
  if (type == "cloglog") predOut <- 1 - exp(0 - exp(object$entropy + link))
  if (type == "logistic") predOut <- 1 / (1 + exp(-object$entropy - link))
  return(predOut)
}

# Helper function for forcing quadratic coefficients to be negative
quadraticUpperLimit0 <- function(m) {
  up <- rep(Inf, ncol(m))
  isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
  to0 <- isquadratic(colnames(m))
  up[to0] <- 0
  return(up)
}

# Helper function for model fitting
cv.maxnet <- function(p.f,
                      data.f,
                      form.f = maxnet::maxnet.formula(p.f, data),
                      regmult = 1,
                      regfun = maxnet::maxnet.default.regularization,
                      myweights=p.f + (1 - p.f) * 100,
                      offset.f = NULL,
                      maxit = 5e5, 
                      standardize = TRUE,
                      foldid = NULL,
                      parallel = 0,
                      quadraticUpperLimit0 = FALSE,
                      lambdaSeq = NULL,
                      ...) {
  # Check for NAs
  keep <- which(complete.cases(data.f))
  data.f <- data.f[keep,]
  p.f <- p.f[keep]
  myweights <- myweights[keep]
  offset.f <- offset.f[keep]
  # Set model matrix and regularization vector
  mm <- model.matrix(form.f, data.f)
  reg <- regfun(p = p.f, mm) * regmult
  # Set glmnet options
  glmnet::glmnet.control(pmin = 1e-08, fdev = 1e-5)
  # Set upper limits on quadratic coefficients to force them to be negative
  if (quadraticUpperLimit0) {
    upper.limits <- quadraticUpperLimit0(mm)
  } else {
    upper.limits <- rep(Inf, ncol(mm))
  }
  # Set lambda sequence
  myLambda <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) * sum(p.f) / sum(myweights)
  # Build cross-validated model
  if (parallel > 0) {
    cl <- makeCluster(parallel)
    registerDoParallel(cl)
  }
  model.cv <- cv.glmnet(x = mm, 
                        y = p.f,
                        family = "binomial",
                        standardize = standardize,
                        penalty.factor = reg,
                        lambda = myLambda,
                        weights = myweights,
                        offset = offset.f,
                        maxit = maxit,
                        foldid = foldid,
                        parallel = ifelse(parallel == 0, F, T),
                        trace.it = 0,
                        upper.limits = upper.limits,
                        type.measure = "deviance",
                        ...) 
  if (parallel) stopCluster(cl)
  
  # Lambda selection
  lam.min.index <- which(model.cv$lambda == model.cv$lambda.min)
  # Use largest lambda with nonzero coefficients if lambda min results in zero coefficients
  if (model.cv$nzero[[lam.min.index]] == 0) {
    model.cv$lambdaRule <- "largestnonzero"
    lambda <- ifelse(any(model.cv$nzero > 0), max(model.cv$lambda[model.cv$nzero > 0]), NA)
    model.cv$lambda.value <- lambda
  } else {
    lam.1se.index <- which(model.cv$lambda == model.cv$lambda.1se)
    if (!model.cv$nzero[[lam.1se.index]] == 0) {
      lambda <- model.cv$lambda.1se
      model.cv$lambda.value <- lambda
      model.cv$lambdaRule <- "1se"
    } else {
      lam.min.index = which(model.cv$lambda == model.cv$lambda.min)
      se <- model.cv$cvsd[lam.min.index]
      lam.cut <- model.cv$lambda.min + 0.5 * se
      which.cut.rev <- findInterval(lam.cut, rev(model.cv$lambda))
      lambda.val.tmp <- rev(model.cv$lambda)[which.cut.rev]
      which.cut <- which(model.cv$lambda ==	lambda.val.tmp)
      if (!model.cv$nzero[[which.cut]] == 0) {
        lambda <- lambda.val.tmp
        model.cv$lambda.value <- lambda
        model.cv$lambdaRule <- "0.5se"
      } else {
        lam.cut <- model.cv$lambda.min + 0.25 * se
        which.cut.rev <- findInterval(lam.cut,rev(model.cv$lambda))
        lambda.val.tmp <- rev(model.cv$lambda)[which.cut.rev]
        which.cut <- which(model.cv$lambda ==	lambda.val.tmp)
        if (!model.cv$nzero[[which.cut]] == 0) {
          lambda <- lambda.val.tmp
          model.cv$lambda.value <- lambda
          model.cv$lambdaRule <- "0.25se"
        } else {
          lambda <- model.cv$lambda.min
          model.cv$lambda.value <- lambda
          model.cv$lambdaRule <- "min"
        }
      }
    }
  }
  model.cv$lambdaKeepIndex <- ifelse(!is.na(model.cv$lambda.value), which(model.cv$lambda == lambda), NA)
  # Store coefficients for selected lambda
  if (!is.na(model.cv$lambdaKeepIndex)) {
    bb <- model.cv$glmnet$beta[, model.cv$lambdaKeepIndex]
    model.cv$betas <- bb[bb != 0] # This may lead to bugs when using hinges - not applicable here.
  } else {
    model.cv$betas <- numeric(0)
  }
  model.cv$alpha <- 0
  # Extract information on predictors
  model.cv$nvarraw <- ncol(data.f)
  model.cv$ncoefraw <- ncol(mm)
  model.cv$varraw_drought <- ifelse(any(grepl("HAS|HAD|HMS|HMD", colnames(data.f))), T, F)
  finalVars <- setNames(as.integer(sapply(colnames(data.f), function(x) any(grepl(x, names(model.cv$betas))))), colnames(data.f))
  model.cv$nvarfit <- sum(finalVars)
  model.cv$varnamesfit <- names(finalVars[finalVars == 1])
  model.cv$ncoeffit <- length(model.cv$betas)
  model.cv$varfit_drought <- ifelse(any(grepl("HAS|HAD|HMS|HMD", names(model.cv$betas))), T, F)
  model.cv$varfit_drought <- ifelse(model.cv$nvarfit == 0, NA, model.cv$varfit_drought)
  # Calculate EPV
  model.cv$EPV <- ifelse(model.cv$nvarfit > 0, sum(p.f) / model.cv$nvarfit, NA)
  # Get entropy and alpha
  if (!is.na(model.cv$lambda.value)) {
    # rr <- predict.cv.maxnet(object = model.cv, newdata = data.frame(data.f[p.f == 0, , drop = FALSE]), type = "exponential", clampVars = F, clampFeatures = F, offset = offset.f[p.f == 0, drop = FALSE])
    # raw <- rr / sum(rr)
    # model.cv$entropy <- -sum(raw * log(raw))
    # model.cv$alpha <- -log(sum(rr))
    # The code above results in overflow when link (see predict.cv.maxnet function code) contains numbers above ~709. Workaround is to directly request the link and rescale it.
    # Additionally, if any values in raw are zero, then log(0) results in -Inf. Since the limit of this is zero, these cases are removed for the purposes of entropy calculation.
    rr <- predict.cv.maxnet(object = model.cv, newdata = data.frame(data.f[p.f == 0, , drop = FALSE]), type = "link", clampVars = F, clampFeatures = F, offset = offset.f[p.f == 0, drop = FALSE])
    rr_max <- max(rr)
    exp_rr <- exp(rr - rr_max)
    raw <- exp_rr / sum(exp_rr)
    model.cv$entropy <- -sum(raw[raw > 0] * log(raw[raw > 0]))
    model.cv$alpha <- -(rr_max + log(sum(exp_rr)))
  } else {
    model.cv$entropy <- NA
    model.cv$alpha <- NA
  }
  # Get penalty factor
  model.cv$penalty.factor <- reg
  # Get variable and feature extremes for clamping
  model.cv$featuremins <- apply(mm, 2, min)
  model.cv$featuremaxs <- apply(mm, 2, max)
  model.cv$featuresd <- apply(mm, 2, sd)
  vv <- (sapply(data.f, class) != "factor")
  model.cv$varmin <- apply(data.f[, vv, drop = FALSE], 2, min)
  model.cv$varmax <- apply(data.f[, vv, drop = FALSE], 2, max)
  model.cv$varminobs <- apply(data.f[p.f == 1, vv, drop = FALSE], 2, min)
  model.cv$varmaxobs <- apply(data.f[p.f == 1, vv, drop = FALSE], 2, max)
  model.cv$varsd <- apply(data.f[, vv, drop = FALSE], 2, sd)
  # Get variable and feature means
  means <- apply(data.f[p.f == 1, vv, drop = FALSE], 2, mean)
  model.cv$samplemedians <- apply(data.f[p.f == 1, vv, drop = FALSE], 2, median)
  majorities <- sapply(names(data.f)[!vv], function(n) which.max(table(data.f[p.f == 1, n, drop = FALSE])))
  names(majorities) <- names(data.f)[!vv]
  model.cv$samplemeans <- c(means, majorities)
  model.cv$levels <- lapply(data.f, levels)
  # Return model
  return(model.cv)
}

# This function calculates BIC of a cv.glmnet model, following Renner et al. (2021).
cv_bic <- function(fit, finalPres) {
  whlm <- which(fit$lambda == fit$lambda.value)
  with(fit$glmnet.fit,
       {
         tLL <- nulldev - nulldev * (1 - dev.ratio)[whlm]
         k <- df[whlm]
         return(log(finalPres) * k - tLL)
       })
}

# This function fits Maxent models
fitModel <- function(pres, preds, bg, sp, tune.args, numCores, folds, varImportance = 0, modelStats = T, statsDistances = NULL, calcThresholds = T) {
  foldCount <- length(unique(folds))
  response <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
  presVals <- extract(preds, pres, ID = F)
  bgVals <- extract(preds, bg, ID = F)
  predVals <- rbind(presVals, bgVals)
  bgFolds <- sample(rep(1:foldCount, length.out = nrow(bg)))
  modelFit <- try(cv.maxnet(p.f = response,
                            data.f = predVals,
                            form.f = maxnet::maxnet.formula(response, predVals, classes = tune.args),
                            regmult = 1,
                            regfun = maxnet::maxnet.default.regularization,
                            myweights = response + (1 - response) * 100,
                            standardize = T,
                            foldid = c(folds, bgFolds),
                            parallel = ifelse(numCores > 1, numCores, 0),
                            quadraticUpperLimit0 = T))
  # Calculate variable importance
  modelFit$varImportance <- NA
  if (!inherits(modelFit, "try-error") & varImportance > 0 & length(modelFit$varnamesfit) > 0) {
    varList <- modelFit$varnamesfit
    importance <- vector("numeric", length(varList))
    names(importance) <- varList
    envData <- predVals
    fullModelVals <- predict.cv.maxnet(object = modelFit, newdata = envData, type = "cloglog", clampVars = F, clampFeatures = F)
    for (thisVar in varList) {
      correls <- vector("numeric", varImportance)
      origVals <- envData[, thisVar]
      for (thisRep in 1:varImportance) {
        permInd <- sample(1:nrow(envData), nrow(envData))
        envData[, thisVar] <- origVals[permInd]
        newModelVals <- predict.cv.maxnet(object = modelFit, newdata = envData, type = "cloglog", clampVars = F, clampFeatures = F)
        correls[thisRep] <- cor(fullModelVals, newModelVals)
      }
      envData[, thisVar] <- origVals
      importance[thisVar] <- mean(correls)
    }
    importance <- 1 - importance
    modelFit$varImportance <- importance / sum(importance)
  }
  # Store model predictions for model statistic or threshold calculation
  if (!inherits(modelFit, "try-error") & (modelStats | calcThresholds) & length(modelFit$varnamesfit) > 0) {
    allPredictions <- as.numeric(predict.cv.maxnet(object = modelFit, newdata = predVals, type = "cloglog", clampVars = F, clampFeatures = F))
    presPredictions <- allPredictions[response == 1]
    absPredictions <- allPredictions[response == 0]
  }
  # Calculate model statistics
  modelFit$AUC <- NA
  modelFit$boyce <- NA
  modelFit$cv_bic <- NA
  if (!inherits(modelFit, "try-error") & modelStats & length(modelFit$varnamesfit) > 0) {
    modelFit$AUC <- as.numeric(auc(roc(response, allPredictions, quiet = T)))
    modelFit$boyce <- ecospat.boyce(fit = allPredictions, obs = presPredictions, PEplot = F)$cor
    modelFit$cv_bic <- cv_bic(modelFit, nrow(presVals))
    # Calculate distance-based AUC and Boyce index
    if (!is.null(statsDistances)) {
      distance_matrix <- distance(bg, pres)
      min_distances <- apply(distance_matrix, 1, min)
      for (distance in statsDistances) {
        background_filtered <- min_distances >= distance
        allPredictions_filtered <- c(presPredictions, absPredictions[background_filtered])
        response_filtered <- c(response[response == 1], rep(0, length(absPredictions[background_filtered])))
        modelFit[[paste0("bgNo_", distance)]] <- sum(background_filtered)
        modelFit[[paste0("AUC_", distance)]] <- as.numeric(auc(roc(response_filtered, allPredictions_filtered, quiet = T)))
        modelFit[[paste0("boyce_", distance)]] <- ecospat.boyce(fit = allPredictions_filtered, obs = presPredictions, PEplot = F)
      }  
    }
  }
  # Calculate thresholds
  modelFit$maxTSS <- NA
  modelFit$eqsensspec <- NA
  modelFit$OR10 <- NA
  if (!inherits(modelFit, "try-error") & calcThresholds & length(modelFit$varnamesfit) > 0) {
    dismoEval <- evaluate(presPredictions, absPredictions)
    modelFit$maxTSS <- dismo::threshold(dismoEval, stat = "spec_sens")
    modelFit$eqsensspec <- dismo::threshold(dismoEval, stat = "equal_sens_spec")
    modelFit$OR10 <- dismo::threshold(dismoEval, stat = "sensitivity", sensitivity = 0.9)
  }
  return(modelFit)
}