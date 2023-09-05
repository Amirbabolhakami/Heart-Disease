Missing_Value <- function(Cleveland) {
  colnames(Cleveland)[14] = 'num'
  A_ca = Cleveland[is.na(Cleveland$ca),]
  B_Thal = Cleveland[is.na(Cleveland$thal),]
  Cleveland1 = Cleveland[!(is.na(Cleveland$ca) | is.na(Cleveland$thal)),]
  #Predictive Model for ca
  Model = svm(Cleveland1$ca~., Cleveland1)
  Pre_Model = predict(Model, A_ca[,-12])
  A_ca$ca = round(Pre_Model)
  Cleveland1 = rbind(Cleveland1, A_ca)
  
  #Predictive Model for Thal
  Model = svm(Cleveland1$thal~., Cleveland1)
  Pre_Model = predict(Model, B_Thal[,-13])
  B_Thal$thal = round(Pre_Model)
  Cleveland1 = rbind(Cleveland1, B_Thal)
  
  return(Cleveland1)
}

Feature_Selection <- function(Cleveland1) {
  
  #///////////////Pearson
  FS_Pearson =cor(Cleveland1, method = "pearson")
  FS_Pearson = round(FS_Pearson, digits = 3)
  png(filename = "FS_Pearson.png")
  print(corrplot::corrplot(FS_Pearson))
  dev.off()
  
  #///////////////Information Gain
  FS_Gain <-information.gain(num~., Cleveland1, unit = "log2")
  print(FS_Gain)
  
  #///////////////PCA
  FS_Pca = princomp(Cleveland1[,-14], 
                      cor = TRUE, scores = TRUE)
  summary(FS_Pca)
  loadings(FS_Pca)
  png(filename = "FS_PCA.png")
  screeplot(FS_Pca, type="line", main="Scree Plot", col = 4)
  dev.off()
  
  png(filename = "FS_PCA1.png")
  fviz_pca_var(FS_Pca,
               col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE)
  dev.off()
  
  write.csv(FS_Pearson, "FS_Pearson.csv", row.names = FALSE)
  write.csv(FS_Gain, "FS_Gain.csv", row.names = FALSE)
  write.csv(FS_Pca$loadings, "FS_PCA.csv", row.names = FALSE)
  
  Cleveland1 = Cleveland1[,-c(2, 4:7)]
  return(Cleveland1)
  
}

Discretization <- function(Cleveland1) {
  
  K_Model = kmeans(Cleveland1$age, 4, 20)
  print(K_Model)
  print(K_Model$size)
  Cleveland1$age = K_Model$cluster

  K_Model = kmeans(Cleveland1$thalach, 4, 20)
  print(K_Model)
  print(K_Model$size)
  Cleveland1$thalach = K_Model$cluster
  
  K_Model = kmeans(Cleveland1$oldpeak, 4, 20)
  print(K_Model)
  print(K_Model$size)
  Cleveland1$oldpeak = K_Model$cluster
  
  for (i in 1:ncol(Cleveland1)) {
    Cleveland1[,i] = as.integer(Cleveland1[,i])
  }
  
  return(Cleveland1)
  
}

Outlier <- function(Cleveland1) {
  
  for (i in 1:ncol(Cleveland1)) {
    print(table(Cleveland1[,i]))
    print(colnames(Cleveland1)[i])
    print("****************")
  }
  
  Cleveland1 = Cleveland1[,-4]
  
  #///////// Column5 - 0 | 1
  Cleveland1[Cleveland1$ca != 0, 6] = 1
  
  #///////// Column6 - 0 | 1
  Cleveland1[Cleveland1$thal == 6, 7] = 7
  
  return(Cleveland1)
}

Predictive_Model_DT <- function(Cleveland1) {
  for (i in 1:ncol(Cleveland1)) {
    Cleveland1[,i] = as.factor(Cleveland1[,i])
  }
  
  Train_J48 <- train(Cleveland1[,-8], Cleveland1$num, "J48", 
                     trControl = trainControl(method = "cv", number = 10))
  PM_J48 <- predict(Train_J48, Cleveland1[,-8])
  
  
  return(Cleveland1)
}

Predictve_Model_ANN <- function(Cleveland2) {
  for (i in 1:7) {
    Cleveland2[,i] = as.numeric(Cleveland2[,i])
  }
  Model1 = neuralnet(Cleveland2$num~., Cleveland2, hidden = c(3))
  plot(Model1)
  comp <- compute(Model1, Cleveland2[-8])
  pred.weights <- comp$net.result
  idx <- apply(pred.weights, 1, which.max)
  Cleveland2$Pred = idx
  Cleveland2[Cleveland2$Pred == 1, 9] = 0
  Cleveland2[Cleveland2$Pred == 2, 9] = 1
  Cleveland2$Pred = as.factor(Cleveland2$Pred)
  
  return(Cleveland2)
}


HybridModel <- function(y, models, lambda = NULL, a.args, e.args, n.args = NULL, 
                        s.args = NULL, t.args = NULL, z.args = NULL,
                        weights, errorMethod, rolling = FALSE, cvHorizon,
                        windowSize, horizonAverage = FALSE,
                        parallel = FALSE, num.cores, verbose = TRUE) {
  
  modelArguments <- list("a" = a.args, "e" = e.args, "f" = NULL, "n" = n.args,
                         "s" = s.args, "t" = t.args, "z" = z.args)
  
  
  y <- prepareTimeseries(y = y)
  
  weights <- match.arg(weights)
  if (weights == "insample.errors") {
    wrnMsg <- paste0("Using insample.error weights is not recommended for accuracy and ",
                     "may be deprecated in the future.")
    warning(wrnMsg)
  }
  errorMethod <- match.arg(errorMethod)
  
  # Match the specified models
  expandedModels <- sort(unique(tolower(unlist(strsplit(models, split = "")))))
  # Check models and data length to ensure enough data: remove models that require more data
  expandedModels <- removeModels(y = y, models = expandedModels)
  
  # Check the parallel arguments
  checkParallelArguments(parallel = parallel, num.cores = num.cores)
  
  # Check a.args/t.args/e.args/n.args/s.args
  checkModelArgs(modelArguments = modelArguments, models = expandedModels)
  
  if (weights == "cv.errors" && errorMethod == "MASE") {
    warning("cv errors currently do not support MASE. Reverting to RMSE.")
    errorMethod <- "RMSE"
  }
  
  ##############################################################################
  # Fit models
  ##############################################################################
  
  if (parallel) {
    if (.Platform$OS.type == "unix") {
      cl <- parallel::makeForkCluster(num.cores)
    } else {
      cl <- parallel::makeCluster(num.cores)
    }
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    currentModel <- NULL
    fitModels <- foreach::foreach(modelCode = expandedModels,
                                  .packages = c("forecast", "forecastHybrid")) %dopar% {
                                    # thetam() currently does not handle arguments
                                    if (modelCode == "f") {
                                      fitModel <- thetam(y)
                                    } else { # All other models handle lambda and additional arguments
                                      argsAdditional <- modelArguments[[modelCode]]
                                      if (is.null(argsAdditional)) {
                                        argsAdditional <- list(lambda = lambda)
                                      } else if (is.null(argsAdditional$lambda)) {
                                        argsAdditional$lambda <- lambda
                                      }
                                      currentModel <- purrr::partial(getModel(modelCode), y = y)
                                      fitModel <- do.call(currentModel, argsAdditional)
                                    }
                                    fitModel
                                  }
    modelResults <- unwrapParallelModels(fitModels, expandedModels)
  } else {# serial execution
    modelResults <- list()
    for (modelCode in expandedModels) {
      modelName <- getModelName(modelCode)
      if (verbose) {
        message("Fitting the ", modelName, " model")
      }
      # thetam() currently does not handle arguments
      if (modelCode == "f") {
        modelResults[[modelName]] <- thetam(y)
      } else { # All other models handle lambda and additional arguments
        argsAdditional <- modelArguments[[modelCode]]
        if (is.null(argsAdditional)) {
          argsAdditional <- list(lambda = lambda)
        } else if (is.null(argsAdditional$lambda)) {
          argsAdditional$lambda <- lambda
        }
        currentModel <- purrr::partial(getModel(modelCode), y = y)
        modelResults[[modelName]] <- do.call(currentModel, argsAdditional)
      }
    }
  }
  
  ##############################################################################
  # Determine model weights
  ##############################################################################
  
  includedModels <- names(modelResults)
  numModels <- length(expandedModels)
  if (weights == "equal") {
    modelResults$weights <- rep(1 / numModels, numModels)
  } else if (weights %in% c("insample.errors", "cv.errors")) {
    modelResults$weights <- rep(0, numModels)
    modResults <- modelResults
    
    if (weights == "cv.errors") {
      for (modelCode in expandedModels) {
        modelName <- getModelName(modelCode)
        currentModel <- getModel(modelCode)
        if (verbose) {
          message("Cross validating the ", modelName, " model")
        }
        modResults[[modelName]] <- cvts(y, FUN = currentModel,
                                        maxHorizon = cvHorizon,
                                        horizonAverage = horizonAverage,
                                        rolling = rolling,
                                        verbose = FALSE,
                                        windowSize = windowSize,
                                        num.cores = num.cores)
      }
    }
    cvHorizon <- ifelse(horizonAverage, 1, cvHorizon)
    cvHorizon <- ifelse(weights != "cv.errors", 1, cvHorizon)
    
    weightFunction <- function(x) accuracy(modResults[[getModelName(x)]])[cvHorizon, errorMethod]
    modelResults$weights <- sapply(expandedModels, weightFunction)
    
    inverseErrors <- 1 / modelResults$weights
    modelResults$weights <- inverseErrors / sum(inverseErrors)
    
  }
  
  if (is.element(NaN, modelResults$weights) & weights %in% c("insample.errors", "cv.errors")) {
    wrnMsg <- paste0("At least one model perfectly fit the series, so accuracy measures cannot",
                     " be used for weights. Reverting to weights = \"equal\".") # nolint
    warning(wrnMsg)
    modelResults$weights <- rep(1 / length(includedModels),
                                length(includedModels))
  }
  names(modelResults$weights) <- includedModels
  
  ##############################################################################
  # Prepare hybridModel object
  ##############################################################################
  
  fits <- sapply(includedModels, FUN = function(x) fitted(modelResults[[x]]))
  fitsWeightsMatrix <- matrix(rep(modelResults$weights[includedModels],
                                  times = nrow(fits)),
                              nrow = nrow(fits), byrow = TRUE)
  fits <- rowSums(fits * fitsWeightsMatrix)
  resid <- y - fits
  tsp(fits) <- tsp(y)
  
  xregs <- list()
  if ("a" %in% expandedModels) {
    xregs$auto.arima <- ifelse("xreg" %in% names(a.args) && !is.null(a.args$xreg), TRUE, FALSE) # nolint
  }
  if ("n" %in% expandedModels) {
    xregs$nnetar <- ifelse("xreg" %in% names(n.args) && !is.null(n.args$xreg), TRUE, FALSE)
  }
  if ("s" %in% expandedModels) {
    methodArima <- "method" %in% names(s.args) && s.args$method == "arima"
    xregs$stlm <- ifelse("xreg" %in% names(s.args) && !is.null(s.args$xreg) && methodArima,
                         TRUE, FALSE)
  }
  
  class(modelResults) <- "hybridModel"
  modelResults$frequency <- frequency(y)
  modelResults$x <- y
  modelResults$xreg <- xregs
  modelResults$models <- includedModels
  modelResults$fitted <- fits
  modelResults$residuals <- resid
  return(modelResults)
}




