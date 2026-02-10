cv.predict.xgboost = function(df,
                              n.folds = NULL,
                              fold.ids,
                              response.col,
                              predictor.cols = NULL,
                              covariate.cols = NULL,
                              data.selection,
                              feature.engineering.col,
                              feature.engineering.row,
                              feature.selection = "none",
                              feature.selection.metric = "sRMSE",
                              feature.selection.metric.threshold = 1,
                              feature.selection.model = "lm",
                              feature.selection.criterion = "relative.gain",
                              feature.selection.include.covariates = TRUE,
                              seed = 12345,
                              n.folds.inner = 5,
                              nrounds = 100,
                              eta.values = c(0.01, 0.1),
                              max.depth.values = c(3, 6),
                              subsample.values = 1,
                              colsample.values = 1,
                              baseline = FALSE,
                              n.cores = 1) {
  # Reproducibility
  set.seed(seed)
  
  # Basic dimensions / predictors
  n <- nrow(df)
  if (is.null(predictor.cols)) {
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  p <- length(pred.names)
  
  pred.names.feature.selection = pred.names[-which(pred.names %in% covariate.cols)]
  
  # Outer folds
  if (!is.null(fold.ids)) {
    n.folds = length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # Containers for predictions and per-fold variable importances
  pred.vec <- rep(NA_real_, n)
  var.imp <- matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) <- pred.names
  
  # --- parallel backend setup ---
  parallel_backend <- FALSE
  n.cores <- as.integer(n.cores)
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    parallel_backend <- TRUE
  }
  # -------------------------------
  
  # Container for variable selection
  feature_table <- data.frame(pred = pred.names.feature.selection, stringsAsFactors = FALSE)
  
  # Outer CV loop
  for (fold in seq_len(n.folds)) {
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    
    # prepare inner fold ids (reproducible per outer fold)
    set.seed(seed + fold)
    n.inner <- nrow(df.train)
    inner.fold.ids <- sample(rep(seq_len(n.folds.inner), length.out = n.inner))
    
    ## Feature selection (performed inside each outer fold)
    if (feature.selection == "none") {
      # No selection: keep all columns except participant_id
      pred.selected = df.train %>%
        select(-participant_id, -all_of(response.col)) %>%
        colnames()
      
      feature.selection.res <- NULL
      
    } else if (feature.selection == "univariate") {
      # Univariate feature selection using inner CV
      feature.selection.res = feature.selection.univariate(
        df = df.train,
        response.col = response.col,
        covariate.cols = covariate.cols,
        predictor.cols = pred.names.feature.selection,
        model = feature.selection.model,
        metric = feature.selection.metric,
        metric.threshold = feature.selection.metric.threshold,
        criterion = feature.selection.criterion,
        include.covariates = feature.selection.include.covariates,
        fold.ids = inner.fold.ids
      )
      pred.selected = c(covariate.cols, feature.selection.res$pred.selected)
    }
    
    # ---- Build fold-specific columns for metrics & selection ----------------
    metric_colname <- paste0("fold_", fold, "_metric")
    sel_colname    <- paste0("fold_", fold, "_selected")
    
    # initialize vectors with default NA / FALSE
    metric_vec <- rep(NA_real_, length(pred.names.feature.selection))
    sel_vec    <- rep(FALSE, length(pred.names.feature.selection))
    
    # fill metric_vec if feature.selection.res$pred.results exists and has data
    if (!is.null(feature.selection.res) &&
        !is.null(feature.selection.res$pred.results) &&
        nrow(feature.selection.res$pred.results) > 0) {
      res <- feature.selection.res$pred.results
      res$pred <- as.character(res$pred)
      
      # match returned results to our master predictor list
      m <- match(res$pred, pred.names.feature.selection)     # NA for any preds not in pred.names
      valid_res <- !is.na(m)
      
      if (any(valid_res)) {
        metric_vec[m[valid_res]] <- as.numeric(res$metric[valid_res])
      }
    }
    
    # fill selection vector: TRUE if predictor is in pred.selected
    sel_vec <- pred.names.feature.selection %in% pred.selected
    
    # attach columns to result dataframe
    feature_table[[metric_colname]] <- metric_vec
    feature_table[[sel_colname]]   <- sel_vec
    
    # Inner hyperparameter tuning
    best.rmse <- Inf
    best.params <- list(
      eta = eta.values[1],
      max_depth = max.depth.values[1],
      subsample = subsample.values[1],
      colsample = colsample.values[1]
    )
    set.seed(seed + fold)
    
    for (eta.val in eta.values) {
      for (md in max.depth.values) {
        for (ss in subsample.values) {
          for (cs in colsample.values) {
            # --- inner CV (parallel if requested) ---
            if (parallel_backend) {
              inner.rmses <- foreach::foreach(
                ifold = seq_len(n.folds.inner),
                .combine = c,
                .packages = "xgboost"
              ) %dopar% {
                set.seed(seed + fold + ifold)
                inner.train <- df.train[inner.fold.ids != ifold, , drop = FALSE]
                inner.test  <- df.train[inner.fold.ids == ifold, , drop = FALSE]
                dtrain <- xgboost::xgb.DMatrix(data = data.matrix(inner.train[, pred.selected, drop = FALSE]),
                                               label = as.numeric(inner.train[[response.col]]))
                dtest  <- xgboost::xgb.DMatrix(data = data.matrix(inner.test[, pred.selected, drop = FALSE]),
                                               label = as.numeric(inner.test[[response.col]]))
                params <- list(
                  objective = "reg:squarederror",
                  eta = eta.val,
                  max_depth = md,
                  subsample = ss,
                  colsample_bytree = cs,
                  seed = seed + fold + ifold,
                  verbosity = 0
                )
                xgb.fit <- xgboost::xgb.train(
                  params = params,
                  data = dtrain,
                  nrounds = nrounds,
                  verbose = 0,
                  nthread = 1
                )
                preds.inner <- as.numeric(predict(xgb.fit, dtest))
                obs.inner <- as.numeric(inner.test[[response.col]])
                sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
              }
            } else {
              inner.rmses <- numeric(n.folds.inner)
              for (ifold in seq_len(n.folds.inner)) {
                inner.train <- df.train[inner.fold.ids != ifold, , drop = FALSE]
                inner.test  <- df.train[inner.fold.ids == ifold, , drop = FALSE]
                dtrain <- xgboost::xgb.DMatrix(data = data.matrix(inner.train[, pred.selected, drop = FALSE]),
                                               label = as.numeric(inner.train[[response.col]]))
                dtest  <- xgboost::xgb.DMatrix(data = data.matrix(inner.test[, pred.selected, drop = FALSE]),
                                               label = as.numeric(inner.test[[response.col]]))
                params <- list(
                  objective = "reg:squarederror",
                  eta = eta.val,
                  max_depth = md,
                  subsample = ss,
                  colsample_bytree = cs,
                  seed = seed + fold + ifold,
                  verbosity = 0
                )
                xgb.fit <- xgboost::xgb.train(
                  params = params,
                  data = dtrain,
                  nrounds = nrounds,
                  verbose = 0,
                  nthread = 1
                )
                preds.inner <- as.numeric(predict(xgb.fit, dtest))
                obs.inner <- as.numeric(inner.test[[response.col]])
                inner.rmses[ifold] <- sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
              }
            }
            # --------------------------------------------
            mean.rmse <- mean(inner.rmses, na.rm = TRUE)
            if (is.finite(mean.rmse) && mean.rmse < best.rmse) {
              best.rmse <- mean.rmse
              best.params <- list(
                eta = eta.val,
                max.depth = md,
                subsample = ss,
                colsample = cs
              )
            }
          }
        }
      }
    }
    
    # Final training & prediction using selected predictors
    dtrain.full <- xgboost::xgb.DMatrix(data = data.matrix(df.train[, pred.selected, drop = FALSE]), label = as.numeric(df.train[[response.col]]))
    dtest.full  <- xgboost::xgb.DMatrix(data = data.matrix(df.test[, pred.selected, drop = FALSE]), label = as.numeric(df.test[[response.col]]))
    final.params <- list(
      objective = "reg:squarederror",
      eta = best.params$eta,
      max_depth = best.params$max.depth,
      subsample = best.params$subsample,
      colsample_bytree = best.params$colsample,
      seed = seed + fold + 1000,
      verbosity = 0
    )
    final.xgb <- xgboost::xgb.train(
      params = final.params,
      data = dtrain.full,
      nrounds = nrounds,
      verbose = 0,
      nthread = 1
    )
    preds.outer <- as.numeric(predict(final.xgb, dtest.full))
    pred.vec[fold.ids == fold] <- preds.outer
    
    # Extract variable importance from the fitted model (features = pred.selected)
    imp.df <- xgboost::xgb.importance(feature_names = pred.selected, model = final.xgb)
    vi <- setNames(rep(0, length(pred.names)), pred.names)
    if (nrow(imp.df) > 0) {
      # use Gain if present, else Cover/Weight
      if ("Gain" %in% colnames(imp.df))
        vals <- imp.df$Gain
      else if ("Cover" %in% colnames(imp.df))
        vals <- imp.df$Cover
      else
        vals <- imp.df$Frequency
      names(vals) <- imp.df$Feature
      vi[names(vals)] <- vals
    }
    vi.max <- max(vi, na.rm = TRUE)
    if (is.finite(vi.max) &&
        vi.max > 0)
      vi <- vi / vi.max
    else
      vi[] <- NA_real_
    for (v in names(vi))
      if (v %in% rownames(var.imp))
        var.imp[v, fold] <- vi[[v]]
  }
  
  # Identify the fold-specific metric and selection columns
  metric_cols <- grep("^fold_\\d+_metric$", names(feature_table), value = TRUE)
  sel_cols    <- grep("^fold_\\d+_selected$", names(feature_table), value = TRUE)
  
  # Compute mean metric across folds (ignoring NAs)
  feature_table$mean_metric <- rowMeans(feature_table[, metric_cols, drop = FALSE], na.rm = TRUE)
  
  # Count how many folds each predictor was selected in
  feature_table$n_selected <- rowSums(feature_table[, sel_cols, drop = FALSE])
  
  # --- parallel backend teardown ---
  if (parallel_backend) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  # -----------------------------------------
  
  # Observed / predicted and metrics
  observed <- df %>% select(all_of(response.col)) %>% as.matrix()
  predictions <- data.frame("observed" = as.numeric(observed),
                            "predicted" = as.numeric(pred.vec))
  colnames(predictions) <- c("observed", "predicted")
  R2 <- cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman <- cor(predictions$observed,
                    predictions$predicted,
                    method = "spearman",
                    use = "complete.obs")
  RMSE <- sqrt(mean((
    predictions$observed - predictions$predicted
  )^2, na.rm = TRUE))
  sRMSE <- RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # Store metrics
  if (baseline) {
    metrics <- data.frame(
      "data.selection" = "baseline",
      "include.covariates" = feature.selection.include.covariates,
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
      "feature.selection.metric" = "baseline",
      "feature.selection.metric.threshold" = ifelse(feature.selection == "none", 0, feature.selection.metric.threshold),
      "feature.selection.model" = "baseline",
      "feature.selection.criterion" = "baseline",
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "xgboost",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  } else {
    metrics <- data.frame(
      "data.selection" = data.selection,
      "include.covariates" = feature.selection.include.covariates,
      "feature.engineering.col" = feature.engineering.col,
      "feature.engineering.row" = feature.engineering.row,
      "feature.selection" = feature.selection,
      "feature.selection.metric" = ifelse(feature.selection == "none", "none", feature.selection.metric),
      "feature.selection.metric.threshold" = ifelse(feature.selection == "none", 0, feature.selection.metric.threshold),
      "feature.selection.model" = ifelse(feature.selection == "none", "none", feature.selection.model),
      "feature.selection.criterion" = ifelse(feature.selection == "none", "none", feature.selection.criterion),
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "xgboost",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Aggregate variable importance across folds (handle NaN / NA sd)
  meanImp <- rowMeans(var.imp, na.rm = TRUE)
  meanImp[!is.finite(meanImp)] <- 0
  
  sdImp <- apply(var.imp, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s))
      0
    else
      s
  })
  
  varImp <- data.frame(
    var = rownames(var.imp),
    meanImp = meanImp,
    sdImp = sdImp,
    varImp.type = "xgboost_importance",
    stringsAsFactors = FALSE
  )
  rownames(varImp) <- NULL
  
  # Prediction plot and return
  prediction.plot <- cv.plot(pred = predictions$predicted, obs = predictions$observed)
  results <- list(
    predictions = predictions,
    metrics = metrics,
    varImp = varImp,
    feature.selection.metrics = feature_table,
    prediction.plot = prediction.plot
  )
  return(results)
}
