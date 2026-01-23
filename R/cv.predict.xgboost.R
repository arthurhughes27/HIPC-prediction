cv.predict.xgboost = function(df,
                              n.folds = NULL,
                              fold.ids,
                              response.col,
                              predictor.cols = NULL,
                              data.selection,
                              feature.engineering.col,
                              feature.engineering.row,
                              feature.selection,
                              seed = 12345,
                              n.folds.inner = 5,
                              nrounds = 100,
                              eta.values = c(0.01, 0.1),
                              max.depth.values = c(3, 6),
                              subsample.values = 1,
                              colsample.values = 1,
                              baseline = FALSE,
                              n.cores = 1) {
  set.seed(seed)
  n <- nrow(df)
  pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  p <- length(pred.names)
  if (!is.null(fold.ids)) {
    n.folds = length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
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
  
  for (fold in seq_len(n.folds)) {
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    
    best.rmse <- Inf
    best.params <- list(
      eta = eta.values[1],
      max_depth = max.depth.values[1],
      subsample = subsample.values[1],
      colsample = colsample.values[1]
    )
    set.seed(seed + fold)
    n.inner <- nrow(df.train)
    inner.fold.ids <- sample(rep(seq_len(n.folds.inner), length.out = n.inner))
    
    for (eta.val in eta.values) {
      for (md in max.depth.values) {
        for (ss in subsample.values) {
          for (cs in colsample.values) {
            # --- inner CV (parallel if requested) ---
            if (parallel_backend) {
              inner.rmses <- foreach::foreach(ifold = seq_len(n.folds.inner),
                                              .combine = c,
                                              .packages = "xgboost") %dopar% {
                                                set.seed(seed + fold + ifold)
                                                inner.train <- df.train[inner.fold.ids != ifold, , drop = FALSE]
                                                inner.test  <- df.train[inner.fold.ids == ifold, , drop = FALSE]
                                                dtrain <- xgboost::xgb.DMatrix(data = data.matrix(inner.train[, pred.names, drop = FALSE]),
                                                                               label = as.numeric(inner.train[[response.col]]))
                                                dtest  <- xgboost::xgb.DMatrix(data = data.matrix(inner.test[, pred.names, drop = FALSE]),
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
                dtrain <- xgboost::xgb.DMatrix(data = data.matrix(inner.train[, pred.names, drop = FALSE]),
                                               label = as.numeric(inner.train[[response.col]]))
                dtest  <- xgboost::xgb.DMatrix(data = data.matrix(inner.test[, pred.names, drop = FALSE]),
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
    
    dtrain.full <- xgboost::xgb.DMatrix(data = data.matrix(df.train[, pred.names, drop = FALSE]), label = as.numeric(df.train[[response.col]]))
    dtest.full  <- xgboost::xgb.DMatrix(data = data.matrix(df.test[, pred.names, drop = FALSE]), label = as.numeric(df.test[[response.col]]))
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
    
    imp.df <- xgboost::xgb.importance(feature_names = pred.names, model = final.xgb)
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
  
  # --- parallel backend teardown ---
  if (parallel_backend) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  # -----------------------------------------
  
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
  
  if (baseline) {
    metrics <- data.frame(
      "data.selection" = "baseline",
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
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
      "feature.engineering.col" = feature.engineering.col,
      "feature.engineering.row" = feature.engineering.row,
      "feature.selection" = feature.selection,
      "model" = "xgboost",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  varImp <- data.frame(
    var = rownames(var.imp),
    meanImp = rowMeans(var.imp, na.rm = TRUE),
    sdImp = apply(var.imp, 1, sd, na.rm = TRUE),
    varImp.type = "xgboost_importance",
    stringsAsFactors = FALSE
  )
  rownames(varImp) <- NULL
  
  prediction.plot <- cv.plot(pred = predictions$predicted, obs = predictions$observed)
  results <- list(
    predictions = predictions,
    metrics = metrics,
    varImp = varImp,
    prediction.plot = prediction.plot
  )
  return(results)
}
