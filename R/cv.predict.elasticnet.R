cv.predict.elasticnet = function(df,
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
                                 alpha.values = seq(0.01, 1, 0.1),
                                 nlambda = 5,
                                 baseline = FALSE,
                                 n.cores = 1) {
  # Reproducibility
  set.seed(seed)
  
  # Number of observations
  n <- nrow(df)
  
  # Predictor names: either provided or derived from df
  if (is.null(predictor.cols)) {
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  p <- length(pred.names)
  
  pred.names.feature.selection = pred.names[-which(pred.names %in% covariate.cols)]
  
  # Outer fold ids: if provided derive n.folds; otherwise sample fold ids
  if (!is.null(fold.ids)) {
    n.folds <- length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # Containers for predictions and per-fold variable importances
  pred.vec <- rep(NA_real_, n)
  var.imp <- matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) <- pred.names
  
  # Container for variable selection
  feature_table <- data.frame(pred = pred.names.feature.selection, stringsAsFactors = FALSE)
  
  # Outer cross-validation loop
  for (fold in seq_len(n.folds)) {
    # Split training / test for this outer fold
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    
    # create inner fold ids (used by feature selection if requested)
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
    
    
    # -------------------------------------------------------------------------
    # Hyperparameter tuning across alpha values (inner CV performed by glmnet)
    # -------------------------------------------------------------------------
    best.rmse <- Inf
    best.params <- list(alpha = alpha.values[1], lambda = NULL)
    set.seed(seed + fold)
    
    # ---------------------------------------------------------------------------
    # Parallel backend setup (PSOCK cluster registered for foreach if n.cores > 1)
    # ---------------------------------------------------------------------------
    parallel_backend <- FALSE
    n.cores <- as.integer(n.cores)
    if (n.cores > 1) {
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
      parallel_backend <- TRUE
    }
    
    if (parallel_backend) {
      res_alpha <- foreach::foreach(
        i = seq_along(alpha.values),
        .combine = rbind,
        .packages = "glmnet"
      ) %dopar% {
        alpha.val <- alpha.values[i]
        set.seed(seed + fold + i)
        x.train <- data.matrix(df.train[, pred.selected, drop = FALSE])
        y.train <- as.numeric(df.train[[response.col]])
        cvfit <- glmnet::cv.glmnet(
          x = x.train,
          y = y.train,
          alpha = alpha.val,
          nfolds = n.folds.inner,
          nlambda = nlambda,
          type.measure = "mse",
          standardize = TRUE,
          intercept = TRUE
        )
        idx <- which.min(abs(cvfit$lambda - cvfit$lambda.min))
        mse <- cvfit$cvm[idx]
        rmse <- sqrt(mse)
        c(
          alpha = alpha.val,
          rmse = rmse,
          lambda = cvfit$lambda.min
        )
      }
      if (is.matrix(res_alpha) && nrow(res_alpha) > 0) {
        for (r in seq_len(nrow(res_alpha))) {
          alpha.val <- as.numeric(res_alpha[r, "alpha"])
          rmse <- as.numeric(res_alpha[r, "rmse"])
          lambda.val <- as.numeric(res_alpha[r, "lambda"])
          if (is.finite(rmse) && rmse < best.rmse) {
            best.rmse <- rmse
            best.params$alpha <- alpha.val
            best.params$lambda <- lambda.val
          }
        }
      }
    } else {
      for (alpha.val in alpha.values) {
        x.train <- data.matrix(df.train[, pred.selected, drop = FALSE])
        y.train <- as.numeric(df.train[[response.col]])
        cvfit <- glmnet::cv.glmnet(
          x = x.train,
          y = y.train,
          alpha = alpha.val,
          nfolds = n.folds.inner,
          nlambda = nlambda,
          type.measure = "mse",
          standardize = TRUE,
          intercept = TRUE
        )
        idx <- which.min(abs(cvfit$lambda - cvfit$lambda.min))
        mse <- cvfit$cvm[idx]
        rmse <- sqrt(mse)
        if (is.finite(rmse) && rmse < best.rmse) {
          best.rmse <- rmse
          best.params$alpha <- alpha.val
          best.params$lambda <- cvfit$lambda.min
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Final model fit on df.train using best params and predict on outer test
    # -------------------------------------------------------------------------
    x.train.full <- data.matrix(df.train[, pred.selected, drop = FALSE])
    y.train.full <- as.numeric(df.train[[response.col]])
    x.test.full  <- data.matrix(df.test[, pred.selected, drop = FALSE])
    
    final.fit <- glmnet::glmnet(
      x = x.train.full,
      y = y.train.full,
      alpha = best.params$alpha,
      lambda = best.params$lambda,
      standardize = TRUE,
      intercept = TRUE,
      nlambda = 1
    )
    
    preds.outer <- as.numeric(predict(final.fit, newx = x.test.full, s = best.params$lambda))
    pred.vec[fold.ids == fold] <- preds.outer
    
    # Extract coefficients (absolute, normalized) and store per-fold variable importance
    coefs <- as.matrix(coef(final.fit, s = best.params$lambda))
    # remove intercept
    if ("(Intercept)" %in% rownames(coefs)) {
      coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
    }
    vi <- abs(as.numeric(coefs))
    names(vi) <- rownames(coefs)
    vi_max <- max(vi, na.rm = TRUE)
    if (is.finite(vi_max) &&
        vi_max > 0)
      vi <- vi / vi_max
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
  
  # ---------------------------------------------------------------------------
  # Parallel backend teardown
  # ---------------------------------------------------------------------------
  if (parallel_backend) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  
  # ---------------------------------------------------------------------------
  # Compute performance metrics and aggregate variable importance across folds
  # ---------------------------------------------------------------------------
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
      "model" = "elasticnet",
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
      "model" = "elasticnet",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Aggregate variable importance across folds:
  # - meanImp: mean importance across folds (variables never selected -> 0)
  # - sdImp:  standard deviation across folds (insufficient values -> 0)
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
    var         = rownames(var.imp),
    meanImp     = meanImp,
    sdImp       = sdImp,
    varImp.type = "coefficient",
    stringsAsFactors = FALSE
  )
  rownames(varImp) <- NULL
  
  # Prediction plot and return object
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
