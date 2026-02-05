cv.predict.ranger = function(df,
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
                             num.trees = 500,
                             mtry.values = NULL,
                             min.node.size.values = NULL,
                             baseline = FALSE,
                             n.cores = 1) {
  # -------------------------
  # Reproducibility
  # -------------------------
  set.seed(seed)
  
  # -------------------------
  # Basic dimensions / predictors
  # -------------------------
  n = nrow(df)
  
  if (is.null(predictor.cols)) {
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  p = length(pred.names)
  
  pred.names.feature.selection = pred.names[-which(pred.names %in% covariate.cols)]
  
  # -------------------------
  # Outer fold ids
  # -------------------------
  if (!is.null(fold.ids)) {
    n.folds = length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # -------------------------
  # Defaults for tuning grids
  # -------------------------
  if (is.null(mtry.values)) {
    vals <- unique(floor(c(1, sqrt(p), p / 5, p / 3, p / 2)))
    vals <- vals[vals >= 1 & vals <= p]
    if (length(vals) == 0) {
      vals <- min(1, p)
    }
    mtry.values <- sort(unique(vals))
  } else {
    mtry.values <- mtry.values[mtry.values >= 1 & mtry.values <= p]
  }
  
  if (is.null(min.node.size.values)) {
    min.node.size.values <- c(1, 3, 5, 10, 20, 50)
  }
  min.node.size.values <- min.node.size.values[min.node.size.values >= 1]
  if (length(min.node.size.values) == 0) {
    min.node.size.values <- 1
  }
  
  # ---------------------------------------------------------------------------
  # Set up parallel backend if requested (PSOCK cluster used for foreach)
  # ---------------------------------------------------------------------------
  parallel_backend <- FALSE
  n.cores <- as.integer(n.cores)
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    parallel_backend <- TRUE
  }
  # ---------------------------------------------------------------------------
  
  # -------------------------
  # Initialise result containers
  # -------------------------
  pred.vec = rep(NA_real_, n)
  var.imp = matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) = pred.names
  
  # Container for variable selection
  feature_table <- data.frame(pred = pred.names.feature.selection, stringsAsFactors = FALSE)
  
  # -------------------------
  # Outer CV loop
  # -------------------------
  for (fold in seq_len(n.folds)) {
    # outer CV: training / testing split
    df.train = df[fold.ids != fold, , drop = FALSE]
    df.test  = df[fold.ids == fold, , drop = FALSE]
    
    # inner CV tuning: initialise best trackers
    best.rmse = Inf
    best.params = list(mtry = mtry.values[1], min.node.size = min.node.size.values[1])
    
    # create inner fold ids (reproducible per outer fold)
    set.seed(seed + fold)
    n.inner = nrow(df.train)
    inner.fold.ids = sample(rep(seq_len(n.folds.inner), length.out = n.inner))
    
    ## Feature selection (performed inside each outer fold)
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
    
    # Update possible mtry values
    mtry.values <- mtry.values[mtry.values >= 1 &
                                 mtry.values <= length(pred.selected)]
    
    # -------------------------------------------------------------------------
    # Inner hyperparameter search (over mtry.values and min.node.size.values)
    # -------------------------------------------------------------------------
    for (mtry.val in mtry.values) {
      for (min.node in min.node.size.values) {
        # compute inner.rmses either in parallel (foreach) or sequentially
        if (parallel_backend) {
          # parallel across inner folds; ensure ranger inside each worker uses a single thread
          inner.rmses <- foreach::foreach(
            ifold = seq_len(n.folds.inner),
            .combine = c,
            .packages = "ranger"
          ) %dopar% {
            set.seed(seed + fold + ifold)
            inner.train = df.train[inner.fold.ids != ifold, , drop = FALSE]
            inner.test  = df.train[inner.fold.ids == ifold, , drop = FALSE]
            # prepare data using selected predictors + response
            inner.train.rf = inner.train[, c(pred.selected, response.col), drop = FALSE]
            inner.test.rf  = inner.test[, c(pred.selected, response.col), drop = FALSE]
            rf.fit = ranger::ranger(
              dependent.variable.name = response.col,
              data = inner.train.rf,
              num.trees = num.trees,
              mtry = mtry.val,
              min.node.size = min.node,
              importance = "permutation",
              write.forest = TRUE,
              seed = seed + fold + ifold,
              num.threads = 1
            )
            preds.inner = as.numeric(predict(rf.fit, data = inner.test.rf)$predictions)
            obs.inner = as.numeric(inner.test.rf[[response.col]])
            sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
          }
        } else {
          # sequential computation
          inner.rmses = numeric(n.folds.inner)
          for (ifold in seq_len(n.folds.inner)) {
            inner.train = df.train[inner.fold.ids != ifold, , drop = FALSE]
            inner.test  = df.train[inner.fold.ids == ifold, , drop = FALSE]
            inner.train.rf = inner.train[, c(pred.selected, response.col), drop = FALSE]
            inner.test.rf  = inner.test[, c(pred.selected, response.col), drop = FALSE]
            rf.fit = ranger::ranger(
              dependent.variable.name = response.col,
              data = inner.train.rf,
              num.trees = num.trees,
              mtry = mtry.val,
              min.node.size = min.node,
              importance = "permutation",
              write.forest = TRUE,
              seed = seed + fold + ifold,
              num.threads = 1
            )
            preds.inner = as.numeric(predict(rf.fit, data = inner.test.rf)$predictions)
            obs.inner = as.numeric(inner.test.rf[[response.col]])
            inner.rmses[ifold] = sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
          }
        }
        
        # Evaluate mean RMSE across inner folds and update best parameters
        mean.rmse = mean(inner.rmses, na.rm = TRUE)
        if (is.finite(mean.rmse) && mean.rmse < best.rmse) {
          best.rmse = mean.rmse
          best.params = list(mtry = mtry.val, min.node.size = min.node)
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Fit final model on df.train using best parameters and predict on df.test
    # -------------------------------------------------------------------------
    train.rf = df.train[, c(pred.selected, response.col), drop = FALSE]
    test.rf  = df.test[, c(pred.selected, response.col), drop = FALSE]
    final.rf = ranger::ranger(
      dependent.variable.name = response.col,
      data = train.rf,
      num.trees = num.trees,
      mtry = best.params$mtry,
      min.node.size = best.params$min.node.size,
      importance = "permutation",
      write.forest = TRUE,
      seed = seed + fold + 1000,
      num.threads = 1
    )
    
    # Predict on outer test and store
    preds.outer = as.numeric(predict(final.rf, data = test.rf)$predictions)
    pred.vec[fold.ids == fold] = preds.outer
    
    # store permutation standardised variable importance for this fold
    vi = final.rf$variable.importance / max(final.rf$variable.importance)
    for (v in names(vi)) {
      if (v %in% rownames(var.imp))
        var.imp[v, fold] = vi[[v]]
    }
  }
  
  # Identify the fold-specific metric and selection columns
  metric_cols <- grep("^fold_\\d+_metric$", names(feature_table), value = TRUE)
  sel_cols    <- grep("^fold_\\d+_selected$", names(feature_table), value = TRUE)
  
  # Compute mean metric across folds (ignoring NAs)
  feature_table$mean_metric <- rowMeans(feature_table[, metric_cols, drop = FALSE], na.rm = TRUE)
  
  # Count how many folds each predictor was selected in
  feature_table$n_selected <- rowSums(feature_table[, sel_cols, drop = FALSE])
  
  # ---------------------------------------------------------------------------
  # Stop cluster if started and unregister foreach backend
  # ---------------------------------------------------------------------------
  if (parallel_backend) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  # ---------------------------------------------------------------------------
  
  # -------------------------
  # Observed / predicted and metrics
  # -------------------------
  observed = df %>% select(all_of(response.col)) %>% as.matrix()
  predictions = data.frame("observed"  = as.numeric(observed),
                           "predicted" = as.numeric(pred.vec))
  colnames(predictions) = c("observed", "predicted")
  
  R2 = cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman = cor(predictions$observed,
                   predictions$predicted,
                   method = "spearman",
                   use = "complete.obs")
  RMSE = sqrt(mean((
    predictions$observed - predictions$predicted
  )^2, na.rm = TRUE))
  sRMSE = RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # -------------------------
  # Store metrics (baseline vs metadata)
  # -------------------------
  if (baseline) {
    metrics <- data.frame(
      "data.selection" = "baseline",
      "include.covariates" = feature.selection.include.covariates,
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
      "feature.selection.metric" = "baseline",
      "feature.selection.metric.threshold" = 0,
      "feature.selection.model" = "baseline",
      "feature.selection.criterion" = "baseline",
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "ranger",
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
      "feature.selection.metric.threshold" = ifelse(feature.selection == "none", "none", feature.selection.metric.threshold),
      "feature.selection.model" = ifelse(feature.selection == "none", "none", feature.selection.model),
      "feature.selection.criterion" = ifelse(feature.selection == "none", "none", feature.selection.criterion),
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "ranger",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  # -------------------------
  # Aggregate variable importance across folds (fix NaNs / NA sd)
  # -------------------------
  meanImp <- rowMeans(var.imp, na.rm = TRUE)
  meanImp[!is.finite(meanImp)] <- 0
  
  sdImp <- apply(var.imp, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s))
      0
    else
      s
  })
  
  varImp = data.frame(
    var = rownames(var.imp),
    meanImp = meanImp,
    sdImp = sdImp,
    varImp.type = "permutation",
    stringsAsFactors = FALSE
  )
  rownames(varImp) = NULL
  
  # -------------------------
  # Prediction plot and return
  # -------------------------
  prediction.plot = cv.plot(pred = predictions$predicted, obs = predictions$observed)
  
  results = list(
    predictions = predictions,
    metrics = metrics,
    varImp = varImp,
    feature.selection.metrics = feature_table,
    prediction.plot = prediction.plot
  )
  return(results)
}
