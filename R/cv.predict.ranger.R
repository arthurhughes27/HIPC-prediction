cv.predict.ranger = function(df,
                             n.folds = NULL,
                             fold.ids,
                             response.col,
                             predictor.cols = NULL,
                             data.selection,
                             feature.engineering.col,
                             feature.engineering.row,
                             feature.selection = "none",
                             feature.selection.metric = "sRMSE",
                             feature.selection.metric.threshold = 1,
                             feature.selection.model = "lm",
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
  
  if (is.null(predictor.cols)){
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  p = length(pred.names)
  
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
    if (length(vals) == 0) { vals <- min(1, p) }
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
  
  # -------------------------
  # Outer CV loop
  # -------------------------
  for (fold in 1:n.folds) {
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
    if (feature.selection == "none") {
      # No selection: keep all columns except participant_id
      pred.selected = df.train %>%
        select(-participant_id, -all_of(response.col)) %>%
        colnames()
      
    } else if (feature.selection == "univariate") {
      # Univariate feature selection using inner CV
      feature.selection.res = feature.selection.univariate(
        df = df.train,
        response.col,
        covariate.cols,
        model = feature.selection.model,
        metric = feature.selection.metric,
        metric.threshold = feature.selection.metric.threshold,
        fold.ids = inner.fold.ids
      )
      pred.selected = c(covariate.cols, feature.selection.res$pred.selected)
    }
    
    # Update possible mtry values
    mtry.values <- mtry.values[mtry.values >= 1 & mtry.values <= length(pred.selected)]
    
    # -------------------------------------------------------------------------
    # Inner hyperparameter search (over mtry.values and min.node.size.values)
    # -------------------------------------------------------------------------
    for (mtry.val in mtry.values) {
      for (min.node in min.node.size.values) {
        # compute inner.rmses either in parallel (foreach) or sequentially
        if (parallel_backend) {
          # parallel across inner folds; ensure ranger inside each worker uses a single thread
          inner.rmses <- foreach::foreach(ifold = seq_len(n.folds.inner),
                                          .combine = c,
                                          .packages = "ranger") %dopar% {
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
  predictions = data.frame(
    "observed"  = as.numeric(observed),
    "predicted" = as.numeric(pred.vec)
  )
  colnames(predictions) = c("observed", "predicted")
  
  R2 = cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman = cor(predictions$observed,
                   predictions$predicted,
                   method = "spearman",
                   use = "complete.obs")
  RMSE = sqrt(mean((predictions$observed - predictions$predicted)^2, na.rm = TRUE))
  sRMSE = RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # -------------------------
  # Store metrics (baseline vs metadata)
  # -------------------------
  if (baseline) {
    metrics <- data.frame(
      "data.selection" = "baseline",
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
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
      "feature.engineering.col" = feature.engineering.col,
      "feature.engineering.row" = feature.engineering.row,
      "feature.selection" = feature.selection,
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
    if (!is.finite(s)) 0 else s
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
    prediction.plot = prediction.plot
  )
  return(results)
}
