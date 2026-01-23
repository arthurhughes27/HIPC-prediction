cv.predict.ranger = function(df,
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
                             num.trees = 500,
                             mtry.values = NULL,
                             min.node.size.values = NULL,
                             baseline = FALSE,
                             n.cores = 1) {        
  # Set seed
  set.seed(seed)
  
  # Number of observations
  n = nrow(df)
  
  # Predictor names
  pred.names = setdiff(colnames(df), c("participant_id", response.col))
  
  # Number of predictors
  p = length(pred.names)
  
  # Outer folds
  if (!is.null(fold.ids)) {
    n.folds = length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # defaults for tuning grid if not provided
  if (is.null(mtry.values)) {
    # candidate mtry values: small, classical, and larger fractions of p
    vals <- unique(floor(c(
      1,
      sqrt(p),
      p / 5,
      p / 3,
      p / 2
    )))
    
    # keep valid range
    vals <- vals[vals >= 1 & vals <= p]
    
    if (length(vals) == 0) {
      vals <- min(1, p)
    }
    
    mtry.values <- sort(unique(vals))
  } else {
    # user-supplied grid: keep only valid values
    mtry.values <- mtry.values[mtry.values >= 1 & mtry.values <= p]
  }
  
  # ---- defaults for min.node.size grid ----
  if (is.null(min.node.size.values)) {
    # broader regularisation grid
    min.node.size.values <- c(1, 3, 5, 10, 20, 50)
  }
  
  # keep only valid values
  min.node.size.values <- min.node.size.values[min.node.size.values >= 1]
  
  if (length(min.node.size.values) == 0) {
    min.node.size.values <- 1
  }
  
  # ---------------------------------------------------------------------------
  # Set up parallel backend if requested (Option B). We create a PSOCK cluster
  # with exactly n.cores workers when n.cores > 1. We will ensure ranger inside
  # each worker uses a single thread to avoid oversubscription.
  # ---------------------------------------------------------------------------
  parallel_backend <- FALSE
  if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1) {
    stop("n.cores must be a positive integer of length 1")
  }
  n.cores <- as.integer(n.cores)
  if (n.cores > 1) {
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for parallel execution (install them).")
    }
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    parallel_backend <- TRUE
  }
  # ---------------------------------------------------------------------------
  
  # Initialise results vectors
  pred.vec = rep(NA_real_, n)
  var.imp = matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) = pred.names
  
  for (fold in 1:n.folds) {
    # outer CV
    df.train = df[fold.ids != fold, , drop = FALSE] # Training data
    df.test  = df[fold.ids == fold, , drop = FALSE] # Testing data
    
    # inner CV tuning on df.train
    best.rmse = Inf # Best root mean squared error value
    best.params = list(mtry = mtry.values[1], min.node.size = min.node.size.values[1]) # Parameters producing best RMSE
    
    # create inner fold ids
    set.seed(seed + fold) # vary per outer fold
    n.inner = nrow(df.train) # number of observations in training data
    inner.fold.ids = sample(rep(seq_len(n.folds.inner), length.out = n.inner)) # Fold ids
    
    for (mtry.val in mtry.values) {
      # For each mtry value
      for (min.node in min.node.size.values) {
        # For each min node value
        # --- compute inner.rmses either in parallel (foreach) or sequentially ---
        if (parallel_backend) {
          # parallel across inner folds; ensure ranger inside each worker uses a single thread
          inner.rmses <- foreach::foreach(ifold = seq_len(n.folds.inner), .combine = c, .packages = "ranger") %dopar% {
            # set a reproducible-ish seed per task (worker RNG is separate)
            set.seed(seed + fold + ifold)
            inner.train = df.train[inner.fold.ids != ifold, , drop = FALSE]
            inner.test  = df.train[inner.fold.ids == ifold, , drop = FALSE]
            # prepare data for ranger: include response and predictors only
            inner.train.rf = inner.train[, c(pred.names, response.col), drop = FALSE]
            inner.test.rf  = inner.test[, c(pred.names, response.col), drop = FALSE]
            # fit ranger (num.threads = 1 inside workers to avoid oversubscription)
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
            # Make predictions on inner test data
            preds.inner = as.numeric(predict(rf.fit, data = inner.test.rf)$predictions)
            # Inner try values
            obs.inner = as.numeric(inner.test.rf[[response.col]])
            # Inner RMSE
            sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
          }
        } else {
          # sequential computation (original approach)
          inner.rmses = numeric(n.folds.inner)
          for (ifold in seq_len(n.folds.inner)) {
            inner.train = df.train[inner.fold.ids != ifold, , drop = FALSE]
            inner.test  = df.train[inner.fold.ids == ifold, , drop = FALSE]
            # prepare data for ranger: include response and predictors only
            inner.train.rf = inner.train[, c(pred.names, response.col), drop = FALSE]
            inner.test.rf  = inner.test[, c(pred.names, response.col), drop = FALSE]
            # fit ranger
            rf.fit = ranger::ranger(
              dependent.variable.name = response.col,
              data = inner.train.rf,
              num.trees = num.trees,
              mtry = mtry.val,
              min.node.size = min.node,
              importance = "permutation",
              write.forest = TRUE,
              seed = seed + fold + ifold,
              num.threads = 1   # keep single-threaded for consistency; change if desired
            )
            # Make predictions on inner test data
            preds.inner = as.numeric(predict(rf.fit, data = inner.test.rf)$predictions)
            # Inner try values
            obs.inner = as.numeric(inner.test.rf[[response.col]])
            # Inner RMSEs
            inner.rmses[ifold] = sqrt(mean((obs.inner - preds.inner)^2, na.rm = TRUE))
          }
        }
        
        # Mean RMSE across folds
        mean.rmse = mean(inner.rmses, na.rm = TRUE)
        if (is.finite(mean.rmse) &&
            mean.rmse < best.rmse) {
          # Find best parameter
          best.rmse = mean.rmse
          best.params = list(mtry = mtry.val, min.node.size = min.node)
        }
      }
    }
    
    # Fit final model using best parameter values
    train.rf = df.train[, c(pred.names, response.col), drop = FALSE]
    test.rf  = df.test[, c(pred.names, response.col), drop = FALSE]
    # Final model (keep ranger single-threaded to avoid interference with cluster workers)
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
    
    # Predict on outer test
    preds.outer = as.numeric(predict(final.rf, data = test.rf)$predictions)
    pred.vec[fold.ids == fold] = preds.outer
    
    # store permutation standardised variable importance
    vi = final.rf$variable.importance / max(final.rf$variable.importance)
    for (v in names(vi)) {
      if (v %in% rownames(var.imp))
        var.imp[v, fold] = vi[[v]]
    }
  }
  
  # If we started a cluster, stop it and unregister
  if (parallel_backend) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  
  # Extract observed values
  observed = df %>% select(all_of(response.col)) %>% as.matrix()
  # Save predicted vs observed values
  predictions = data.frame("observed" = as.numeric(observed),
                           "predicted" = as.numeric(pred.vec))
  colnames(predictions) = c("observed", "predicted")
  
  # Compute metrics
  R2 = cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman = cor(predictions$observed,
                   predictions$predicted,
                   method = "spearman",
                   use = "complete.obs")
  RMSE = sqrt(mean((
    predictions$observed - predictions$predicted
  )^2, na.rm = TRUE))
  sRMSE = RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # Store metrics
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
  
  # Store mean variable importance
  varImp = data.frame(
    var = rownames(var.imp),
    meanImp = rowMeans(var.imp, na.rm = TRUE),
    sdImp = apply(var.imp, 1, sd, na.rm = TRUE),
    varImp.type = "permutation",
    stringsAsFactors = FALSE
  )
  rownames(varImp) = NULL
  
  # Prediction plot
  prediction.plot = cv.plot(pred = predictions$predicted, obs = predictions$observed)
  
  # Return all results
  results = list(
    predictions = predictions,
    metrics = metrics,
    varImp = varImp,
    prediction.plot = prediction.plot
  )
  return(results)
}
