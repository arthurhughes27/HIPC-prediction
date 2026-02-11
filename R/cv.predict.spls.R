cv.predict.spls = function(df,
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
                           alpha.values = seq(0.01, 1, 0.05),
                           nlambda = 5,
                           baseline = FALSE,
                           n.cores = 1,
                           gender.select = NULL) {
  set.seed(seed)
  n <- nrow(df)
  if (is.null(predictor.cols)) {
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  p <- length(pred.names)
  pred.names.feature.selection = pred.names[-which(pred.names %in% covariate.cols)]
  if (!is.null(fold.ids)) {
    n.folds <- length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  pred.vec <- rep(NA_real_, n)
  var.imp <- matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) <- pred.names
  feature_table <- data.frame(pred = pred.names.feature.selection, stringsAsFactors = FALSE)
  
  for (fold in seq_len(n.folds)) {
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    set.seed(seed + fold)
    n.inner <- nrow(df.train)
    inner.fold.ids <- sample(rep(seq_len(n.folds.inner), length.out = n.inner))
    
    if (feature.selection == "none") {
      pred.selected = df.train %>%
        select(-participant_id, -all_of(response.col)) %>%
        colnames()
      feature.selection.res <- NULL
    } else if (feature.selection == "univariate") {
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
    
    metric_colname <- paste0("fold_", fold, "_metric")
    sel_colname    <- paste0("fold_", fold, "_selected")
    metric_vec <- rep(NA_real_, length(pred.names.feature.selection))
    sel_vec    <- rep(FALSE, length(pred.names.feature.selection))
    if (!is.null(feature.selection.res) &&
        !is.null(feature.selection.res$pred.results) &&
        nrow(feature.selection.res$pred.results) > 0) {
      res <- feature.selection.res$pred.results
      res$pred <- as.character(res$pred)
      m <- match(res$pred, pred.names.feature.selection)
      valid_res <- !is.na(m)
      if (any(valid_res)) {
        metric_vec[m[valid_res]] <- as.numeric(res$metric[valid_res])
      }
    }
    sel_vec <- pred.names.feature.selection %in% pred.selected
    feature_table[[metric_colname]] <- metric_vec
    feature_table[[sel_colname]]   <- sel_vec
    
    # ----------------- sPLS hyperparameter tuning & final fit -----------------
    best.params <- list(ncomp = 1, keepX = NULL)
    set.seed(seed + fold)
    
    # prepare data matrices
    x.train <- as.matrix(df.train[, pred.selected, drop = FALSE])
    y.train <- as.matrix(as.numeric(df.train[[response.col]]))
    x.test  <- as.matrix(df.test[, pred.selected, drop = FALSE])
    
    # tuning grids
    max.ncomp <- min(5, ncol(x.train), nrow(x.train) - 1)
    if (max.ncomp < 1) max.ncomp <- 1
    ncomp.grid <- seq_len(max.ncomp)
    # sensible keepX grid: up to p selected vars, at most 6 candidates
    psel <- ncol(x.train)
    keepX.grid <- unique(as.integer(round(seq(1, max(1, psel), length.out = min(6, psel)))))
    
    # parallel backend for tune.spls via cpus argument (mixOmics supports cpus)
    cpus.arg <- ifelse(as.integer(n.cores) > 1, as.integer(n.cores), 1)
    tune.res <- mixOmics::tune.spls(
      x = x.train,
      y = y.train,
      ncomp = ncomp.grid,
      test.keepX = keepX.grid,
      validation = "Mfold",
      folds = n.folds.inner,
      nrepeat = 1,
      measure = "MSE",
      cpus = cpus.arg,
      progressBar = FALSE
    )
    # extract choices
    if (!is.null(tune.res$choice.ncomp)) best.params$ncomp <- tune.res$choice.ncomp else best.params$ncomp <- ncomp.grid[1]
    if (!is.null(tune.res$choice.keepX)) {
      best.keepX <- tune.res$choice.keepX
      # ensure keepX length equals ncomp
      if (length(best.keepX) < best.params$ncomp) {
        best.keepX <- rep(best.keepX[1], best.params$ncomp)
      } else if (length(best.keepX) > best.params$ncomp) {
        best.keepX <- best.keepX[seq_len(best.params$ncomp)]
      }
      best.params$keepX <- best.keepX
    } else {
      best.params$keepX <- rep(min(keepX.grid), best.params$ncomp)
    }
    
    final.fit <- mixOmics::spls(
      x = x.train,
      y = y.train,
      ncomp = best.params$ncomp,
      keepX = best.params$keepX,
      scale = TRUE,
      mode = "regression"
    )
    
    preds.arr <- predict(final.fit, newdata = x.test)$predict
    # preds.arr dims: samples x response_vars x ncomp; take predictions at chosen ncomp
    if (length(dim(preds.arr)) == 3) {
      preds.outer <- as.numeric(preds.arr[, 1, best.params$ncomp])
    } else {
      preds.outer <- as.numeric(preds.arr[, 1])
    }
    pred.vec[fold.ids == fold] <- preds.outer
    
    # variable importance: use sum absolute loadings across components
    loadx <- final.fit$loadings$X
    if (is.null(dim(loadx))) {
      vi <- abs(loadx)
      names(vi) <- colnames(x.train)
    } else {
      vi <- rowSums(abs(loadx[, seq_len(best.params$ncomp), drop = FALSE]), na.rm = TRUE)
      names(vi) <- rownames(loadx)
    }
    vi_max <- max(vi, na.rm = TRUE)
    if (is.finite(vi_max) && vi_max > 0) vi <- vi / vi_max else vi[] <- NA_real_
    for (v in names(vi))
      if (v %in% rownames(var.imp))
        var.imp[v, fold] <- vi[[v]]
  }
  
  metric_cols <- grep("^fold_\\d+_metric$", names(feature_table), value = TRUE)
  sel_cols    <- grep("^fold_\\d+_selected$", names(feature_table), value = TRUE)
  feature_table$mean_metric <- rowMeans(feature_table[, metric_cols, drop = FALSE], na.rm = TRUE)
  feature_table$n_selected <- rowSums(feature_table[, sel_cols, drop = FALSE])
  
  # no extra cluster handling needed (tune.spls used cpus)
  
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
      "gender.select" = ifelse(is.null(gender.select), "none", gender.select),
      "include.covariates" = feature.selection.include.covariates, 
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
      "feature.selection.metric" = "baseline",
      "feature.selection.metric.threshold" = ifelse(feature.selection == "none", 0, feature.selection.metric.threshold),
      "feature.selection.model" = "baseline",
      "feature.selection.criterion" = "baseline",
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "spls",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  } else {
    metrics <- data.frame(
      "data.selection" = data.selection,
      "gender.select" = ifelse(is.null(gender.select), "none", gender.select),
      "include.covariates" = feature.selection.include.covariates, 
      "feature.engineering.col" = feature.engineering.col,
      "feature.engineering.row" = feature.engineering.row,
      "feature.selection" = feature.selection,
      "feature.selection.metric" = ifelse(feature.selection == "none", "none", feature.selection.metric),
      "feature.selection.metric.threshold" = ifelse(feature.selection == "none", 0, feature.selection.metric.threshold),
      "feature.selection.model" = ifelse(feature.selection == "none", "none", feature.selection.model),
      "feature.selection.criterion" = ifelse(feature.selection == "none", "none", feature.selection.criterion),
      "feature.selection.include.covariates" = feature.selection.include.covariates,
      "model" = "spls",
      "R2" = R2,
      "R.spearman" = R.spearman,
      "RMSE" = RMSE,
      "sRMSE" = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  meanImp <- rowMeans(var.imp, na.rm = TRUE)
  meanImp[!is.finite(meanImp)] <- 0
  sdImp <- apply(var.imp, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s)) 0 else s
  })
  varImp <- data.frame(
    var         = rownames(var.imp),
    meanImp     = meanImp,
    sdImp       = sdImp,
    varImp.type = "loading",
    stringsAsFactors = FALSE
  )
  rownames(varImp) <- NULL
  
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