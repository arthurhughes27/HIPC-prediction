cv.predict.lm = function(df,
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
                         n.folds.inner = 5,
                         seed = 12345,
                         baseline = FALSE) {
  # Reproducibility
  set.seed(seed)
  
  # Number of observations
  n <- nrow(df)
  
  # Determine predictor names (either provided or derived from df)
  if (is.null(predictor.cols)) {
    pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  } else {
    pred.names = predictor.cols
  }
  
  # Number of predictors
  p <- length(pred.names)
  
  # Outer fold ids: if provided, derive n.folds; otherwise sample fold ids
  if (!is.null(fold.ids)) {
    n.folds = length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # Containers for predictions and per-fold variable importances (coefficients)
  pred.vec <- rep(NA_real_, n)
  var.imp <- matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) <- pred.names
  
  # Outer cross-validation loop
  for (fold in seq_len(n.folds)) {
    # Split training / test for this outer fold
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    
    # Create inner fold ids for feature-selection (reproducible per outer fold)
    set.seed(seed + fold) # vary per outer fold
    n.inner = nrow(df.train) # number of observations in training data
    inner.fold.ids = sample(rep(seq_len(n.folds.inner), length.out = n.inner)) # Fold ids
    
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
    
    # Build formula safely quoting predictor names to allow non-syntactic names
    predictor.cols.safe <- paste0("`", pred.selected, "`")
    form <- as.formula(
      paste0(
        "`", response.col, "` ~ ",
        paste(predictor.cols.safe, collapse = " + ")
      )
    )
    
    # Fit linear model on training data for this outer fold
    fit <- lm(form, data = df.train)
    
    # Predict on outer test data and store
    preds.outer <- predict(fit, newdata = df.test)
    pred.vec[fold.ids == fold] <- as.numeric(preds.outer)
    
    # Extract coefficients (excluding intercept) and record into var.imp matrix
    coefs <- coef(fit)
    coefs <- coefs[names(coefs) != "(Intercept)"]
    var.names <- rownames(var.imp)
    for (v in names(coefs)) {
      v.clean <- gsub("`", "", v)
      if (v.clean %in% var.names)
        var.imp[v.clean, fold] <- coefs[v]
    }
  }
  
  # Observed / predicted summary
  observed <- df %>% select(all_of(response.col)) %>% as.matrix()
  predictions <- data.frame(
    "observed"  = as.numeric(observed),
    "predicted" = as.numeric(pred.vec)
  )
  colnames(predictions) <- c("observed", "predicted")
  
  # Performance metrics
  R2 <- cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman <- cor(predictions$observed,
                    predictions$predicted,
                    method = "spearman",
                    use = "complete.obs")
  RMSE <- sqrt(mean((predictions$observed - predictions$predicted)^2, na.rm = TRUE))
  sRMSE <- RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # Store metrics (baseline vs supplied metadata)
  if (baseline) {
    metrics <- data.frame(
      "data.selection"             = "baseline",
      "feature.engineering.col"    = "baseline",
      "feature.engineering.row"    = "baseline",
      "feature.selection"          = "baseline",
      "model"                      = "lm",
      "R2"                         = R2,
      "R.spearman"                 = R.spearman,
      "RMSE"                       = RMSE,
      "sRMSE"                      = sRMSE,
      stringsAsFactors = FALSE
    )
  } else {
    metrics <- data.frame(
      "data.selection"             = data.selection,
      "feature.engineering.col"    = feature.engineering.col,
      "feature.engineering.row"    = feature.engineering.row,
      "feature.selection"          = feature.selection,
      "model"                      = "lm",
      "R2"                         = R2,
      "R.spearman"                 = R.spearman,
      "RMSE"                       = RMSE,
      "sRMSE"                      = sRMSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Aggregate variable importance across folds:
  # - meanImp: mean coefficient across folds (treat variables never selected as 0)
  # - sdImp:  standard deviation across folds (treat insufficient values as 0)
  meanImp <- rowMeans(var.imp, na.rm = TRUE)
  # when a variable was never selected rowMeans(...) gives NaN -> set to 0
  meanImp[!is.finite(meanImp)] <- 0
  
  sdImp <- apply(var.imp, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    # sd() returns NA if fewer than 2 non-NA values; treat that as 0
    if (!is.finite(s)) 0 else s
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
    predictions     = predictions,
    metrics         = metrics,
    varImp          = varImp,
    prediction.plot = prediction.plot
  )
  
  return(results)
}
