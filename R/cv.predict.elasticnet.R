cv.predict.elasticnet = function(df,
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
                                 alpha.values = seq(0,1,0.1),
                                 nlambda = 10,
                                 baseline = FALSE) {
  set.seed(seed)
  n <- nrow(df)
  pred.names <- setdiff(colnames(df), c("participant_id", response.col))
  p <- length(pred.names)
  if (!is.null(fold.ids)) {
    n.folds <- length(unique(fold.ids))
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  pred.vec <- rep(NA_real_, n)
  var.imp <- matrix(NA_real_, nrow = p, ncol = n.folds)
  rownames(var.imp) <- pred.names
  
  for (fold in seq_len(n.folds)) {
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    
    best.rmse <- Inf
    best.params <- list(alpha = alpha.values[1], lambda = NULL)
    set.seed(seed + fold)
    
    for (alpha.val in alpha.values) {
      x.train <- data.matrix(df.train[, pred.names, drop = FALSE])
      y.train <- as.numeric(df.train[[response.col]])
      cvfit <- glmnet::cv.glmnet(x = x.train,
                                 y = y.train,
                                 alpha = alpha.val,
                                 nfolds = n.folds.inner,
                                 nlambda = nlambda,
                                 type.measure = "mse",
                                 standardize = TRUE,
                                 intercept = TRUE)
      idx <- which.min(abs(cvfit$lambda - cvfit$lambda.min))
      mse <- cvfit$cvm[idx]
      rmse <- sqrt(mse)
      if (is.finite(rmse) && rmse < best.rmse) {
        best.rmse <- rmse
        best.params$alpha <- alpha.val
        best.params$lambda <- cvfit$lambda.min
      }
    }
    
    x.train.full <- data.matrix(df.train[, pred.names, drop = FALSE])
    y.train.full <- as.numeric(df.train[[response.col]])
    x.test.full  <- data.matrix(df.test[, pred.names, drop = FALSE])
    
    final.fit <- glmnet::glmnet(x = x.train.full,
                                y = y.train.full,
                                alpha = best.params$alpha,
                                lambda = best.params$lambda,
                                standardize = TRUE,
                                intercept = TRUE,
                                nlambda = 1)
    
    preds.outer <- as.numeric(predict(final.fit, newx = x.test.full, s = best.params$lambda))
    pred.vec[fold.ids == fold] <- preds.outer
    
    coefs <- as.matrix(coef(final.fit, s = best.params$lambda))
    # remove intercept
    if ("(Intercept)" %in% rownames(coefs)) {
      coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
    }
    vi <- abs(as.numeric(coefs))
    names(vi) <- rownames(coefs)
    vi_max <- max(vi, na.rm = TRUE)
    if (is.finite(vi_max) && vi_max > 0) vi <- vi / vi_max else vi[] <- NA_real_
    for (v in names(vi)) if (v %in% rownames(var.imp)) var.imp[v, fold] <- vi[[v]]
  }
  
  observed <- df %>% select(all_of(response.col)) %>% as.matrix()
  predictions <- data.frame("observed" = as.numeric(observed), "predicted" = as.numeric(pred.vec))
  colnames(predictions) <- c("observed", "predicted")
  R2 <- cor(predictions$observed, predictions$predicted, use = "complete.obs")^2
  R.spearman <- cor(predictions$observed, predictions$predicted, method = "spearman", use = "complete.obs")
  RMSE <- sqrt(mean((predictions$observed - predictions$predicted)^2, na.rm = TRUE))
  sRMSE <- RMSE / sd(predictions$observed, na.rm = TRUE)
  
  # Store metrics
  if (baseline) {
    metrics <- data.frame(
      "data.selection" = "baseline",
      "feature.engineering.col" = "baseline",
      "feature.engineering.row" = "baseline",
      "feature.selection" = "baseline",
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
      "feature.engineering.col" = feature.engineering.col,
      "feature.engineering.row" = feature.engineering.row,
      "feature.selection" = feature.selection,
      "model" = "elasticnet",
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
    varImp.type = "elasticnet_coefficient",
    stringsAsFactors = FALSE
  )
  rownames(varImp) <- NULL
  
  prediction.plot <- cv.plot(pred = predictions$predicted, obs = predictions$observed)
  results <- list(predictions = predictions, metrics = metrics, varImp = varImp, prediction.plot = prediction.plot)
  return(results)
}
