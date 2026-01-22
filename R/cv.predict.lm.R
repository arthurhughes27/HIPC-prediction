cv.predict.lm = function(df,
                         n.folds = NULL,
                         fold.ids,
                         response.col,
                         predictor.cols = NULL,
                         data.selection,
                         feature.engineering.col,
                         feature.engineering.row,
                         feature.selection,
                         seed = 12345,
                         baseline = FALSE) {
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
  for (fold in seq_len(n.folds)) {
    df.train <- df[fold.ids != fold, , drop = FALSE]
    df.test  <- df[fold.ids == fold, , drop = FALSE]
    predictor_names <- intersect(pred.names, colnames(df.train))
    predictor.cols.safe <- paste0("`", predictor_names, "`")
    form <- as.formula(paste0(
      "`",
      response.col,
      "` ~ ",
      paste(predictor.cols.safe, collapse = " + ")
    ))
    fit <- lm(form, data = df.train)
    preds.outer <- predict(fit, newdata = df.test)
    pred.vec[fold.ids == fold] <- as.numeric(preds.outer)
    coefs <- coef(fit)
    coefs <- coefs[names(coefs) != "(Intercept)"]
    var.names <- rownames(var.imp)
    for (v in names(coefs)) {
      v.clean <- gsub("`", "", v)
      if (v.clean %in% var.names)
        var.imp[v.clean, fold] <- coefs[v]
    }
  }
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
      "model" = "lm",
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
      "model" = "lm",
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
    varImp.type = "coefficient",
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
