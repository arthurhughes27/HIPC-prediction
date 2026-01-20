# Function implementing cross-validation for linear regression

cv.predict.lm = function(df,
                         n.folds = NULL,
                         fold.ids, 
                         response.col,
                         predictor.cols,
                         data.selection,
                         feature.engineering.col,
                         feature.engineering.row,
                         feature.selection,
                         model,
                         seed = 12345){
  
  set.seed(seed)
  
  n <- nrow(df)
  p = ncol(df) - 2
  
  if (!is.null(n.folds)) {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  } else {
    n.folds = length(fold.ids)
  }
  
  pred.vec = rep(0, n)
  var.imp = matrix(0, nrow = p, ncol = n.folds)
  rownames(var.imp) = colnames(df %>% 
                                 select(-participant_id, 
                                        -all_of(response.col)))
  for (fold in 1:n.folds) {
    df.train = df[fold.ids != fold, ]
    df.test =   df[fold.ids == fold, ]
    predictor.cols.safe <- paste0("`", setdiff(colnames(df.train), c("participant_id", response.col)), "`")
    
    form <- as.formula(paste0(
      "`",
      response.col,
      "` ~ ",
      paste(predictor.cols.safe, collapse = " + ")
    ))
    
    # Fit the linear model
    fit <- lm(form, data = df.train)
    
    # Predict on test data
    pred.vec[fold.ids == fold] = predict(fit, newdata = df.test)
    
    # store coefficients in var.imp
    coefs <- coef(fit)
    
    # remove intercept if present
    coefs <- coefs[names(coefs) != "(Intercept)"]
    
    # store in var.imp: match by name
    var.names <- rownames(var.imp)
    for (v in names(coefs)) {
      # remove backticks to match rownames
      v.clean <- gsub("`", "", v)
      if (v.clean %in% var.names) {
        var.imp[v.clean, fold] <- coefs[v]
      }
    }
  }
  
  observed = df %>%
    select(all_of(response.col)) %>%
    as.matrix()
  
  predictions = data.frame("observed" = observed, "predicted" = pred.vec)
  
  colnames(predictions) = c("observed", "predicted")
  
  # Metrics
  R2 = cor(predictions$observed, predictions$predicted)^2
  R.spearman = cor(predictions$observed, predictions$predicted, method = "spearman")
  RMSE = rmse(predictions$observed, predictions$predicted)
  sRMSE = rmse(predictions$observed, predictions$predicted)/sd(predictions$observed)
  
  metrics = data.frame("data.selection" = data.selection,
                       "feature.engineering.col" = feature.engineering.col,
                       "feature.engineering.row" = feature.engineering.row,
                       "feature.selection" = feature.selection,
                       "model" = model,
                       "R2" = R2,
                       "R.spearman" = R.spearman,
                       "RMSE" = RMSE,
                       "sRMSE" = sRMSE)
  
  varImp <- data.frame(
    var         = rownames(var.imp),
    meanImp     = rowMeans(var.imp, na.rm = TRUE),
    sdImp       = apply(var.imp, 1, sd, na.rm = TRUE),
    varImp.type = "coefficient",
    stringsAsFactors = FALSE
  )
  
  # remove row names
  rownames(varImp) <- NULL
  
  prediction.plot = cv.plot(pred = predictions$predicted, obs = predictions$observed)
  
  results = list("predictions" = predictions, 
                 "metrics" = metrics,
                 "varImp" = varImp,
                 "prediction.plot" = prediction.plot)
  
  return(results)
  
}
