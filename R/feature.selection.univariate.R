# Function to perform univariate feature selection using cross-validated predictiveness
feature.selection.univariate = function(df,
                                        response.col,
                                        covariate.cols,
                                        model = "lm",
                                        metric = "sRMSE",
                                        metric.threshold = 1,
                                        n.folds = 5,
                                        fold.ids) {
  
  # Identify predictor columns: exclude participant_id, response, and covariates
  pred.cols <- df %>%
    select(-participant_id, -all_of(c(response.col, covariate.cols))) %>%
    colnames()
  
  n <- nrow(df)          # Number of observations
  p <- length(pred.cols) # Number of candidate predictors
  
  # Create fold IDs if not supplied
  if (!is.null(fold.ids)) {
    n.folds <- NULL  # Use supplied fold IDs
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # Initialize result dataframe to store metric for each predictor
  pred.res <- data.frame(
    pred = pred.cols,
    metric = 0,
    stringsAsFactors = FALSE
  )
  
  for (pred in pred.cols) {
    if (model == "lm") {
      # Initialize vector for storing predicted values
      pred.vec <- rep(NA_real_, nrow(df))
      
      # Loop over unique folds
      for (fold in unique(fold.ids)) {
        train.idx <- which(fold.ids != fold)
        test.idx  <- which(fold.ids == fold)
        
        # Prepare training data frame with correct column name
        train.df <- data.frame(y = df[[response.col]][train.idx], x = df[[pred]][train.idx])
        
        # Fit simple linear model: y ~ x
        lm.fit <- lm(y ~ x, data = train.df)
        
        # Prepare test data frame with same column name
        test.df <- data.frame(x = df[[pred]][test.idx])
        
        # Predict on test fold
        pred.vec[test.idx] <- predict(lm.fit, newdata = test.df)
      }
      
      # Compute the metric
      obs <- df[[response.col]]
      if (metric == "sRMSE") {
        RMSE <- sqrt(mean((obs - pred.vec)^2, na.rm = TRUE))
        pred.res$metric[pred.res$pred == pred] <- RMSE / sd(obs, na.rm = TRUE)
      } else if (metric == "RMSE") {
        pred.res$metric[pred.res$pred == pred] <- sqrt(mean((obs - pred.vec)^2, na.rm = TRUE))
      } else if (metric == "R2") {
        pred.res$metric[pred.res$pred == pred] <- cor(obs, pred.vec)^2
      } else if (metric == "R.spearman") {
        pred.res$metric[pred.res$pred == pred] <- cor(obs, pred.vec, method = "spearman")
      }
    }
  }
  
  
  # Select predictors based on threshold
  if (metric %in% c("sRMSE","RMSE")) {
    # Lower sRMSE is better
    pred.selected <- pred.res %>%
      filter(metric < metric.threshold) %>%
      pull(pred)
  } else {
    # For metrics where higher is better
    pred.selected <- pred.res %>%
      filter(metric > metric.threshold) %>%
      pull(pred)
  }
  
  # Return list with selected predictors and metrics for all predictors
  res.list <- list(
    pred.selected = pred.selected,
    pred.results = pred.res
  )
  
  return(res.list)
}
