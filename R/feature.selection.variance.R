# Function to perform feature selection based on predictor variance
feature.selection.variance = function(df,
                                      response.col,
                                      covariate.cols,
                                      predictor.cols,
                                      criterion = "threshold",
                                      metric.threshold = 0) {
  # Basic bookkeeping
  n <- nrow(df)
  p <- length(predictor.cols)
  
  # Prepare results data.frame
  pred.res <- data.frame(
    pred = predictor.cols,
    metric = NA_real_,            # store variance here (same column name as in univariate)
    stringsAsFactors = FALSE
  )
  
  # Compute variance for each predictor (NA's ignored)
  for (pred in predictor.cols) {
    vec <- df[[pred]]
    # Use var(..., na.rm = TRUE). If all NA or length < 2, var yields NA; keep that.
    pred.var <- tryCatch(var(vec, na.rm = TRUE), error = function(e) NA_real_)
    pred.res$metric[pred.res$pred == pred] <- pred.var
  }
  
  # Handle selection criteria
  if (criterion == "threshold") {
    # Select predictors with variance strictly greater than threshold
    pred.selected <- pred.res %>%
      dplyr::filter(metric > metric.threshold) %>%
      dplyr::pull(pred)
    
  } else if (criterion == "topN") {
    # threshold interpreted as integer N (top N by variance, most -> least)
    N <- as.integer(metric.threshold)
    if (is.na(N) || N <= 0) {
      pred.selected <- character(0)
    } else {
      pred.selected <- pred.res %>%
        dplyr::arrange(dplyr::desc(metric)) %>%
        dplyr::slice_head(n = N) %>%
        dplyr::pull(pred)
    }
    
  } else {
    stop("Unknown criterion. Use 'threshold' or 'topN'.")
  }
  
  # Return list (same structure as feature.selection.univariate)
  res.list <- list(pred.selected = pred.selected, pred.results = pred.res)
  return(res.list)
}