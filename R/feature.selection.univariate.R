# Function to perform univariate feature selection using cross-validated predictiveness
feature.selection.univariate = function(df,
                                        response.col,
                                        covariate.cols,
                                        predictor.cols,
                                        model = "lm",
                                        metric = "sRMSE",
                                        criterion = "relative.gain",
                                        metric.threshold = 1,
                                        include.covariates = TRUE,
                                        n.folds = 5,
                                        fold.ids = NULL) {
  n <- nrow(df)          # Number of observations
  p <- length(predictor.cols) # Number of candidate predictors
  
  # Create fold IDs if not supplied
  if (!is.null(fold.ids)) {
    n.folds <- NULL  # Use supplied fold IDs
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  # Initialize result dataframe to store metric for each predictor
  pred.res <- data.frame(pred = predictor.cols,
                         metric = 0,
                         stringsAsFactors = FALSE)
  
  if (criterion == "relative.gain") {
    pred.vec <- rep(NA_real_, nrow(df))
    for (fold in sort(unique(fold.ids))) {
      train.idx <- which(fold.ids != fold)
      test.idx  <- which(fold.ids == fold)
      
      if (model == "lm") {
        if (include.covariates) {
          train.df <- df[train.idx, c(response.col, covariate.cols), drop = FALSE]
          test.df  <- df[test.idx, covariate.cols, drop = FALSE]
          
          lm.fit <- lm(reformulate(covariate.cols, response = response.col),
                       data = train.df)
          
          pred.vec[test.idx] <- predict(lm.fit, newdata = test.df)
          
        } else {
          train.df <- df[train.idx, response.col, drop = FALSE]
          
          # Intercept-only model
          lm.fit <- lm(reformulate(character(0), response = response.col),
                       data = train.df)
          
          # No predictors needed
          pred.vec[test.idx] <- coef(lm.fit)[1]
        }
      }
      
      if (model == "ranger") {
        if (include.covariates) {
          train.df <- df[train.idx, c(response.col, covariate.cols), drop = FALSE]
          test.df  <- df[test.idx, covariate.cols, drop = FALSE]
          
          rf.fit <- ranger::ranger(
            formula = reformulate(covariate.cols, response = response.col),
            data = train.df
          )
          
          pred.vec[test.idx] <- predict(rf.fit, data = test.df)$predictions
          
        } else {
          train.df <- df[train.idx, response.col, drop = FALSE]
          pred.vec[test.idx] <- mean(train.df[[response.col]], na.rm = TRUE)
        }
      }
    }
    
    obs <- df[[response.col]]
    if (metric == "sRMSE") {
      RMSE <- sqrt(mean((obs - pred.vec)^2, na.rm = TRUE))
      pred.res.baseline <- RMSE / sd(obs, na.rm = TRUE)
    } else if (metric == "RMSE") {
      pred.res.baseline <- sqrt(mean((obs - pred.vec)^2, na.rm = TRUE))
    } else if (metric == "R2") {
      pred.res.baseline <- cor(obs, pred.vec)^2
    } else if (metric == "R.spearman") {
      pred.res.baseline <- cor(obs, pred.vec, method = "spearman")
    }
  }
  
  for (pred in predictor.cols) {
    pred.vec <- rep(NA_real_, nrow(df))
    
    for (fold in unique(fold.ids)) {
      train.idx <- which(fold.ids != fold)
      test.idx  <- which(fold.ids == fold)
      
      if (model == "lm") {
        if (include.covariates) {
          train.df <- df[train.idx, c(response.col, pred, covariate.cols), drop = FALSE]
          test.df  <- df[test.idx, c(pred, covariate.cols), drop = FALSE]
          metric
          rhs.cols <- c(pred, covariate.cols)
          rhs.cols <- paste0("`", rhs.cols, "`")  # backticks for formula
          
          lm.fit <- lm(reformulate(rhs.cols, response = response.col),
                       data = train.df)
          pred.vec[test.idx] <- predict(lm.fit, newdata = test.df)
          
        } else {
          train.df <- df[train.idx, c(response.col, pred), drop = FALSE]
          test.df  <- df[test.idx, c(pred), drop = FALSE]
          
          rhs.cols <- c(pred)
          rhs.cols <- paste0("`", rhs.cols, "`")  # backticks for formula
          
          lm.fit <- lm(reformulate(rhs.cols, response = response.col),
                       data = train.df)
          pred.vec[test.idx] <- predict(lm.fit, newdata = test.df)
        }
      } else if (model == "ranger") {
        if (include.covariates) {
          train.df <- df[train.idx, c(response.col, pred, covariate.cols), drop = FALSE]
          test.df  <- df[test.idx, c(pred, covariate.cols), drop = FALSE]
          
          rf.fit <- ranger::ranger(dependent.variable.name = response.col,
                                   data = train.df)
          
          pred.vec[test.idx] <- predict(rf.fit, data = test.df)$predictions
          
        } else {
          train.df <- df[train.idx, c(response.col, pred), drop = FALSE]
          test.df  <- df[test.idx, c(pred), drop = FALSE]
          
          rf.fit <- ranger::ranger(dependent.variable.name = response.col,
                                   data = train.df)
          
          pred.vec[test.idx] <- predict(rf.fit, data = test.df)$predictions
        }
      }
    }
    
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
  
  
  if (criterion == "threshold") {
    pred.res = pred.res %>% 
      mutate(baseline.metric = pred.res.baseline,
             selected = ifelse(pred %in% pred.selected, TRUE, FALSE)) %>% 
      dplyr::select(pred, baseline.metric, metric, selected)
    
    # Select predictors based on threshold
    if (metric %in% c("sRMSE", "RMSE")) {
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
  } else if (criterion %in% c("relative.gain", "relative.gain.absolute")) {
    
    if (criterion == "relative.gain") {
      if (metric %in% c("sRMSE", "RMSE")) {
        pred.res = pred.res %>%
          mutate(relative.gain = 100 * (pred.res.baseline - metric) / pred.res.baseline)
      } else {
        pred.res = pred.res %>%
          mutate(relative.gain = -100 * (pred.res.baseline - metric) / pred.res.baseline)
      }
    } else if (criterion == "relative.gain.absolute") {
      if (metric %in% c("sRMSE", "RMSE")) {
        pred.res = pred.res %>%
          mutate(relative.gain = pred.res.baseline - metric)
      } else {
        pred.res = pred.res %>%
          mutate(relative.gain = pred.res.baseline - metric)
      }
    }
    
    pred.selected <- pred.res %>%
      filter(relative.gain > metric.threshold) %>%
      pull(pred)
    
    pred.res = pred.res %>%
      mutate(metric.original = metric,
             baseline.metric.original = pred.res.baseline,
             selected = ifelse(pred %in% pred.selected, TRUE, FALSE)) %>%
      mutate(metric = relative.gain) %>%
      dplyr::select(pred, baseline.metric.original, metric.original, metric, selected)
  }
  
  
  # Return list with selected predictors and metrics for all predictors
  res.list <- list(pred.selected = pred.selected, pred.results = pred.res)
  
  return(res.list)
}
