# Function to derive cross-validation predictions according to a given specification

cv.predict.baseline = function(df,
                               predictor.cols,
                               response.col,
                               model = "lm",
                               n.folds = NULL,
                               fold.ids = NULL,
                               seed = 12345,
                               gender.select = NULL) {
  set.seed(seed)
  
  if (!is.null(gender.select)) {
    if (gender.select == "Male") {
      df = df %>%
        filter(genderMale == 1)
    } else if (gender.select == "Female") {
      df = df %>%
        filter(genderMale == 0)
    }
  }
  
  ### CROSS-VALIDATION ###
  
  # Randomise folds
  
  n <- nrow(df)
  p = length(predictor.cols)
  
  if (!is.null(fold.ids)) {
    n.folds = NULL
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  
  if (model == "lm") {
    res = cv.predict.lm(
      df = df,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = "placeholder",
      feature.engineering.col = "placeholder",
      feature.engineering.row = "placeholder",
      feature.selection = "none",
      seed = seed,
      baseline = TRUE
    )
  } else if (model == "ranger") {
    res = cv.predict.ranger(
      df = df,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = "placeholder",
      feature.engineering.col = "placeholder",
      feature.engineering.row = "placeholder",
      feature.selection = "none",
      seed = seed,
      baseline = TRUE
    )
  } else if (model == "xgboost") {
    res = cv.predict.xgboost(
      df = df,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = "placeholder",
      feature.engineering.col = "placeholder",
      feature.engineering.row = "placeholder",
      feature.selection = "none",
      seed = seed,
      baseline = TRUE
    )
  } else if (model == "elasticnet") {
    res = cv.predict.elasticnet(
      df = df,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = "placeholder",
      feature.engineering.col = "placeholder",
      feature.engineering.row = "placeholder",
      feature.selection = "none",
      seed = seed,
      baseline = TRUE
    )
  }
  
  return(res)
}