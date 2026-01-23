# Function to derive cross-validation predictions according to a given specification

cv.predict = function(df.predictor.list,
                      df.clinical,
                      covariate.cols,
                      response.col,
                      data.selection = "d0",
                      feature.engineering.col = "none",
                      feature.engineering.row = "mean",
                      feature.selection = "none",
                      model = "lm",
                      n.folds = NULL,
                      fold.ids = NULL,
                      seed = 12345,
                      n.cores = 1) {
  set.seed(seed)
  
  ### DATA SELECTION & FEATURE ENGINEERING ###
  
  # Select data from the engineered list according to the user parameter
  if (data.selection %in% names(df.predictor.list)) {
    df.predictors = df.predictor.list[[data.selection]][[feature.engineering.col]][[feature.engineering.row]] %>%
      select(-study_time_collected)
    
    # Select relevant columns from the clinical data
    df.clinical = df.clinical %>%
      select(participant_id,
             all_of(covariate.cols),
             all_of(response.col))
    
    # Merge these into one dataframe
    df.all = right_join(x = df.clinical, y = df.predictors, by = "participant_id") %>%
      distinct()
    
  } else {
    pids_to_retain = df.predictor.list[[1]][[1]][[1]] %>%
      pull(participant_id) %>%
      unique()
    
    df.all <- df.clinical %>%
      filter(participant_id %in% pids_to_retain) %>%
      select(participant_id,
             all_of(covariate.cols),
             all_of(response.col)) %>%
      distinct() %>%
      drop_na()
    
  }
  
  
  ### FEATURE SELECTION
  
  ### CROSS-VALIDATION ###
  
  # Randomise folds
  
  n <- nrow(df.all)
  p = ncol(df.all %>%
             select(-participant_id))
  
  if (!is.null(fold.ids)){
    n.folds = NULL
  } else {
    fold.ids <- sample(rep(seq_len(n.folds), length.out = n))
  }
  
  
  if (model == "lm") {
    res = cv.predict.lm(
      df = df.all,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = data.selection,
      feature.engineering.col = feature.engineering.col,
      feature.engineering.row = feature.engineering.row,
      feature.selection = feature.selection,
      seed = seed
    )
  } else if (model == "ranger") {
    res = cv.predict.ranger(
      df = df.all,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = data.selection,
      feature.engineering.col = feature.engineering.col,
      feature.engineering.row = feature.engineering.row,
      feature.selection = feature.selection,
      seed = seed,
      n.cores = n.cores
    )
  } else if (model == "xgboost") {
    res = cv.predict.xgboost(
      df = df.all,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = data.selection,
      feature.engineering.col = feature.engineering.col,
      feature.engineering.row = feature.engineering.row,
      feature.selection = feature.selection,
      seed = seed
    )
  } else if (model == "elasticnet"){
    res = cv.predict.elasticnet(
      df = df.all,
      fold.ids = fold.ids,
      response.col = response.col,
      data.selection = data.selection,
      feature.engineering.col = feature.engineering.col,
      feature.engineering.row = feature.engineering.row,
      feature.selection = feature.selection,
      seed = seed
    )
  }
  
  return(res)
}