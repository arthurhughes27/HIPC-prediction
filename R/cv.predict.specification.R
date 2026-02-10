cv.predict.specification = function(study_of_interest,
                                    dataset,
                                    data.sel,
                                    include.cov,
                                    covariate.cols,
                                    response.col,
                                    seed,
                                    K,
                                    mod,
                                    feat.eng.col,
                                    feat.eng.row,
                                    feat.select,
                                    feat.select.metric,
                                    feat.select.metric.threshold,
                                    feat.select.model,
                                    feat.select.criterion,
                                    n.cores,
                                    gender.select = NULL) {
  # Libraries
  library(tidyverse)
  library(fs)
  library(Metrics)
  library(ranger)
  library(xgboost)
  library(foreach)
  library(parallel)
  
  # Source internal functions
  sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)
  
  # Data folder
  processed_data_folder <- "data"
  # Output folder to save results
  output_folder = fs::path("output", "results")
  
  # Path for predictor sets
  df.predictor.list.path = fs::path(
    processed_data_folder,
    paste0(
      "engineered_dataframes_influenzain_",
      dataset,
      "_",
      study_of_interest,
      ".rds"
    )
  )
  
  # Path for clinical
  df.clinical.path = fs::path(processed_data_folder,
                              paste0("hipc_merged_clinical_immresp_", dataset, ".rds"))
  
  # Load the data
  df.predictor.list = readRDS(df.predictor.list.path)
  df.clinical = readRDS(df.clinical.path)
  
  if (!is.null(gender.select)) {
    if (gender.select == "Male") {
      df.clinical = df.clinical %>%
        filter(genderMale == 1)
    } else if (gender.select == "Female") {
      df.clinical = df.clinical %>%
        filter(genderMale == 0)
    }
  }
  
  # Generate the fold on the union of participant ids
  ids_exp = df.predictor.list[["d0"]][["none"]][["none"]] %>%
    pull(participant_id)
  
  ids_clinical = df.clinical %>% 
    pull(participant_id)
  
  ids = intersect(ids_exp, ids_clinical)
  
  # Number of participants
  n = length(ids)
  
  # Randomly sampled fold ids
  fold_id <- sample(rep(1:K, length.out = n))
  
  # Attribute fold ids to individuals
  fold_df <- tibble(participant_id = ids, fold = fold_id)
  
  in.time = Sys.time()
  
  # Extract the correct dataframe
  df.temp = df.predictor.list[[data.sel]][[feat.eng.col]][[feat.eng.row]]
  
  # Extract the relevant participant identifiers
  pids.temp = df.temp %>%
    pull(participant_id)
  
  # Extract the relevant folds
  fold.ids = fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)
  
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = data.sel,
    feature.engineering.col = feat.eng.col,
    feature.engineering.row = feat.eng.row,
    feature.selection = feat.select,
    feature.selection.metric = feat.select.metric,
    feature.selection.metric.threshold = feat.select.metric.threshold,
    feature.selection.model = feat.select.model,
    feature.selection.criterion = feat.select.criterion,
    include.covariates = include.cov,
    model = mod,
    fold.ids = fold.ids,
    seed = seed,
    n.cores = n.cores,
    gender.select = gender.select
  )
  
  out.time = Sys.time()
  
  diff.time <- as.numeric(difftime(out.time, in.time, units = "mins"))
  
  print(message(sprintf("Completed in %.1f mins.", diff.time)))
  
  print(res$metrics$R2)
  print(res$prediction.plot)
  return(res)
}