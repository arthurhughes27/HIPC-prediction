# Script to derive cross-validated predictions for inactivated influenza vaccine

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
df.predictor.list.path = fs::path(processed_data_folder,
                                  "engineered_dataframes_influenzain_young_norm.rds")
# Path for clinical
df.clinical.path = fs::path(processed_data_folder, "hipc_merged_clinical_immresp_young_norm.rds")

# Load the data
df.predictor.list = readRDS(df.predictor.list.path)
df.clinical = readRDS(df.clinical.path)

# Define the covariates to always include
covariate.cols = c("genderMale",
                   "age_imputed",
                   "immResp_MFC_hai_log2_pre_value")

# Define the response variable to predict
response.col = "immResp_MFC_hai_log2_post_value"

# Fix the feature selection algorithm
feature.selection = "none"
# Fix the predictive model
model = "elasticnet"
# Set the seed
seed = 21012026

# Generate the fold on the union of participant ids
ids = df.predictor.list[["d0"]][["none"]][["none"]] %>%
  pull(participant_id)
# Number of participants
n = length(ids)
# Number of folds
K <- 10
# Randomly sampled fold ids
fold_id <- sample(rep(1:K, length.out = n))
# Attribute fold ids to individuals
fold_df <- tibble(participant_id = ids, fold = fold_id)

# Calculate total number of iterations to do
total.combinations = length(df.predictor.list) * length(df.predictor.list[["d0"]]) *
  length(df.predictor.list[["d0"]][["none"]]) * 3

# Initialise results list
res.list <- list()
# Set counter
i <- 1
for (mod in c("xgboost", "ranger", "elasticnet")) {
  for (data.sel in names(df.predictor.list)) {
    # For each combination of data selection and transformation
    for (feat.eng.col in names(df.predictor.list[[data.sel]])) {
      for (feat.eng.row in names(df.predictor.list[[data.sel]][[feat.eng.col]][-1])) {
        # Print a progress message
        message(
          sprintf(
            "Running data.selection = %s | feat.eng.col = %s | feat.eng.row = %s | iteration = %d of %d | model = %s",
            data.sel,
            feat.eng.col,
            feat.eng.row,
            i,
            total.combinations,
            mod
          )
        )
        
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
          feature.selection = feature.selection,
          model = mod,
          fold.ids = fold.ids,
          seed = seed,
          n.cores = 11
        )
        
        out.time = Sys.time()
        
        diff.time = out.time - in.time
        
        message(round(diff.time, 1), " mins elapsed.")
        
        
        # Store the metrics
        res.list[[i]] <- res[["metrics"]]
        # Iterate the counter
        i <- i + 1
        
        # Clear the temporary memory
        gc()
      }
    }
  }
}

# Bind the metrics into a dataframe
metrics.df <- bind_rows(res.list)

# Derive baseline results with just clinical information
df = df.clinical %>%
  filter(participant_id %in% ids) %>%
  select(participant_id, all_of(covariate.cols), all_of(response.col)) %>%
  distinct()

# Get the relevant participant ids
pids.temp = df %>%
  pull(participant_id)

# Extract the relevant folds
fold.ids = fold_df %>%
  filter(participant_id %in% pids.temp) %>%
  pull(fold)

# Compute the cross-validation results
for (mod in c("lm","xgboost", "ranger", "elasticnet")){
  baseline_results = cv.predict.baseline(
    df = df,
    predictor.cols = covariate.cols,
    response.col = response.col,
    model = mod,
    fold.ids = fold.ids,
    seed = seed,
    n.folds = NULL
  )
  
  # Bind the baseline results to the metrics
  metrics.df = bind_rows(metrics.df, baseline_results$metrics)
}

# Path to save the results
p_save <- fs::path(output_folder, "metrics_transformations_allmodels_young_norm.rds")
# Save the results
saveRDS(metrics.df, p_save)
