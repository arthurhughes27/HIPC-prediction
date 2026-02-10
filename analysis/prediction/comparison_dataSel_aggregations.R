# R script to perform prediction on a given study dataset and compare data selection and geneset aggregation approaches

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

# Study of interest
study_of_interest = "SDY1276"

# Dataset of interest 
dataset_of_interest = "all_noNorm"

# Assay of interest
assay_of_interest = "hai"

# Path for predictor sets
df.predictor.list.path = fs::path(
  processed_data_folder,
  paste0(
    "engineered_dataframes_influenzain_",
    dataset_of_interest,
    "_",
    study_of_interest,
    ".rds"
  )
)

# Path for clinical
df.clinical.path = fs::path(
  processed_data_folder,
  paste0("hipc_merged_clinical_immresp_", dataset_of_interest, ".rds")
)

# Load the data
df.predictor.list = readRDS(df.predictor.list.path)
df.clinical = readRDS(df.clinical.path)

# Define the covariates to always include
covariate.cols = c("genderMale", "age_imputed"
                   , paste0("immResp_mean_",assay_of_interest,"_pre_value")
                   )

# Define the response variable to predict
response.col = paste0("immResp_mean_",assay_of_interest,"_post_value")

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
total.combinations = length(names(df.predictor.list)) * length(names(df.predictor.list[["d0"]][["none"]]))

# Initialise results list
res.list <- list()
# Set counter
i <- 1

# Fix the model to random forest
mod = "elasticnet"

# Fix the feature engineering pre-transformation to none
feat.eng.col = "none"

# Fix the feature selection parameters
feature.selection = "none"
feature.selection.metric = "sRMSE"
feature.selection.metric.threshold = 1
feature.selection.model = "lm"
feature.selection.criterion = "relative.gain"

# Set the seed
seed = 10022026

# Fix the gender
gender.select = "Female"

for (data.sel in names(df.predictor.list)) {
  for (feat.eng.row in names(df.predictor.list[[data.sel]][[feat.eng.col]])) {
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
    
    if (!is.null(gender.select)) {
      if (gender.select == "Male") {
        df.clinical = df.clinical %>%
          filter(genderMale == 1)
      } else if (gender.select == "Female") {
        df.clinical = df.clinical %>%
          filter(genderMale == 0)
      }
    }
    
    # Extract the correct dataframe
    df.temp = df.predictor.list[[data.sel]][[feat.eng.col]][[feat.eng.row]]
    
    # Extract the relevant participant identifiers
    pids.expr = df.temp %>%
      pull(participant_id)
    
    pids.clinical = df.clinical %>%
      pull(participant_id)
    
    pids.temp = intersect(pids.expr, pids.clinical)
    
    
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
      feature.selection.metric = feature.selection.metric,
      feature.selection.metric.threshold = feature.selection.metric.threshold,
      feature.selection.model = feature.selection.model,
      feature.selection.criterion = feature.selection.criterion,
      model = mod,
      fold.ids = fold.ids,
      seed = seed,
      n.cores = 1,
      gender.select = gender.select
    )
    
    out.time = Sys.time()
    
    diff.time <- as.numeric(difftime(out.time, in.time, units = "mins"))
    
    message(sprintf("Completed in %.1f mins.", diff.time))
    
    
    # Store the metrics
    res.list[[i]] <- res[["metrics"]]
    # Iterate the counter
    i <- i + 1
    
    print(res[["metrics"]]$R2)
    
    # Clear the temporary memory
    gc()
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
pids.temp = df.clinical %>%
  pull(participant_id)

# Extract the relevant folds
fold.ids = fold_df %>%
  filter(participant_id %in% pids.temp) %>%
  pull(fold)

# Compute the cross-validation results
for (mod in c("lm")) {
  baseline_results = cv.predict.baseline(
    df = df,
    predictor.cols = covariate.cols,
    response.col = response.col,
    model = mod,
    fold.ids = fold.ids,
    seed = seed,
    n.folds = NULL,
    gender.select = gender.select
  )
  
  # Bind the baseline results to the metrics
  metrics.df = bind_rows(metrics.df, baseline_results$metrics)
}

metrics.df = metrics.df %>% 
  arrange(desc(R2))

# Path to save the results
file_name = paste0(
  "metrics_transformations_",
  mod,
  "_",
  dataset_of_interest,
  "_",
  study_of_interest,
  "_",
  assay_of_interest,
  "_",
  gender.select,
  ".rds"
)

p_save <- fs::path(output_folder, file_name)
# Save the results
saveRDS(metrics.df, p_save)
