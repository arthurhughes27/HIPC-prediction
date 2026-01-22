# Script to derive cross-validated predictions for inactivated influenza vaccine

# Libraries
library(tidyverse)
library(fs)
library(Metrics)
library(ranger)
library(xgboost)

# Source internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

processed_data_folder <- "data"
output_folder = fs::path("output", "results")
df.predictor.list.path = fs::path(processed_data_folder,
                                  "engineered_dataframes_influenzain.rds")

df.clinical.path = fs::path(processed_data_folder, "hipc_merged_clinical_immresp.rds")

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
model = "ranger"
# Set the seed
seed = 21012026

# Generate the fold on the union of participant ids
ids = df.predictor.list[["d0"]][["none"]][["none"]] %>%
  pull(participant_id)
n = length(ids)
K <- 10
fold_id <- sample(rep(1:K, length.out = n))
fold_df <- tibble(participant_id = ids, fold = fold_id)

total.combinations = length(df.predictor.list) * length(df.predictor.list[["d0"]]) *
  length(df.predictor.list[["d0"]][["none"]]) + 1


res.list <- list()
i <- 1
for (data.sel in names(df.predictor.list)) {
  for (feat.eng.col in names(df.predictor.list[[data.sel]])) {
    for (feat.eng.row in names(df.predictor.list[[data.sel]][[feat.eng.col]])[-1]) {
      message(
        sprintf(
          "Running data.selection = %s | feat.eng.col = %s | feat.eng.row = %s | iteration = %d of %d" ,
          data.sel,
          feat.eng.col,
          feat.eng.row,
          i,
          total.combinations
        )
      )
      
      df.temp = df.predictor.list[[data.sel]][[feat.eng.col]][[feat.eng.row]]
      
      pids.temp = df.temp %>%
        pull(participant_id)
      
      fold.ids = fold_df %>%
        filter(participant_id %in% pids.temp) %>%
        pull(fold)
      
      res = cv.predict(
        df.predictor.list = df.predictor.list,
        df.clinical = df.clinical,
        covariate.cols = covariate.cols,
        response.col = response.col,
        data.selection = data.sel,
        feature.engineering.col = feat.eng.col,
        feature.engineering.row = feat.eng.row,
        feature.selection = feature.selection,
        model = model,
        fold.ids = fold.ids,
        seed = seed
      )
      
      res.list[[i]] <- res[["metrics"]]
      i <- i + 1
      
      gc()
    }
  }
}

metrics.df <- bind_rows(res.list)

# Derive baseline results
df = df.clinical %>%
  filter(participant_id %in% ids) %>%
  select(participant_id, all_of(covariate.cols), all_of(response.col)) %>%
  distinct()

pids.temp = df %>%
  pull(participant_id)

fold.ids = fold_df %>%
  filter(participant_id %in% pids.temp) %>%
  pull(fold)

baseline_results = cv.predict.baseline(
  df = df,
  predictor.cols = covariate.cols,
  response.col = response.col,
  model = "ranger",
  fold.ids = fold.ids,
  seed = seed,
  n.folds = NULL
)

metrics.df = bind_rows(metrics.df, baseline_results$metrics)

p_save <- fs::path(output_folder, "metrics_transformations_ranger.rds")

# Save the data
saveRDS(metrics.df, p_save)
