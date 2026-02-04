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
                                  "engineered_dataframes_influenzain_all_norm.rds")
# Path for clinical
df.clinical.path = fs::path(processed_data_folder,
                            "hipc_merged_clinical_immresp_all_norm.rds")

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
feature.selection = "univariate"
# Fix the predictive model
model = "lm"
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

# Initialise results list
res.list <- list()
# Set counter
i <- 1

feat.select = "univariate"
mod = "elasticnet"
data.sel = "d1"
feat.eng.col = "none"
feat.eng.row = "mean"

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

data.selection = data.sel
feature.engineering.col = feat.eng.col
feature.engineering.row = feat.eng.row
feature.selection = feat.select
feature.selection.metric = "sRMSE"
feature.selection.criterion = "relative.gain"
feature.selection.metric.threshold = 0.75
feature.selection.model = "lm"
model = mod
n.cores = 12
include.covariates = TRUE

# Cross-validation
res = cv.predict(
  df.predictor.list = df.predictor.list,
  df.clinical = df.clinical,
  covariate.cols = covariate.cols,
  response.col = response.col,
  data.selection = data.sel,
  feature.engineering.col = feat.eng.col,
  feature.engineering.row = feat.eng.row,
  feature.selection = "none",
  feature.selection.metric = "sRMSE",
  feature.selection.metric.threshold = 1,
  feature.selection.model = "lm",
  feature.selection.criterion = "relative.gain",
  include.covariates = F,
  model = mod,
  fold.ids = fold.ids,
  seed = seed,
  n.cores = 1
)

out.time = Sys.time()

diff.time <- as.numeric(difftime(out.time, in.time, units = "mins"))

message(sprintf("Completed in %.1f mins.", diff.time))



