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
df.predictor.list.path = fs::path(processed_data_folder,
                                  "engineered_dataframes_influenzain.rds")

df.clinical.path = fs::path(processed_data_folder, "hipc_merged_clinical_immresp.rds")

# Load the data
df.predictor.list = readRDS(df.predictor.list.path)
df.clinical = readRDS(df.clinical.path)
covariate.cols = c("genderMale",
                   "age_imputed",
                   "immResp_MFC_hai_log2_pre_value")
response.col = "immResp_MFC_hai_log2_post_value"
data.selection = "d0"
feature.engineering.col = "none"
feature.engineering.row = "mean"
feature.selection = "none"
model = "xgboost"
seed = 12345
n = nrow(df.predictor.list[["d1"]][["none"]][["none"]])
n.folds = NULL
fold.ids = sample(rep(seq_len(10), length.out = n))

res = cv.predict(
  df.predictor.list = df.predictor.list,
  df.clinical = df.clinical,
  covariate.cols = covariate.cols,
  response.col = response.col,
  data.selection = data.selection,
  feature.engineering.col = feature.engineering.col,
  feature.engineering.row = feature.engineering.row,
  feature.selection = feature.selection,
  model = model,
  n.folds = 10,
  # fold.ids = fold.ids,
  seed = 20012026
)

