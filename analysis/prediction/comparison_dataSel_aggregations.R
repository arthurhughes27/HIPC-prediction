# R script to perform prediction on a given study dataset and compare data selection and geneset aggregation approaches
# Libraries
library(tidyverse)
library(fs)
library(Metrics)
library(ranger)
library(xgboost)
library(foreach)
library(parallel)

seed = 04032026

# Source internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

# Data folder
processed_data_folder <- "data"

# Output folder to save results
output_folder = fs::path("output", "results")

# Figures folder to store graphics
figures_folder = fs::path("output", "figures", "diagnosis")

# Study of interest
study_of_interest = "SDY1276"

# Dataset of interest
dataset_of_interest = "all_noNorm"

# Assay of interest
assay_of_interest = "nAb"

# Gender of interest
gender_of_interest = "none"

# Model of interest
model_of_interest = "elasticnet"

# Vaccine of interest 
vaccine_of_interest = "Influenza (IN)"

# Response transformation of interest
response_transformation_of_interest = "mean"

# Response value of interest
response_value_of_interest = "post_value"

# Standardised response prior to aggregation?
response_standardised = TRUE

# Path for predictor sets
df.predictor.list.path = fs::path(
  processed_data_folder,
  paste0(
    "engineered_dataframes_",
    gsub("[[:space:]()]", "", tolower(vaccine_of_interest)),
    "_",
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
df.clinical = readRDS(df.clinical.path) %>% 
  mutate(ethnicityHisp = ifelse(ethnicity == "Hispanic or Latino", 1, 0),
         ethnicityOther = ifelse(ethnicity == "Other", 1, 0),)

# Define the response variable to predict
response.col = paste0("immResp_",
                      response_transformation_of_interest,
                      "_", 
                      assay_of_interest, 
                      ifelse(response_standardised, "_std_", "_log2_"),
                      response_value_of_interest)

response.col.pre = paste0("immResp_",
                          response_transformation_of_interest,
                          "_", 
                          assay_of_interest, 
                          ifelse(response_standardised, "_std_", "_log2_"),
                          "pre_value")

# Define the covariates to always include
covariate.cols = c(
  "genderMale",
  "age_imputed"#,
  # "ethnicityHisp",
  # "ethnicityOther",
  # response.col.pre
)

# Generate the fold on the union of participant ids
ids = df.predictor.list[["d0"]][["none"]][["none"]] %>%
  pull(participant_id)

# Number of participants
n = length(ids)

# Number of folds
K <- 10

# Get concerned participants
pid_df <- df.clinical %>%
  filter(participant_id %in% ids) %>%  # ensure same set
  distinct(participant_id, !!rlang::sym(response.col), .keep_all = TRUE)

# Generate folds with internal function balanced on covariates and response

fold_df = balance_folds(df = pid_df, 
                        ind.col = "participant_id",
                        covariate.cols = c(response.col),
                        n.folds = K,
                        n.continuous.split = 5)

fold.ids = fold_df$fold

# Calculate total number of iterations to do
total.combinations = length(names(df.predictor.list)) * length(names(df.predictor.list[["d0"]][["none"]])) * length(names(df.predictor.list[["d0"]]))

# Initialise results list
res.list <- list()
# Set counter
i <- 1

# Fix the feature engineering pre-transformation to none
# feat.eng.col = "none"

# Fix the feature selection parameters
feature.selection = "none"
feature.selection.metric = "sRMSE"
feature.selection.metric.threshold = 1
feature.selection.model = "lm"
feature.selection.criterion = "relative.gain"

# Set the seed
seed = 10022026

for (data.sel in names(df.predictor.list)) {
  for (feat.eng.col in names(df.predictor.list[[data.sel]])) {
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
          model_of_interest
        )
      )
      in.time = Sys.time()
      
      if (!is.null(gender_of_interest)) {
        if (gender_of_interest == "Male") {
          df.clinical = df.clinical %>%
            filter(genderMale == 1)
        } else if (gender_of_interest == "Female") {
          df.clinical = df.clinical %>%
            filter(genderMale == 0)
        }
      }
      
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
        model = model_of_interest,
        fold.ids = fold.ids,
        seed = seed,
        n.cores = 1,
        gender.select = gender_of_interest
      )
      
      out.time = Sys.time()
      
      diff.time <- as.numeric(difftime(out.time, in.time, units = "mins"))
      
      message(sprintf("Completed in %.1f mins.", diff.time))
      
      
      # Store the metrics
      res.list[[i]] <- res[["metrics"]]
      # Iterate the counter
      i <- i + 1
      
      print(res[["metrics"]]$R2)
      print(res$prediction.plot)
      
      # Clear the temporary memory
      gc()
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
    gender.select = gender_of_interest
  )
  
  # Bind the baseline results to the metrics
  metrics.df = bind_rows(metrics.df, baseline_results$metrics)
}

metrics.df = metrics.df %>%
  arrange(desc(R2))

# Path to save the results
file_name = paste0(
  "metrics_transformations_",
  gsub("[[:space:]()]", "", tolower(vaccine_of_interest)),
  "_",
  dataset_of_interest,
  "_",
  study_of_interest,
  "_",
  assay_of_interest,
  "_",
  gender_of_interest,
  "_",
  model_of_interest,
  ".rds"
)

p_save <- fs::path(output_folder, file_name)
# Save the results
saveRDS(metrics.df, p_save)
