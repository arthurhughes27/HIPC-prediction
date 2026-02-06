# Script to investigate univariate predictiveness measures for feature selection

library(dplyr)

# Source internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

# Data folder
processed_data_folder <- "data"
# Output folder to save results
output_folder = fs::path("output", "results")

study_of_interest = "SDY80"
# Path for predictor sets
df.predictor.list.path = fs::path(processed_data_folder,
                                  paste0("engineered_dataframes_influenzain_all_noNorm_",study_of_interest,".rds"))
# Path for clinical
df.clinical.path = fs::path(processed_data_folder,
                            "hipc_merged_clinical_immresp_all_noNorm.rds")

# Load the data
df.predictor.list = readRDS(df.predictor.list.path)
df.clinical = readRDS(df.clinical.path)

# Define the covariates to always include
covariate.cols = c("genderMale",
                   "age_imputed",
                   "immResp_MFC_nAb_log2_pre_value")

# Define the response variable to predict
response.col = "immResp_MFC_nAb_log2_post_value"

# Number of folds
K <- 10

# Generate the fold on the union of participant ids
ids = df.predictor.list[["d0"]][["none"]][["none"]] %>%
  pull(participant_id)

# Number of participants
n = length(ids)

# Randomly sampled fold ids
fold_id <- sample(rep(1:K, length.out = n))

# Attribute fold ids to individuals
fold_df <- tibble(participant_id = ids, fold = fold_id)

# Set the seed
seed = 21012026

# Initialise results list
res.list <- list()
feat.eng.col = "none"
feat.select = "univariate"
data.sel = "d3"
feat.eng.row = "mean"

predictor.cols = df.predictor.list[["d0"]][["none"]][["mean"]] %>%
  select(-any_of(
    c(
      "participant_id",
      "study_time_collected",
      response.col,
      covariate.cols
    )
  )) %>%
  colnames()

feature.selection.results.df = data.frame(
  "predictor" = predictor.cols,
  "data.selection" = NA_real_,
  "feature.engineering.col" = NA_real_,
  "feature.engineering.row" = NA_real_,
  "metric" = NA_real_
)


# collect every run's pred-level results here
results_list <- list()
ctr <- 1

for (data.sel in names(df.predictor.list)) {
  in.time <- Sys.time()
  
  # Extract the correct dataframe
  df.temp <- df.predictor.list[[data.sel]][[feat.eng.col]][[feat.eng.row]]
  
  df.clinical <- df.clinical %>%
    select(participant_id,
           all_of(covariate.cols),
           all_of(response.col))
  
  df.all <- right_join(x = df.clinical, y = df.temp, by = "participant_id") %>%
    distinct()
  
  pids.temp <- df.temp %>% pull(participant_id)
  
  fold.ids <- fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)
  
  feature.selection.results <- feature.selection.univariate(
    df = df.all,
    response.col = response.col,
    covariate.cols = covariate.cols,
    predictor.cols = predictor.cols,
    model = "lm",
    metric = "R2",
    criterion = "relative.gain",
    metric.threshold = 2,
    include.covariates = TRUE,
    fold.ids = fold.ids
  )
  
  out.time <- Sys.time()
  message(sprintf("Completed in %.1f mins.", as.numeric(
    difftime(out.time, in.time, units = "mins")
  )))
  
  # build the predictor × specification rows for this run
  pred_res <- feature.selection.results[["pred.results"]] %>%
    select(pred, metric) %>%
    rename(predictor = pred) %>%
    mutate(
      data.selection = data.sel,
      feature.engineering.row = feat.eng.row,
      feature.engineering.col = feat.eng.col
    )
  
  # store
  results_list[[ctr]] <- pred_res
  ctr <- ctr + 1
}

# final stacked dataframe: one row per predictor × specification
feature.selection.results.all.df <- bind_rows(results_list) %>%
  select(
    predictor,
    data.selection,
    feature.engineering.col,
    feature.engineering.row,
    metric
  )

# Path to save the results
p_save <- fs::path(output_folder, "feature.selection.lm.mean.all.rds")
# Save the results
saveRDS(feature.selection.results.all.df, p_save)
