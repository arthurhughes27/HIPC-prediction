# R script to perform prediction on a given study dataset and compare data selection and geneset aggregation approaches
# Libraries
library(tidyverse)
library(fs)
library(Metrics)
library(ranger)
library(xgboost)
library(foreach)
library(parallel)

for (assay_of_interest in c("hai", "nAb")) {
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
  # assay_of_interest = "hai"
  
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
    paste0(
      "hipc_merged_clinical_immresp_",
      dataset_of_interest,
      ".rds"
    )
  )
  
  # Load the data
  df.predictor.list = readRDS(df.predictor.list.path)
  df.clinical = readRDS(df.clinical.path) %>%
    mutate(
      ethnicityHisp = ifelse(ethnicity == "Hispanic or Latino", 1, 0),
      ethnicityOther = ifelse(ethnicity == "Other", 1, 0),
    )
  
  # Define the response variable to predict
  response.col = paste0(
    "immResp_",
    response_transformation_of_interest,
    "_",
    assay_of_interest,
    ifelse(response_standardised, "_std_", "_log2_"),
    response_value_of_interest
  )
  
  response.col.pre = paste0(
    "immResp_",
    response_transformation_of_interest,
    "_",
    assay_of_interest,
    ifelse(response_standardised, "_std_", "_log2_"),
    "pre_value"
  )
  
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
  
  fold_df = balance_folds(
    df = pid_df,
    ind.col = "participant_id",
    covariate.cols = c(response.col),
    n.folds = K,
    n.continuous.split = 5,
    seed = seed
  )
  
  fold.ids = fold_df$fold
  
  # Derive baseline results with just clinical information
  df = df.clinical %>%
    filter(participant_id %in% ids) %>%
    select(participant_id,
           all_of(covariate.cols),
           all_of(response.col)) %>%
    distinct()
  
  # Compute the cross-validation results
  baseline_results = cv.predict.baseline(
    df = df,
    predictor.cols = covariate.cols,
    response.col = response.col,
    model = "lm",
    fold.ids = fold.ids,
    seed = seed,
    n.folds = NULL,
    gender.select = gender_of_interest
  )
  
  p1 <- baseline_results$prediction.plot
  
  plot_df <- p1$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p1$data <- plot_df
  
  p1 <- p1 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p1
  
  ggsave(
    filename = paste0("predictions_baseline_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p1,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now add gene expression, day 1, elastic net model
  
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d0",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p2 <- res$prediction.plot
  
  plot_df <- p2$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p2$data <- plot_df
  
  p2 <- p2 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p2
  
  ggsave(
    filename = paste0("predictions_baselineGE_d0_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p2,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now add gene expression, day 1, elastic net model
  
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d1",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p3 <- res$prediction.plot
  
  plot_df <- p3$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p3$data <- plot_df
  
  p3 <- p3 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p3
  
  ggsave(
    filename = paste0("predictions_baselineGE_d1_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p3,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  
  # Now add gene expression, day 3, elastic net model
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p4 = res$prediction.plot
  
  p4 <- res$prediction.plot
  
  plot_df <- p4$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p4$data <- plot_df
  
  p4 <- p4 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p4
  
  ggsave(
    filename = paste0("predictions_baselineGE_d3_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p4,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now add gene expression, day 0+1, elastic net model
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d0+d1",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p5 = res$prediction.plot
  
  p5 <- res$prediction.plot
  
  plot_df <- p5$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p5$data <- plot_df
  
  p5 <- p5 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p5
  
  ggsave(
    filename = paste0("predictions_baselineGE_d0d1_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p5,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now add gene expression, day 0+3, elastic net model
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d0+d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p6 = res$prediction.plot
  
  p6 <- res$prediction.plot
  
  plot_df <- p6$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p6$data <- plot_df
  
  p6 <- p6 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p6
  
  ggsave(
    filename = paste0("predictions_baselineGE_d0d3_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p6,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  
  
  # Now add gene expression, day 1+3, elastic net model
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d1+d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p7 = res$prediction.plot
  
  p7 <- res$prediction.plot
  
  plot_df <- p7$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p7$data <- plot_df
  
  p7 <- p7 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p7
  
  ggsave(
    filename = paste0("predictions_baselineGE_d1d3_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p7,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  
  
  
  # Now add gene expression, day 1+3, elastic net model
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d0+d1+d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest
  )
  
  p8 = res$prediction.plot
  
  p8 <- res$prediction.plot
  
  plot_df <- p8$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p8$data <- plot_df
  
  p8 <- p8 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p8
  
  ggsave(
    filename = paste0(
      "predictions_baselineGE_d0d1d3_",
      assay_of_interest,
      ".pdf"
    ),
    path = figures_folder,
    plot = p8,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  
  
  # Now d1 gene expression only, elastic net model
  
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d1",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest,
    include.covariates = F
  )
  
  p9 <- res$prediction.plot
  
  plot_df <- p9$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p9$data <- plot_df
  
  p9 <- p9 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p9
  
  ggsave(
    filename = paste0("predictions_GE_d1_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p9,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now d3 gene expression only, elastic net model
  
  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = gender_of_interest,
    include.covariates = F
  )
  
  p10 <- res$prediction.plot
  
  plot_df <- p10$data %>%
    mutate(participant_id = ind) %>%
    left_join(df.clinical %>%
                select(participant_id, gender) %>%
                distinct(),
              by = "participant_id")
  
  # Replace data inside ggplot object
  p10$data <- plot_df
  
  p10 <- p10 +
    aes(color = gender) +
    labs(color = "Gender")
  
  p10
  
  ggsave(
    filename = paste0("predictions_GE_d3_", assay_of_interest, ".pdf"),
    path = figures_folder,
    plot = p10,
    width = 22,
    height = 15,
    units = "cm"
  )
  
  # Now do baseline model for both males and females

  # Derive baseline results with just clinical information
  df = df.clinical %>%
    filter(participant_id %in% ids,
           genderMale == 1) %>%
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
  baseline_results = cv.predict.baseline(
    df = df,
    predictor.cols = covariate.cols,
    response.col = response.col,
    model = "lm",
    fold.ids = fold.ids,
    seed = seed,
    n.folds = NULL,
    gender.select =  "Male"
  )

  p11 = baseline_results$prediction.plot

  p11

  ggsave(
    filename = paste0("predictions_baselineMale_", assay_of_interest,".pdf"),
    path = figures_folder,
    plot = p11,
    width = 22,
    height = 15,
    units = "cm"
  )

  # Derive baseline results with just clinical information
  df = df.clinical %>%
    filter(participant_id %in% ids,
           genderMale == 0) %>%
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
  baseline_results = cv.predict.baseline(
    df = df,
    predictor.cols = covariate.cols,
    response.col = response.col,
    model = "lm",
    fold.ids = fold.ids,
    seed = seed,
    n.folds = NULL,
    gender.select =  "Female"
  )

  p12 = baseline_results$prediction.plot

  p12

  ggsave(
    filename = paste0("predictions_baselineFemale_", assay_of_interest,".pdf"),
    path = figures_folder,
    plot = p12,
    width = 22,
    height = 15,
    units = "cm"
  )

  # Males, d1 gene expression
  df.clinical.filtered = df.clinical %>%
    filter(genderMale == 1)

  # Extract the correct dataframe
  df.temp = df.predictor.list[["d1"]][["none"]][["mean"]]

  # Extract the relevant participant identifiers
  pids.expr = df.temp %>%
    pull(participant_id)

  pids.clinical = df.clinical.filtered %>%
    pull(participant_id)

  pids.temp = intersect(pids.expr, pids.clinical)

  # Extract the relevant folds
  fold.ids = fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)

  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical.filtered,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d1",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = "Male"
  )

  p13 = res$prediction.plot

  p13

  ggsave(filename = paste0("predictions_GE_d1_males_", assay_of_interest,".pdf"),
    path = figures_folder,
    plot = p13,
    width = 22,
    height = 15,
    units = "cm"
  )

  # Females, d1 gene expression
  df.clinical.filtered = df.clinical %>%
    filter(genderMale == 0)

  # Extract the correct dataframe
  df.temp = df.predictor.list[["d1"]][["none"]][["mean"]]

  # Extract the relevant participant identifiers
  pids.expr = df.temp %>%
    pull(participant_id)

  pids.clinical = df.clinical.filtered %>%
    pull(participant_id)

  pids.temp = intersect(pids.expr, pids.clinical)

  # Extract the relevant folds
  fold.ids = fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)

  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical.filtered,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d1",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = "Female"
  )

  p14 = res$prediction.plot

  p14

  ggsave(
    filename = paste0("predictions_GE_d1_females_", assay_of_interest,".pdf"),
    path = figures_folder,
    plot = p14,
    width = 22,
    height = 15,
    units = "cm"
  )



  # Males, d3 gene expression
  df.clinical.filtered = df.clinical %>%
    filter(genderMale == 1)

  # Extract the correct dataframe
  df.temp = df.predictor.list[["d3"]][["none"]][["mean"]]

  # Extract the relevant participant identifiers
  pids.expr = df.temp %>%
    pull(participant_id)

  pids.clinical = df.clinical.filtered %>%
    pull(participant_id)

  pids.temp = intersect(pids.expr, pids.clinical)

  # Extract the relevant folds
  fold.ids = fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)

  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical.filtered,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = "Male"
  )

  p15 = res$prediction.plot

  p15

  ggsave(filename = paste0("predictions_GE_d3_males_", assay_of_interest,".pdf"),
         path = figures_folder,
         plot = p15,
         width = 22,
         height = 15,
         units = "cm"
  )

  # Females, d3 gene expression
  df.clinical.filtered = df.clinical %>%
    filter(genderMale == 0)

  # Extract the correct dataframe
  df.temp = df.predictor.list[["d3"]][["none"]][["mean"]]

  # Extract the relevant participant identifiers
  pids.expr = df.temp %>%
    pull(participant_id)

  pids.clinical = df.clinical.filtered %>%
    pull(participant_id)

  pids.temp = intersect(pids.expr, pids.clinical)

  # Extract the relevant folds
  fold.ids = fold_df %>%
    filter(participant_id %in% pids.temp) %>%
    pull(fold)

  # Cross-validation
  res = cv.predict(
    df.predictor.list = df.predictor.list,
    df.clinical = df.clinical.filtered,
    covariate.cols = covariate.cols,
    response.col = response.col,
    data.selection = "d3",
    feature.engineering.col = "none",
    feature.engineering.row = "mean",
    feature.selection = "none",
    feature.selection.metric = "none",
    feature.selection.metric.threshold = "none",
    feature.selection.model = "none",
    feature.selection.criterion = "none",
    model = "elasticnet",
    fold.ids = fold.ids,
    seed = seed,
    n.cores = 1,
    gender.select = "Female"
  )

  p16 = res$prediction.plot

  p16

  ggsave(
    filename = paste0("predictions_GE_d3_females_", assay_of_interest,".pdf"),
    path = figures_folder,
    plot = p16,
    width = 22,
    height = 15,
    units = "cm"
  )
}

# Males, d1 gene expression
df.clinical.filtered = df.clinical %>%
  filter(genderMale == 1)

# Extract the correct dataframe
df.temp = df.predictor.list[["d1+d3"]][["none"]][["mean"]]

# Extract the relevant participant identifiers
pids.expr = df.temp %>%
  pull(participant_id)

pids.clinical = df.clinical.filtered %>%
  pull(participant_id)

pids.temp = intersect(pids.expr, pids.clinical)

# Extract the relevant folds
fold.ids = fold_df %>%
  filter(participant_id %in% pids.temp) %>%
  pull(fold)

# Cross-validation
res = cv.predict(
  df.predictor.list = df.predictor.list,
  df.clinical = df.clinical.filtered,
  covariate.cols = covariate.cols,
  response.col = response.col,
  data.selection = "d1+d3",
  feature.engineering.col = "none",
  feature.engineering.row = "none",
  feature.selection = "none",
  feature.selection.metric = "none",
  feature.selection.metric.threshold = "none",
  feature.selection.model = "none",
  feature.selection.criterion = "none",
  model = "elasticnet",
  fold.ids = fold.ids,
  seed = seed,
  n.cores = 1,
  gender.select = "Male"
)

p17 = res$prediction.plot

p17

ggsave(filename = paste0("predictions_GE_d1d3_males_", assay_of_interest,"_best.pdf"),
       path = figures_folder,
       plot = p17,
       width = 22,
       height = 15,
       units = "cm"
)




# Females, best
df.clinical.filtered = df.clinical %>%
  filter(genderMale == 0)

# Extract the correct dataframe
df.temp = df.predictor.list[["d0"]][["none"]][["mean"]]

# Extract the relevant participant identifiers
pids.expr = df.temp %>%
  pull(participant_id)

pids.clinical = df.clinical.filtered %>%
  pull(participant_id)

pids.temp = intersect(pids.expr, pids.clinical)

# Extract the relevant folds
fold.ids = fold_df %>%
  filter(participant_id %in% pids.temp) %>%
  pull(fold)

# Cross-validation
res = cv.predict(
  df.predictor.list = df.predictor.list,
  df.clinical = df.clinical.filtered,
  covariate.cols = covariate.cols,
  response.col = response.col,
  data.selection = "d0",
  feature.engineering.col = "none",
  feature.engineering.row = "pc1",
  feature.selection = "none",
  feature.selection.metric = "none",
  feature.selection.metric.threshold = "none",
  feature.selection.model = "none",
  feature.selection.criterion = "none",
  model = "elasticnet",
  fold.ids = fold.ids,
  seed = seed,
  n.cores = 1,
  gender.select = "Female"
)

p18 = res$prediction.plot

p18

ggsave(
  filename = paste0("predictions_GE_d0_females_", assay_of_interest,"_best.pdf"),
  path = figures_folder,
  plot = p18,
  width = 22,
  height = 15,
  units = "cm"
)

