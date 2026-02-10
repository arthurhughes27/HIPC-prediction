# File to cross check response pre-processing against theodora's work

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)
library(Metrics)  # for rmse function

tp = 1
tp_ab = 30

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_nAb <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical_all_norm.rds")

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# Clinical data for filtration of participants
hipc_clinical_all_norm = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = hipc_clinical_all_norm %>%
  filter(study_accession == "SDY1294",
         # gender == "Female",
         study_time_collected %in% c(0, tp)) %>%
  group_by(participant_id) %>% # keep only participants with Both 0 and 1
  filter(n_distinct(study_time_collected) == 2) %>%
  ungroup() %>%
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_nAb = response_nAb %>%
  filter(participant_id %in% participants,
         study_time_collected %in% c(0, tp_ab)) %>%
  dplyr::select(participant_id, study_time_collected, virus, value_preferred)

# First get the studies and other clinical data corresponding to each participant id
hipc_studies = hipc_clinical_all_norm %>%
  filter(participant_id %in% participants) %>%
  dplyr::select(
    participant_id,
    age_imputed,
    gender,
    study_accession,
    vaccine,
    vaccine_type,
    pathogen
  ) %>%
  distinct()

# Merge the study names into the immune response data (it is not directly given)
raw_response_influenzain = merge(x = response_nAb,
                                 y = hipc_studies,
                                 by = "participant_id",
                                 all = F) %>%
  dplyr::select(
    participant_id,
    study_accession,
    age_imputed,
    gender,
    virus,
    study_time_collected,
    value_preferred
  ) %>%
  filter(!is.na(value_preferred))

pre_vax <- raw_response_influenzain %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, virus) %>%
  slice_max(study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(response_MFC_pre_value = value_preferred,
         response_MFC_pre_time = study_time_collected) %>%
  dplyr::select(
    participant_id,
    study_accession,
    virus,
    response_MFC_pre_time,
    response_MFC_pre_value
  )

post_vax <- raw_response_influenzain %>%
  filter(study_time_collected > 0) %>%
  rename(response_MFC_post_value = value_preferred,
         response_MFC_post_time = study_time_collected) %>%
  dplyr::select(
    participant_id,
    study_accession,
    virus,
    response_MFC_post_time,
    response_MFC_post_value
  )

merged_vax <- full_join(
  pre_vax,
  post_vax,
  by = c("participant_id", "study_accession", "virus"),
  relationship = "many-to-many"
) %>%
  dplyr::select(
    participant_id,
    study_accession,
    virus,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value
  ) %>%
  filter(!is.na(response_MFC_pre_value),!is.na(response_MFC_post_value))

# Derive the log2 fold changes.
response_MFC_df <- merged_vax %>%
  mutate(
    # Define fold change
    response_MFC = response_MFC_post_value / response_MFC_pre_value,
    # If pre-vaccination value is 0, set fold change to post-vaccination value
    response_MFC = ifelse(
      is.infinite(response_MFC),
      response_MFC_post_value,
      response_MFC
    ),
    # Set fold change to 1 if both response_MFC_pre_value and response_MFC_post_value are 0
    response_MFC = ifelse(
      response_MFC_pre_value == 0 &
        response_MFC_post_value == 0,
      1,
      response_MFC
    )
  ) %>%
  # Calculate log2 fold change
  mutate(
    response_log2_MFC = ifelse(is.infinite(log2(response_MFC)), NA, log2(response_MFC)),
    response_log2_MFC_pre_value = log2(response_MFC_pre_value),
    response_log2_MFC_post_value = log2(response_MFC_post_value)
  ) %>%
  # Select relevant columns
  dplyr::select(
    participant_id,
    study_accession,
    virus,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value,
    response_MFC,
    response_log2_MFC_pre_value,
    response_log2_MFC_post_value,
    response_log2_MFC
  )

mean_response_nAb <- response_MFC_df %>%
  group_by(participant_id, study_accession) %>%
  summarise(
    immResp_mean_nAb_pre_value  = mean(response_MFC_pre_value, na.rm = TRUE),
    immResp_mean_nAb_post_value = mean(response_MFC_post_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    immResp_mean_nAb_log2_pre_value = log2(immResp_mean_nAb_pre_value),
    immResp_mean_nAb_log2_post_value = log2(immResp_mean_nAb_post_value),
    immResp_mean_nAb_FC = immResp_mean_nAb_post_value / immResp_mean_nAb_pre_value,
    immResp_mean_nAb_log2_FC = log2(immResp_mean_nAb_FC)
  )


clinical_df = hipc_studies %>% 
  dplyr::select(participant_id, 
         gender,
         age_imputed)

response_df = merge(x = clinical_df,
           y = mean_response_nAb,
           by = "participant_id") %>% 
  dplyr::select(participant_id, age_imputed, immResp_mean_nAb_pre_value, immResp_mean_nAb_post_value)

# File to cross-check expression data with theodora's work
library(tidyverse)
library(readxl)
# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data paths robustly
p_load_all_norm <- fs::path(raw_data_folder, "all_norm_eset.rds")
# Read in the rds file
all_norm_eset <- readRDS(p_load_all_norm)

# Load the expression data
all_norm_expr = all_norm_eset@assayData[["exprs"]] %>% 
  t() %>% 
  as.data.frame()

# Make the column names lowercase
colnames(all_norm_expr) = all_norm_expr %>% 
  colnames() %>% 
  tolower()

# Now extract information for the first two columns (participant_id and study_time_collected)
sample_info_all_norm = rownames(all_norm_expr)

# Use `stringr` and regular expressions to extract the participant_id and study_time_collected from the unique identifiers
matches_all_norm <- str_match(sample_info_all_norm , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")

# Participant ids
participant_id_all_norm <- matches_all_norm[, 2]

# Study times (numeric, rounded)
study_time_collected_all_norm <- matches_all_norm[, 3] %>% 
  as.numeric() %>% 
  round(2)

# Insert the identifying information as the first two columns 
all_norm_expr <- all_norm_expr %>%
  mutate(participant_id = participant_id_all_norm, study_time_collected = study_time_collected_all_norm) %>%
  dplyr::select(participant_id, study_time_collected, everything())

pids = response_df %>% 
  pull(participant_id) %>% 
  unique()

all_norm_expr_filtered = all_norm_expr %>% 
  filter(study_time_collected == tp,
         participant_id %in% pids)

all_merged = left_join(x = response_df,
                       y = all_norm_expr_filtered,
                       by = "participant_id") %>% 
  dplyr::select(-participant_id,
         -study_time_collected)

# theodora_genes_path <- fs::path("/", "home","ah3","Desktop","Work","PhD","PROJECT 3",
#                                 "Internship files TG","DAVIDGeneLists",
#                                 "sPLS_prediction_genes_d1_SDY1276.txt")
# 
# theodora_genes <- readr::read_csv(theodora_genes_path) %>% 
#   colnames() %>% 
#   tolower()

gene_names = all_norm_expr_filtered %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()

pred_df = all_merged %>% 
  mutate(response = log2(immResp_mean_nAb_post_value)) %>% 
  dplyr::select(age_imputed, 
         response,
         # immResp_mean_nAb_pre_value,
         any_of(gene_names))

pred_df  <- pred_df[, colSums(is.na(pred_df)) == 0]

set.seed(1)
library(caret)
library(ggplot2)

data <- na.omit(pred_df)  # drop rows with any NA

ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")

fit <- train(
  response ~ .,
  data = data,
  method = "glmnet",
  preProcess = c("center","scale"),
  trControl = ctrl,
  tuneLength = 10
)

# Extract CV predictions and observed
preds <- fit$pred$pred
obs   <- fit$pred$obs

# Standardized RMSE (divide by SD of observed)
std_rmse <- rmse(obs, preds) / sd(obs)

# R-squared
r2 <- cor(obs, preds)^2

# Plot with annotation
ggplot(fit$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed", y = "CV predicted") +
  theme_minimal() + 
  annotate(
    "text", x = 5, y = 8, 
    label = paste0("Std RMSE = ", round(std_rmse, 2), 
                   "\nR² = ", round(r2, 2)),
    hjust = 0
  )


