# File to cross check response pre-processing against theodora's work

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)


# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_nAb <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical_all_noNorm.rds")

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# Clinical data for filtration of participants
hipc_clinical_all_noNorm = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = hipc_clinical_all_noNorm %>%
  filter(study_accession == "SDY1276",
         gender == "Female",
         study_time_collected %in% c(0, 1)) %>%
  group_by(participant_id) %>% # keep only participants with Both 0 and 1
  filter(n_distinct(study_time_collected) == 2) %>%
  ungroup() %>%
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_nAb = response_nAb %>%
  filter(participant_id %in% participants,
         study_time_collected %in% c(0, 28)) %>%
  select(participant_id, study_time_collected, virus, value_preferred)

# First get the studies and other clinical data corresponding to each participant id
hipc_studies = hipc_clinical_all_noNorm %>%
  filter(participant_id %in% participants) %>%
  select(
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
  select(
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
  select(
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
  select(
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
  select(
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
  select(
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

check_df1 = mean_response_nAb %>%
  mutate(response = immResp_mean_nAb_post_value) %>% 
  select(participant_id,
         response) %>%
  arrange(participant_id)

check_df3 = mean_response_nAb %>%
  mutate(response = immResp_mean_nAb_pre_value) %>% 
  select(participant_id,
         response) %>%
  arrange(participant_id)

# Loading the data
## Define the paths
raw_data_path = fs::path("data-raw")
path_ge = fs::path(raw_data_path, "all_norm_withResponse_eset.rds")
path_nab = fs::path(raw_data_path, "neut_ab_titer_2025-01-10_01-13-22.xlsx")

# Load the data
data = readRDS(path_ge)
neut.ab.assay <- read_excel(path_nab)

# Extract the clinical data
rawdata <- as.data.frame(data@phenoData@data)
# Extract the gene expression data
ge.data <- as.data.frame(data@assayData[["exprs"]])
rm(data)

# Preprocessing
## Select only females with 2008 TIV vaccine with gene expression info at baseline and 1 day, from SDY1276
rawdata <- rawdata %>%
  filter(
    pathogen == "Influenza",
    vaccine == "TIV (2008)",
    gender == "Female",
    study_time_collected %in% c(0, 1),
    study_accession == "SDY1276"
  ) %>%
  rename(expression_time = study_time_collected) %>%
  group_by(participant_id) %>% # keep only participants with Both 0 and 1
  filter(n_distinct(expression_time) == 2) %>%
  ungroup()

# Extract participant ids from filtered data
filtered_ids = rawdata %>%
  pull(participant_id) %>%
  unique()

# Filter particpants from neutralizing antibody dataframe
## Only choose participants with antibody levels at day 0 and 28 without missing values
nab_filtered <- neut.ab.assay %>%
  filter(
    `Participant ID` %in% filtered_ids,
    `Study Time Collected` %in% c(0, 28),
    !is.na(`Value Preferred`)
  )

# Grab the average of the antibody response
nab_composite <- nab_filtered %>%
  group_by(`Participant ID`, `Study Time Collected`) %>%
  summarize(response = mean(`Value Preferred`, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(`Participant ID`) %>%
  filter(n_distinct(`Study Time Collected`) == 2) %>%
  ungroup()

check_df2 = nab_composite %>% 
  filter(`Study Time Collected` == 28) %>% 
  arrange(`Participant ID`) %>% 
  select(`Participant ID`,
         response)

check_df4 = nab_composite %>% 
  filter(`Study Time Collected` == 0) %>% 
  arrange(`Participant ID`) %>% 
  select(`Participant ID`,
         response)

# Check a basic predictive lm model : age, pre-vaccine antibodies on post-vaccine antibodies

clinical_df = hipc_studies %>% 
  select(participant_id, 
         gender,
         age_imputed)

df = merge(x = clinical_df,
           y = mean_response_nAb,
           by = "participant_id")

set.seed(1)
library(caret)
library(ggplot2)

ctrl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final"
)

fit <- train(
  immResp_mean_nAb_post_value ~ age_imputed + immResp_mean_nAb_pre_value,
  data = df,
  method = "lm",
  trControl = ctrl
)

# CV predictions vs observed
ggplot(fit$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "CV predicted", y = "Observed") +
  theme_minimal() + 
  coord_fixed() + 
  xlim(4,12) + 
  ylim(4,12)

