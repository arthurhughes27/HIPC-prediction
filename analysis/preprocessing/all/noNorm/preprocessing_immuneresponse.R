# File for preprocessing the raw immune response data files 
# We do this for all flu studies 
# Output is a dataframe in long format giving the participant, study, strain, assay (HAI or nAb),
# timepoint and value

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
p_load_hai <- fs::path(raw_data_folder, "hai_2025-01-10_01-13-41.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical_all_noNorm.rds")

# HAI response data
response_hai = read_excel(p_load_hai) %>%
  clean_names()

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# Clinical data for filtration of participants
hipc_clinical_all_noNorm = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = hipc_clinical_all_noNorm %>%
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_hai = response_hai %>%
  filter(participant_id %in% participants)

response_nAb = response_nAb %>%
  filter(participant_id %in% participants)

# Depending on the assay, there may be multiple viral strains or analytes measured.
# We rename some columns and add a column for the assay name, so that we can merge these data together.

response_hai = response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai")

response_nAb = response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb")

# Now merge all the raw response data
raw_response_influenzain = bind_rows(response_nAb, response_hai) %>%
  select(-gender) %>% 
  arrange(participant_id) %>%
  distinct()

# First get the studies and other clinical data corresponding to each participant id
hipc_studies = hipc_clinical_all_noNorm %>%
  select(participant_id,
         age_imputed,
         gender,
         study_accession,
         vaccine,
         vaccine_type,
         pathogen) %>%
  distinct()

# Merge the study names into the immune response data (it is not directly given)
raw_response_influenzain = merge(x = raw_response_influenzain,
                                    y = hipc_studies,
                                    by = "participant_id",
                                    all = F) %>%
  select(
    participant_id,
    study_accession,
    age_imputed, 
    gender,
    response_strain_analyte,
    assay,
    study_time_collected,
    value_preferred
  )


# There are a series of study-specific errors which have been inferred from the "immuneResponseCallGeneration.R" script in the ImmuneSignatures2 GitHub page.
# We fix these here.
# First, "SDY1276" is apparently scaled in log4.

raw_response_influenzain = raw_response_influenzain %>%
  mutate(value_preferred = if_else(
    study_accession == "SDY1276",
    4^value_preferred,
    value_preferred
  ))

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "raw_response_influenzain_all_noNorm.rds")

# Save dataframe
saveRDS(raw_response_influenzain, file = p_save)

# Now derive maximum fold-change values for each individual

# We intend to predict on immune response data at day 28 (plus or minus 7 days), so we filter to get the immune data at these timepoints or pre-vaccination.

raw_response_influenzain = raw_response_influenzain %>%
  filter((study_time_collected < 36 &
            study_time_collected > 20) | study_time_collected <= 0)


# Now derive MFC values, taking most recent measurement as baseline in the case where there are two pre-vaccination measurements

pre_vax <- raw_response_influenzain %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, assay, response_strain_analyte) %>%
  slice_max(study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(response_MFC_pre_value = value_preferred,
         response_MFC_pre_time = study_time_collected) %>% 
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_pre_time, response_MFC_pre_value)

post_vax <- raw_response_influenzain %>%
  filter(study_time_collected > 0) %>%
  rename(response_MFC_post_value = value_preferred,
         response_MFC_post_time = study_time_collected) %>% 
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_post_time, response_MFC_post_value)

# Now we merge pre- and post- values to get a dataframe in wide format with one line per unique combination of participant/assay/strain-analyte
merged_vax <- full_join(
  pre_vax,
  post_vax,
  by = c("participant_id", "study_accession","response_strain_analyte", "assay"),
  relationship = "many-to-many"
) %>%
  # mutate(genderMale = genderMale.x,
  #        age_imputed = age_imputed.x,
  #        study_accession = study_accession.x) %>%
  select(
    participant_id,
    study_accession,
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value
  )

# Now we can exclude any rows with NA values in either pre- or post-.
merged_vax = merged_vax %>%
  filter(!is.na(response_MFC_pre_value),
         !is.na(response_MFC_post_value))

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
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value,
    response_MFC, 
    response_log2_MFC_pre_value,
    response_log2_MFC_post_value,
    response_log2_MFC
  )



# Now, for each participant, we can derive the maximum of these fold changes.
## Across all assays
max_response_MFC_df_any <- response_MFC_df %>%
  group_by(participant_id) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup() %>%
  rename(
    immResp_MFC_anyAssay_response_strain_analyte = response_strain_analyte,
    immResp_MFC_anyAssay_assay = assay,
    immResp_MFC_anyAssay_pre_time = response_MFC_pre_time,
    immResp_MFC_anyAssay_post_time = response_MFC_post_time,
    immResp_MFC_anyAssay_pre_value = response_MFC_pre_value,
    immResp_MFC_anyAssay_post_value = response_MFC_post_value,
    immResp_MFC_anyAssay_MFC = response_MFC,
    immResp_MFC_anyAssay_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_anyAssay_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_anyAssay_log2_MFC = response_log2_MFC
  )

## And within each assay
max_response_MFC_df_each <- response_MFC_df %>%
  group_by(participant_id, assay) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup()

max_response_MFC_df_nAb = max_response_MFC_df_each %>%
  filter(assay == "nAb") %>%
  rename(
    immResp_MFC_nAb_response_strain_analyte = response_strain_analyte,
    immResp_MFC_nAb_pre_time = response_MFC_pre_time,
    immResp_MFC_nAb_post_time = response_MFC_post_time,
    immResp_MFC_nAb_pre_value = response_MFC_pre_value,
    immResp_MFC_nAb_post_value = response_MFC_post_value,
    immResp_MFC_nAb_MFC = response_MFC,
    immResp_MFC_nAb_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_nAb_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_nAb_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

max_response_MFC_df_hai = max_response_MFC_df_each %>%
  filter(assay == "hai") %>%
  rename(
    immResp_MFC_hai_response_strain_analyte = response_strain_analyte,
    immResp_MFC_hai_pre_time = response_MFC_pre_time,
    immResp_MFC_hai_post_time = response_MFC_post_time,
    immResp_MFC_hai_pre_value = response_MFC_pre_value,
    immResp_MFC_hai_post_value = response_MFC_post_value,
    immResp_MFC_hai_MFC = response_MFC,
    immResp_MFC_hai_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_hai_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_hai_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

# Now merge those dataframes together to get the final immune response data
hipc_immResp <- list(
  max_response_MFC_df_any,
  max_response_MFC_df_nAb,
  max_response_MFC_df_hai
) %>%
  reduce(full_join, by = c("participant_id", "study_accession")) %>%
  arrange(participant_id)

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "hipc_immResp.rds")

# Save dataframe
saveRDS(hipc_immResp, file = p_save)

rm(list = ls())
