# File for preprocessing the raw immune response data files
# We do this for all flu studies
# Output is a dataframe in long format giving the participant, study, strain, assay (HAI or nAb),
# timepoint and value

# ---------------------------
# Load packages
# ---------------------------
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)
library(lme4)

# ---------------------------
# Paths / input file locations
# ---------------------------
raw_data_folder = "data-raw"
processed_data_folder = "data"

p_load_nAb <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_hai <- fs::path(raw_data_folder, "hai_2025-01-10_01-13-41.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical_all_noNorm.rds")

# ---------------------------
# Read raw data
# ---------------------------
response_hai = read_excel(p_load_hai) %>%
  clean_names()

response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

hipc_clinical_all_noNorm = readRDS(p_load_clinical)

# ---------------------------
# Filter to participants with gene-expression data
# ---------------------------
participants = hipc_clinical_all_noNorm %>%
  pull(participant_id) %>%
  unique()

response_hai = response_hai %>%
  filter(participant_id %in% participants)

response_nAb = response_nAb %>%
  filter(participant_id %in% participants)

# ---------------------------
# Standardize column names, tag assay, and merge HAI + nAb into one table
# ---------------------------
response_hai = response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai")

response_nAb = response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb")

raw_response_influenzain = bind_rows(response_nAb, response_hai) %>%
  select(-gender) %>%
  arrange(participant_id) %>%
  distinct()

# ---------------------------
# Add study/clinical metadata (age, gender, study_accession, etc.)
# ---------------------------
hipc_studies = hipc_clinical_all_noNorm %>%
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

raw_response_influenzain = merge(
  x = raw_response_influenzain,
  y = hipc_studies,
  by = "participant_id",
  all = F
) %>%
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

# Specific study adjustment: transform SDY1276 values (same as original)
raw_response_influenzain = raw_response_influenzain %>%
  mutate(value_preferred = if_else(
    study_accession == "SDY1276",
    4 ^ value_preferred,
    value_preferred
  ))

# ---------------------------
# Save merged raw responses (unfiltered timepoints)
# ---------------------------
p_save <- fs::path(processed_data_folder, "raw_response_influenzain_all_noNorm.rds")
saveRDS(raw_response_influenzain, file = p_save)

# ---------------------------
# Filter to pre-vax and ~day-28 post-vax (day 28 ± 7 days)
# ---------------------------
raw_response_influenzain = raw_response_influenzain %>%
  filter((study_time_collected < 36 & study_time_collected > 20) | study_time_collected <= 0)

# ---------------------------
# Prepare pre- and post- measurement tables
# - pre_vax: most recent pre (study_time_collected <= 0), per participant/assay/strain
# - post_vax: post measurements (study_time_collected > 0)
# ---------------------------
pre_vax <- raw_response_influenzain %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, assay, response_strain_analyte) %>%
  slice_max(study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    response_MFC_pre_value = value_preferred,
    response_MFC_pre_time = study_time_collected
  ) %>%
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_pre_time, response_MFC_pre_value)

post_vax <- raw_response_influenzain %>%
  filter(study_time_collected > 0) %>%
  rename(
    response_MFC_post_value = value_preferred,
    response_MFC_post_time = study_time_collected
  ) %>%
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_post_time, response_MFC_post_value)

# ---------------------------
# Merge pre/post into wide-format per (participant, assay, strain)
# ---------------------------
merged_vax <- full_join(
  pre_vax,
  post_vax,
  by = c("participant_id", "study_accession", "response_strain_analyte", "assay"),
  relationship = "many-to-many"
) %>%
  select(
    participant_id,
    study_accession,
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value
  ) %>%
  filter(!is.na(response_MFC_pre_value), !is.na(response_MFC_post_value))

# ---------------------------
# Compute fold change and log2 versions (same edge-case handling as original)
# ---------------------------
response_MFC_df <- merged_vax %>%
  mutate(
    response_MFC = response_MFC_post_value / response_MFC_pre_value,
    # If pre is zero (division -> Inf), set MFC to post value
    response_MFC = ifelse(is.infinite(response_MFC), response_MFC_post_value, response_MFC),
    # If both pre and post are zero, set fold change to 1
    response_MFC = ifelse(response_MFC_pre_value == 0 & response_MFC_post_value == 0, 1, response_MFC)
  ) %>%
  mutate(
    response_log2_MFC = ifelse(is.infinite(log2(response_MFC)), NA, log2(response_MFC)),
    response_log2_MFC_pre_value = log2(response_MFC_pre_value),
    response_log2_MFC_post_value = log2(response_MFC_post_value)
  ) %>%
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

# ---------------------------
# Summaries per participant
# - max across assays (any assay)
# - max within assay
# - mean within assay (across strains)
# ---------------------------
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

mean_response_nAb <- response_MFC_df %>%
  filter(assay == "nAb") %>%
  group_by(participant_id, study_accession) %>%
  summarise(
    immResp_mean_nAb_pre_value  = mean(response_MFC_pre_value,  na.rm = TRUE),
    immResp_mean_nAb_post_value = mean(response_MFC_post_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    immResp_mean_nAb_log2_pre_value = log2(immResp_mean_nAb_pre_value),
    immResp_mean_nAb_log2_post_value = log2(immResp_mean_nAb_post_value),
    immResp_mean_nAb_FC = immResp_mean_nAb_post_value/immResp_mean_nAb_pre_value,
    immResp_mean_nAb_log2_FC = log2(immResp_mean_nAb_FC)
  )

mean_response_hai <- response_MFC_df %>%
  filter(assay == "hai") %>%
  group_by(participant_id, study_accession) %>%
  summarise(
    immResp_mean_hai_pre_value  = mean(response_MFC_pre_value,  na.rm = TRUE),
    immResp_mean_hai_post_value = mean(response_MFC_post_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    immResp_mean_hai_log2_pre_value = log2(immResp_mean_hai_pre_value),
    immResp_mean_hai_log2_post_value = log2(immResp_mean_hai_post_value),
    immResp_mean_hai_FC = immResp_mean_hai_post_value/immResp_mean_hai_pre_value,
    immResp_mean_hai_log2_FC = log2(immResp_mean_hai_FC)
  )

# ---------------------------
# Combine participant-level immune-response summaries
# ---------------------------
hipc_immResp <- list(
  max_response_MFC_df_any,
  max_response_MFC_df_nAb,
  max_response_MFC_df_hai,
  mean_response_nAb,
  mean_response_hai
) %>%
  reduce(full_join, by = c("participant_id", "study_accession")) %>%
  arrange(participant_id)

# ---------------------------
# Prepare dataset for LMM estimation (keep non-missing values)
# ---------------------------
base_df <- raw_response_influenzain %>%
  filter(!is.na(value_preferred))

# ---------------------------
# Function to compute participant-level LMM-derived slope (per study)
# ---------------------------
compute_LMM <- function(df_sub) {
  df_sub <- df_sub %>%
    rename(
      participant_id = participant_id,
      strain = response_strain_analyte,
      day = study_time_collected,
      y = value_preferred
    ) %>%
    mutate(t28 = day / 28)
  
  # remove participants with fewer than 2 measurements
  keep_pids <- df_sub %>%
    group_by(participant_id) %>%
    summarise(n = n()) %>%
    filter(n >= 2) %>%
    pull(participant_id)
  
  df_sub <- df_sub %>% filter(participant_id %in% keep_pids)
  
  # skip studies with <2 participants or <2 strains
  if(length(unique(df_sub$participant_id)) < 2) return(NULL)
  if(length(unique(df_sub$strain)) < 2) return(NULL)
  
  fit <- lmer(y ~ strain * t28 + (1 + t28 | participant_id) + (1 | participant_id:strain),
              data = df_sub, REML = TRUE)
  
  re_pid <- ranef(fit)$participant_id
  beta_hat <- fixef(fit)["t28"]
  
  data.frame(
    participant_id = rownames(re_pid),
    study_accession = unique(df_sub$study_accession),
    immResp_LMM = as.numeric(beta_hat + re_pid[,"t28"]),
    row.names = NULL
  )
}

# ---------------------------
# Apply LMM per study for HAI and nAb assays
# ---------------------------
hai_df <- base_df %>%
  filter(assay == "hai") %>%
  group_split(study_accession) %>%
  purrr::map_dfr(compute_LMM) %>%
  rename(immResp_LMM_hai = immResp_LMM)

nab_df <- base_df %>%
  filter(assay == "nAb") %>%
  group_split(study_accession) %>%
  purrr::map_dfr(compute_LMM) %>%
  rename(immResp_LMM_nAb = immResp_LMM)

immResp_LMM_df <- full_join(hai_df, nab_df, by = c("participant_id", "study_accession"))

hipc_immResp = hipc_immResp %>%
  full_join(y = immResp_LMM_df, by = c("participant_id", "study_accession"))

# ---------------------------
# Compute strain-wise z-scores, then participant-level mean/max of z-scores
# Produces columns:
# immResp_mean_nAb_std_pre_value, immResp_mean_nAb_std_post_value,
# immResp_max_nAb_std_pre_value,  immResp_max_nAb_std_post_value,
# immResp_mean_hai_std_pre_value, immResp_mean_hai_std_post_value,
# immResp_max_hai_std_pre_value,  immResp_max_hai_std_post_value
# ---------------------------
library(tidyr)

# 1) pivot pre/post into long time column and z-score within (study x assay x strain x time)
std_long <- response_MFC_df %>%
  select(participant_id, study_accession, assay, response_strain_analyte,
         response_MFC_pre_value, response_MFC_post_value) %>%
  pivot_longer(
    cols = c(response_MFC_pre_value, response_MFC_post_value),
    names_to = "time",
    values_to = "value"
  ) %>%
  mutate(time = ifelse(time == "response_MFC_pre_value", "pre", "post")) %>%
  group_by(study_accession, assay, response_strain_analyte, time) %>%
  mutate(
    grp_mean = mean(value, na.rm = TRUE),
    grp_sd   = sd(value, na.rm = TRUE),
    value_z  = ifelse(is.na(value), NA_real_,
                      ifelse(is.na(grp_sd) | grp_sd == 0, 0, (value - grp_mean) / grp_sd))
  ) %>%
  ungroup() %>%
  select(participant_id, study_accession, assay, response_strain_analyte, time, value_z)

# 2) summarise per participant/study/assay/time: mean and max of z-scores across strains
std_summary <- std_long %>%
  group_by(participant_id, study_accession, assay, time) %>%
  summarise(
    mean_std = if(all(is.na(value_z))) NA_real_ else mean(value_z, na.rm = TRUE),
    max_std  = if(all(is.na(value_z))) NA_real_ else max(value_z, na.rm = TRUE),
    .groups = "drop"
  )

# 3) split into separate tibbles for each (transformation x assay x time) with requested names
immResp_mean_nAb_std_pre  <- std_summary %>%
  filter(assay == "nAb", time == "pre") %>%
  rename(immResp_mean_nAb_std_pre_value  = mean_std)

immResp_mean_nAb_std_post <- std_summary %>%
  filter(assay == "nAb", time == "post") %>%
  rename(immResp_mean_nAb_std_post_value = mean_std)

immResp_max_nAb_std_pre   <- std_summary %>%
  filter(assay == "nAb", time == "pre") %>%
  rename(immResp_max_nAb_std_pre_value   = max_std)

immResp_max_nAb_std_post  <- std_summary %>%
  filter(assay == "nAb", time == "post") %>%
  rename(immResp_max_nAb_std_post_value  = max_std)

immResp_mean_hai_std_pre  <- std_summary %>%
  filter(assay == "hai", time == "pre") %>%
  rename(immResp_mean_hai_std_pre_value  = mean_std)

immResp_mean_hai_std_post <- std_summary %>%
  filter(assay == "hai", time == "post") %>%
  rename(immResp_mean_hai_std_post_value = mean_std)

immResp_max_hai_std_pre   <- std_summary %>%
  filter(assay == "hai", time == "pre") %>%
  rename(immResp_max_hai_std_pre_value   = max_std)

immResp_max_hai_std_post  <- std_summary %>%
  filter(assay == "hai", time == "post") %>%
  rename(immResp_max_hai_std_post_value  = max_std)

# 4) join the eight tibbles into one table keyed by participant_id + study_accession,
#    keeping only the final eight value columns plus keys
std_cols <- list(
  immResp_mean_nAb_std_pre,
  immResp_mean_nAb_std_post,
  immResp_max_nAb_std_pre,
  immResp_max_nAb_std_post,
  immResp_mean_hai_std_pre,
  immResp_mean_hai_std_post,
  immResp_max_hai_std_pre,
  immResp_max_hai_std_post
) %>%
  purrr::reduce(function(x, y) full_join(x, y, by = c("participant_id", "study_accession"))) %>%
  select(
    participant_id,
    study_accession,
    immResp_mean_nAb_std_pre_value,
    immResp_mean_nAb_std_post_value,
    immResp_max_nAb_std_pre_value,
    immResp_max_nAb_std_post_value,
    immResp_mean_hai_std_pre_value,
    immResp_mean_hai_std_post_value,
    immResp_max_hai_std_pre_value,
    immResp_max_hai_std_post_value
  )

# 5) merge into hipc_immResp (preserve existing cols; append only the 8 standardized columns)
hipc_immResp <- hipc_immResp %>%
  left_join(std_cols, by = c("participant_id", "study_accession")) %>%
  arrange(participant_id)

# ---------------------------
# Save final dataset
# ---------------------------
p_save <- fs::path(processed_data_folder, "hipc_immResp.rds")
saveRDS(hipc_immResp, file = p_save)

# cleanup
rm(list = ls())