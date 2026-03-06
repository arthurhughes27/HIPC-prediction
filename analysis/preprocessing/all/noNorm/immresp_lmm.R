# File for preprocessing the raw immune response data files 
# We do this for all flu studies 
# Output is a dataframe in long format giving the participant, study, strain, assay (HAI or nAb),
# timepoint and value

library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)
library(lme4)

raw_data_folder = "data-raw"
processed_data_folder = "data"

p_load_nAb <- fs::path(processed_data_folder, "raw_response_influenzain_all_noNorm.rds")
raw_response_influenzain = readRDS(p_load_nAb)

base_df <- raw_response_influenzain %>% 
  filter(!is.na(value_preferred))

compute_LMM <- function(df_sub) {
  df_sub <- df_sub %>%
    rename(participant_id = participant_id,
           strain = response_strain_analyte,
           day = study_time_collected,
           y = value_preferred) %>%
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
