# File to pre-process the HIPC IS2 clinical characteristics data for use
## This includes some data engineering on columns (e.g. collapse "unknown" and "Not Specified" gender to same value)
## and also the creation of some data for downstream analysis purposes (e.g. the associating a colour with each vaccine)

# Load necessary packages
library(fs)
library(Biobase)
library(dplyr)
library(forcats)
library(stringr)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data path robustly
p_load_young_norm <- fs::path(raw_data_folder, "young_norm_eset.rds")

# Read in the rds file
young_norm_eset <- readRDS(p_load_young_norm)

# Extract the clinical data as a dataframe
hipc_clinical_young_norm = young_norm_eset@phenoData@data %>% 
  as.data.frame()

# Data engineering
# "collapsing" the values from certain important columns
# For gender and race, "Not Specified" and "Unknown" represent the same thing
# For ethnicity, the values "Not Hispanic or Latino" and "Other" represent the same thing (the only other values are
# "Not Specified" and "Hispanic or Latino").
# In addition, we round numerical time columns to avoid long recurring numbers
hipc_clinical_young_norm <- hipc_clinical_young_norm %>%
  mutate(
    gender    = fct_collapse(factor(gender), Unknown = c("Not Specified", "Unknown")),
    race      = fct_collapse(factor(race), Unknown = c("Not Specified", "Unknown")),
    ethnicity = fct_collapse(factor(ethnicity), Other   = c("Not Hispanic or Latino", "Other")),
    study_time_collected = round(study_time_collected, 2),
    time_post_last_vax = round(time_post_last_vax, 2)
  )

# Abbreviate the vaccine type column
hipc_clinical_young_norm$vaccine_type <- hipc_clinical_young_norm$vaccine_type %>%
  as.factor() %>%
  recode_factor(
    "Conjugate" = "CJ",
    "Inactivated" = "IN",
    "Inactivated/Recombinant protein" = "IN/RP",
    "Live virus" = "LV",
    "Polysaccharide" = "PS",
    "Recombinant Viral Vector" = "RVV",
    "Recombinant protein" = "RP"
  )

# Create vaccine name column by combining pathogen and vaccine type
hipc_clinical_young_norm <- hipc_clinical_young_norm %>%
  mutate(vaccine_name = str_c(pathogen, " (", vaccine_type, ")"))

# Create a shortened vaccine name variable
hipc_clinical_young_norm$vaccine_name_short <- recode(
  hipc_clinical_young_norm$vaccine_name,
  "Ebola (RVV)"           = "Ebola (RVV)",
  "Yellow Fever (LV)"     = "Y.F. (LV)",
  "Smallpox (LV)"         = "Smallpox (LV)",
  "Tuberculosis (RVV)"    = "T.B. (RVV)",
  "Hepatitis A/B (IN/RP)" = "Hep. (IN/RP)",
  "Meningococcus (CJ)"    = "Men. (CJ)",
  "Meningococcus (PS)"    = "Men. (PS)",
  "Malaria (RP)"          = "Malaria (RP)",
  "Influenza (IN)"        = "Inf. (IN)",
  "Influenza (LV)"        = "Inf. (LV)",
  "Pneumococcus (PS)"     = "Pneumo. (PS)",
  "Varicella Zoster (LV)" = "Varicella (LV)",
  "HIV (RVV)"             = "HIV (RVV)"
)

# Define an ordering for the vaccines (this is for later to make figures consistent)
conditions_order <- c(
  "Tuberculosis (RVV)",
  "Varicella Zoster (LV)",
  "Yellow Fever (LV)",
  "Ebola (RVV)",
  "Hepatitis A/B (IN/RP)",
  "HIV (RVV)",
  "Influenza (IN)",
  "Influenza (LV)",
  "Malaria (RP)",
  "Meningococcus (CJ)",
  "Meningococcus (PS)",
  "Pneumococcus (PS)",
  "Smallpox (LV)"
)

# Assign this order to the vaccine names
hipc_clinical_young_norm <- hipc_clinical_young_norm %>% 
  mutate(vaccine_name = factor(vaccine_name, levels = conditions_order))

# Same for the shortened names
conditions_order_short <- c(
  "T.B. (RVV)",
  "Varicella (LV)",
  "Y.F. (LV)",
  "Ebola (RVV)",
  "Hep. (IN/RP)",
  "HIV (RVV)",
  "Inf. (IN)",
  "Inf. (LV)",
  "Malaria (RP)",
  "Men. (CJ)",
  "Men. (PS)",
  "Pneumo. (PS)",
  "Smallpox (LV)"
)


# Assign this order to the shortened vaccine names
hipc_clinical_young_norm <- hipc_clinical_young_norm %>% 
  mutate(vaccine_name_short = factor(vaccine_name_short, levels = conditions_order_short))

# Define a colour for each vaccine
## This colour palette was chosen to maximise visual distinctiveness for 13 vaccines
## Using the "iwanthue" tool (https://medialab.github.io/iwanthue/)
color_palette_vaccine = c(
  "#b94a73",
  "#c6aa3c",
  "#6f71d9",
  "#64c46a",
  "#be62c2",
  "#7d973c",
  "#563382",
  "#4ea76e",
  "#bc69b0",
  "#33d4d1",
  "#bb4c41",
  "#6a87d3",
  "#b57736"
)

# Write a helper function to assign the colours
assign_color <- function(vaccine_name) {
  return(color_palette_vaccine[match(hipc_clinical_young_norm$vaccine_name,
                             levels(hipc_clinical_young_norm$vaccine_name))])
}

# Assign the colours to the vaccine names
hipc_clinical_young_norm$vaccine_colour <-
  assign_color(hipc_clinical_young_norm$vaccine_name)

# Define an ordering for the studies (this is for later to make figures consistent)
study_descriptions_young_norm = hipc_clinical_young_norm %>% 
  dplyr::select(vaccine_name, study_accession) %>% 
  distinct() %>% 
  arrange(vaccine_name)

# There are three studies which include multiple vaccines 
# SDY1260 and SDY1325 for Meningococcus CJ and PS, 
# SDY269 for influenza IN and LV, 
# SDY180 for influenza IN and Pneumococcus (PS)
# Define a new variable renaming these to "SDYXa" (CJ, IN) and "SDYXb" (PS, LV) where X is replaced with the appropriate study number
hipc_clinical_young_norm = hipc_clinical_young_norm %>%
  mutate(study_accession_unique = study_accession) 

hipc_clinical_young_norm = hipc_clinical_young_norm %>% 
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1260" &
        vaccine_name == "Meningococcus (CJ)",
      "SDY1260a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1260" &
        vaccine_name == "Meningococcus (PS)",
      "SDY1260b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1325" &
        vaccine_name == "Meningococcus (CJ)",
      "SDY1325a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1325" &
        vaccine_name == "Meningococcus (PS)",
      "SDY1325b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY269" &
        vaccine_name == "Influenza (IN)",
      "SDY269a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY269" &
        vaccine_name == "Influenza (LV)",
      "SDY269b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY180" &
        vaccine_name == "Influenza (IN)",
      "SDY180a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY180" &
        vaccine_name == "Pneumococcus (PS)",
      "SDY180b",
      study_accession_unique
    )
  )

# Redefine the study order
study_descriptions_young_norm = hipc_clinical_young_norm %>% 
  dplyr::select(vaccine_name, study_accession_unique) %>% 
  distinct() %>% 
  arrange(vaccine_name)

study_order <- study_descriptions_young_norm$study_accession_unique

# Assign this order to the vaccine names
hipc_clinical_young_norm <- hipc_clinical_young_norm %>% 
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order))

# Define a colour for each study
## This colour palette was chosen to represent 30 studies across 13 vaccines
## The study colours are chosen to be a gradient within each vaccine based on the previously chosen vaccine colours 
## Using the "supercolorpalette" tool (https://supercolorpalette.com)
color_palette_study = c(
  "#b94a73", # Tuberculosis (RVV)
  "#c6aa3c", # Varicella Zoster (LV)
  "#F1D986",
  "#4042FB",  # Yellow Fever (LV)
  "#6077FC",
  "#889FFD",
  "#B4C4FD",
  "#64c46a", # Ebola (RVV)
  "#be62c2", # Hepatitis A/B (IN/RP)
  "#7d973c", # HIV (RVV)
  "#170931",# Influenza (IN)
  "#230C48",
  "#2E0D5F", 
  "#3A0A77",
  "#460090",
  "#5100A5",
  "#5C00BA",
  "#6701D0", 
  "#7301E5", #
  "#7E01FB",
  "#8632FF",
  "#8D4DFF",
  "#9562FF",
  "#9D74FF",
  "#A684FF", 
  "#AF94FF",
  "#4ea76e", # Influenza (LV)
  "#bc69b0", # Malaria (RP) 
  "#33d4d1", # Meningococcus (CJ)
  "#B2FBF8",
  "#bb4c41", # Meningococcus (PS)
  "#F39C90",
  "#6a87d3", # Pneumococcus (PS)
  "#b57736"  # Smallpox (LV)
)

# Write a helper function to assign the colours
assign_color <- function(study_accession_unique) {
  return(color_palette_study[match(hipc_clinical_young_norm$study_accession_unique,
                                     levels(hipc_clinical_young_norm$study_accession_unique))])
}

# Assign the colours to the vaccine names
hipc_clinical_young_norm$study_colour <-
  assign_color(hipc_clinical_young_norm$study_accession_unique)

# Final dataframe to be saved has samples as rows and variables as columns
dim(hipc_clinical_young_norm)

# Save processed dataframe

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save_young_norm <- fs::path(processed_data_folder, "hipc_clinical_young_norm.rds")

# Save dataframe
saveRDS(hipc_clinical_young_norm, file = p_save_young_norm)

rm(list = ls())
