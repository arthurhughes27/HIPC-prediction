# File for preprocessing gene expression data from IS2 for downstream use
# The end result should be a dataframe with columns as variables and rows as samples
# The first two columns will give the participant identifier and time of sample, allowing for unique identification


# Load necessary packages
library(fs)
library(dplyr)
library(stringr)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data paths robustly
p_load_young_norm <- fs::path(raw_data_folder, "young_norm_eset.rds")

# Read in the rds file
young_norm_eset <- readRDS(p_load_young_norm)

# Load the expression data
young_norm_expr = young_norm_eset@assayData[["exprs"]] %>% 
  t() %>% 
  as.data.frame()

# Make the column names lowercase
colnames(young_norm_expr) = young_norm_expr %>% 
  colnames() %>% 
  tolower()

# Now extract information for the first two columns (participant_id and study_time_collected)
sample_info_young_norm = rownames(young_norm_expr)

# Use `stringr` and regular expressions to extract the participant_id and study_time_collected from the unique identifiers
matches_young_norm <- str_match(sample_info_young_norm , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")

# Participant ids
participant_id_young_norm <- matches_young_norm[, 2]

# Study times (numeric, rounded)
study_time_collected_young_norm <- matches_young_norm[, 3] %>% 
  as.numeric() %>% 
  round(2)

# Insert the identifying information as the first two columns 
young_norm_expr <- young_norm_expr %>%
  mutate(participant_id = participant_id_young_norm, study_time_collected = study_time_collected_young_norm) %>%
  select(participant_id, study_time_collected, everything())

# Save processed dataframes

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save_young_norm <- fs::path(processed_data_folder, "young_norm_expr.rds")

# Save dataframe
saveRDS(young_norm_expr, file = p_save_young_norm)

rm(list = ls())
