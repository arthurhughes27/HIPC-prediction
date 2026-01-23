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
p_load_young_noNorm <- fs::path(raw_data_folder, "young_noNorm_eset.rds")

# Read in the rds file
young_noNorm_eset <- readRDS(p_load_young_noNorm)

# Load the expression data
young_noNorm_expr = young_noNorm_eset@assayData[["exprs"]] %>% 
  t() %>% 
  as.data.frame()

# Make the column names lowercase
colnames(young_noNorm_expr) = young_noNorm_expr %>% 
  colnames() %>% 
  tolower()

# Now extract information for the first two columns (participant_id and study_time_collected)
sample_info_young_noNorm = rownames(young_noNorm_expr)

# Use `stringr` and regular expressions to extract the participant_id and study_time_collected from the unique identifiers
matches_young_noNorm <- str_match(sample_info_young_noNorm , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")

# Participant ids
participant_id_young_noNorm <- matches_young_noNorm[, 2]

# Study times (numeric, rounded)
study_time_collected_young_noNorm <- matches_young_noNorm[, 3] %>% 
  as.numeric() %>% 
  round(2)

# Insert the identifying information as the first two columns 
young_noNorm_expr <- young_noNorm_expr %>%
  mutate(participant_id = participant_id_young_noNorm, study_time_collected = study_time_collected_young_noNorm) %>%
  select(participant_id, study_time_collected, everything())

# Save processed dataframes

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save_young_noNorm <- fs::path(processed_data_folder, "young_noNorm_expr.rds")

# Save dataframe
saveRDS(young_noNorm_expr, file = p_save_young_noNorm)

rm(list = ls())
