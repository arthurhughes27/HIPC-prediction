# R script to interpret results from the prediction task for HAI in SDY1276 in both males and females
library(dplyr)
library(ggplot2)
library(tidyr)

# Path to save the results
results_folder = fs::path("output", "results")
figures_folder = fs::path("output", "figures", "comparison")
p_load_male <- fs::path(results_folder, "metrics_transformations_lm_all_noNorm_SDY1276_hai_Male.rds")
p_load_female <- fs::path(results_folder, "metrics_transformations_lm_all_noNorm_SDY1276_hai_Female.rds")

# Read the results
metrics_df_male = readRDS(p_load_male) %>% 
  mutate(gender = "Male")
metrics_df_female = readRDS(p_load_female) %>% 
  mutate(gender = "Female")

metrics_df_all= bind_rows(metrics_df_male, metrics_df_female)
