# R script to run all data engineering files

source_path = fs::path("analysis", "engineering")

source(fs::path(source_path, "data_engineering_all_noNorm.R"))

source(fs::path(source_path, "data_engineering_all_norm.R"))

source(fs::path(source_path, "data_engineering_young_noNorm.R"))

source(fs::path(source_path, "data_engineering_young_norm.R"))

