# R script to run preprocessing files in order

source_path = fs::path("analysis", "preprocessing", "young", "noNorm")

source(fs::path(source_path, "preprocessing_clinical.R"))

source(fs::path(source_path, "preprocessing_expression.R"))

source(fs::path(source_path, "preprocessing_GSA.R"))

source(fs::path(source_path, "preprocessing_immuneresponse.R"))

source(fs::path(source_path, "preprocessing_merging.R"))

