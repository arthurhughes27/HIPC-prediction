# R script to run preprocessing files in order

source_path = fs::path("analysis", "preprocessing", "young", "norm")

source(fs::path(source_path, "preprocessing_clinical.R"))

source_path = fs::path("analysis", "preprocessing", "young", "norm")

source(fs::path(source_path, "preprocessing_expression.R"))

source_path = fs::path("analysis", "preprocessing", "young", "norm")

source(fs::path(source_path, "preprocessing_GSA.R"))

source_path = fs::path("analysis", "preprocessing", "young", "norm")

source(fs::path(source_path, "preprocessing_immuneresponse.R"))

source_path = fs::path("analysis", "preprocessing", "young", "norm")

source(fs::path(source_path, "preprocessing_merging.R"))


