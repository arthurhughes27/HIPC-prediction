# File to pre-process the Blood Transcriptional Modules (BTMs) for downstream use
## This includes removing empty genesets and genesets with no genes in our data
## It also includes the engineering of new "geneset aggregates" which are groups of similar genesets

# Load necessary packages
library(fs)
library(dplyr)
library(stringr)
library(GSA)
library(readr)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data paths robustly
p_load_BTM <- fs::path(raw_data_folder, "BTM.gmt")
p_load_aggregates <- fs::path(raw_data_folder, "BTM_functional_groups.txt")
p_load_expression <- fs::path(raw_data_folder, "all_noNorm_eset.rds")

# BTMs need to be loaded in with `GSA.read.gmt()` from the `GSA` package
BTM = GSA.read.gmt(p_load_BTM)

# Aggregates are a .txt file so need to be loaded with `read.delim()` from `readr` package
BTM_functional_groups <-
  read.delim(p_load_aggregates, comment.char = "#")

# Expression data can be loaded with `readRDS()`
all_noNorm_eset <- readRDS(p_load_expression)

# First take the expression data
expr_data = all_noNorm_eset@assayData[["exprs"]]

# Now extract the gene names from the expression matrix - we are only interested in keeping these
gene_names = rownames(expr_data) %>%
  tolower()

# Now in the geneset data - make the geneset names lowercase
BTM[["genesets"]] = lapply(
  BTM[["genesets"]],
  FUN = function(x) {
    tolower(x)
  }
)

# Remove any empty elements
BTM[["genesets"]]  = lapply(
  BTM[["genesets"]],
  FUN = function(x) {
    x[x != ""]
  }
)

# Use only genesets which contain at least one of the commonly aligned genes in IS2
## Find the indices of the genesets containing at least one relevant gene
gs_membership <-
  sapply(BTM[["genesets"]], function(gs)
    any(gs %in% gene_names))

# Filter the BTM object to remove genesets without IS2 member
BTM[["genesets"]] <- BTM[["genesets"]][gs_membership]

BTM[["geneset.names"]] <- BTM[["geneset.names"]][gs_membership]

BTM[["geneset.descriptions"]] <-
  BTM[["geneset.descriptions"]][gs_membership]


# Now we preprocess the aggregate data

# Put everything in lowercase
BTM_functional_groups$SUBGROUP <-
  str_to_title(tolower(BTM_functional_groups$SUBGROUP))

# Manually replace specific terms to stay in uppercase (i.e. where uppercase represents an acronym)
BTM_functional_groups$SUBGROUP <-
  gsub("\\bEcm\\b", "ECM", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <-
  gsub("\\bDc\\b", "DC", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <-
  gsub("\\bTlr\\b", "TLR", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <-
  gsub("\\bNk\\b", "NK", BTM_functional_groups$SUBGROUP)

BTM[["geneset.aggregates"]] <-
  BTM_functional_groups$SUBGROUP[match(tolower(BTM[["geneset.names"]]), tolower(BTM_functional_groups$NAME))]

# Where a geneset does not have an aggregate, replace these empty mappings with NA
BTM[["geneset.aggregates"]] <-
  ifelse(is.na(BTM[["geneset.aggregates"]]) |
           BTM[["geneset.aggregates"]] == "", "NA", BTM[["geneset.aggregates"]])

# Rename the column names.descriptions to names
BTM[["geneset.names.descriptions"]] = BTM[["geneset.names"]]

# Remove label in parentheses from the descriptions (S or M)
BTM[["geneset.descriptions"]] <-
  sub(" \\(S[0-9]+.*\\)|\\(M[0-9]+.*\\)", "", BTM[["geneset.names"]])

# Extract label (S or M) and keep it as name
BTM[["geneset.names"]] <-
  sub(".*(\\(S[0-9]+.*\\)|\\(M[0-9]+.*\\))", "\\1", BTM[["geneset.names"]])

# Remove spaces at the end of names
BTM[["geneset.names"]] <- gsub("[()]", "", BTM[["geneset.names"]])

# Remove trailing white spaces from descriptions
BTM[["geneset.descriptions"]] <-
  trimws(BTM[["geneset.descriptions"]])

# Define an ordering for consistent figures in downstream analysis
## Roughly, this ordering is "innate immunity", "cellular processes", and "adaptive immunity"
category_order <- c(
  "Antigen Presentation",
  "Inflammatory/TLR/Chemokines",
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "DC Activation",
  "Neutrophils",
  "NK Cells",
  "Signal Transduction",
  "ECM And Migration",
  "Energy Metabolism",
  "Cell Cycle",
  "Platelets",
  "T Cells",
  "B Cells",
  "Plasma Cells",
  "NA"
)


# Assign aggregated categories as factors to genesets
BTM[["geneset.aggregates"]] <-
  factor(BTM[["geneset.aggregates"]], levels = category_order)

BTM_processed = BTM

# Save the processed geneset object

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "BTM_processed.rds")

# Save dataframe
saveRDS(BTM_processed, file = p_save)


rm(list = ls())
