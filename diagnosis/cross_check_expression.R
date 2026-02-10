# File to cross-check expression data with theodora's work
library(tidyverse)
library(readxl)
# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data paths robustly
p_load_all_norm <- fs::path(raw_data_folder, "all_norm_eset.rds")
p_load_all_norm2 <- fs::path(raw_data_folder, "all_norm_withResponse_eset.rds")
# Read in the rds file
all_norm_eset <- readRDS(p_load_all_norm)
all_norm_eset2 <- readRDS(p_load_all_norm2)

# Load the expression data
all_norm_expr = all_norm_eset@assayData[["exprs"]] %>% 
  t() %>% 
  as.data.frame()

# Load the expression data
all_norm_expr2 = all_norm_eset2@assayData[["exprs"]] %>% 
  t() %>% 
  as.data.frame()

# Make the column names lowercase
colnames(all_norm_expr) = all_norm_expr %>% 
  colnames() %>% 
  tolower()

# Make the column names lowercase
colnames(all_norm_expr2) = all_norm_expr2 %>% 
  colnames() %>% 
  tolower()

# Now extract information for the first two columns (participant_id and study_time_collected)
sample_info_all_norm = rownames(all_norm_expr)
sample_info_all_norm2 = rownames(all_norm_expr2)

# Use `stringr` and regular expressions to extract the participant_id and study_time_collected from the unique identifiers
matches_all_norm <- str_match(sample_info_all_norm , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")
matches_all_norm2 <- str_match(sample_info_all_norm2 , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")

# Participant ids
participant_id_all_norm <- matches_all_norm[, 2]
participant_id_all_norm2 <- matches_all_norm2[, 2]

# Study times (numeric, rounded)
study_time_collected_all_norm <- matches_all_norm[, 3] %>% 
  as.numeric() %>% 
  round(2)

study_time_collected_all_norm2 <- matches_all_norm2[, 3] %>% 
  as.numeric() %>% 
  round(2)

# Insert the identifying information as the first two columns 
all_norm_expr <- all_norm_expr %>%
  mutate(participant_id = participant_id_all_norm, study_time_collected = study_time_collected_all_norm) %>%
  select(participant_id, study_time_collected, everything())

all_norm_expr2 <- all_norm_expr2 %>%
  mutate(participant_id = participant_id_all_norm2, study_time_collected = study_time_collected_all_norm2) %>%
  select(participant_id, study_time_collected, everything())

pids = all_norm_expr2 %>% 
  pull(participant_id) %>% 
  unique()

all_norm_expr = all_norm_expr %>% 
  filter(participant_id %in% pids)

all_norm_expr  <- all_norm_expr[, colSums(is.na(all_norm_expr)) == 0]
all_norm_expr2 <- all_norm_expr2[, colSums(is.na(all_norm_expr2)) == 0]

sum(all_norm_expr[,5] != all_norm_expr2[,5])
# Response dataframe and non-response data frame have equal values

# Now check against theodora's 

# Loading the data
## Define the paths
raw_data_path = fs::path("data-raw")
path_ge = fs::path(raw_data_path, "all_norm_withResponse_eset.rds")
path_nab = fs::path(raw_data_path, "neut_ab_titer_2025-01-10_01-13-22.xlsx")

# Load the data
data = readRDS(path_ge)
neut.ab.assay <- read_excel(path_nab)

# Extract the clinical data
rawdata <- as.data.frame(data@phenoData@data)
# Extract the gene expression data
ge.data <- as.data.frame(data@assayData[["exprs"]])
rm(data)

# Preprocessing
## Select only females with 2008 TIV vaccine with gene expression info at baseline and 1 day, from SDY1276
rawdata <- rawdata %>%
  filter(
    pathogen == "Influenza",
    vaccine == "TIV (2008)",
    gender == "Female",
    study_time_collected %in% c(0, 1),
    study_accession == "SDY1276"
  ) %>%
  rename(expression_time = study_time_collected) %>%
  group_by(participant_id) %>% # keep only participants with Both 0 and 1
  filter(n_distinct(expression_time) == 2) %>%
  ungroup()

# Extract participant ids from filtered data
filtered_ids = rawdata %>%
  pull(participant_id) %>%
  unique()

# Filter particpants from neutralizing antibody dataframe
## Only choose participants with antibody levels at day 0 and 28 without missing values
nab_filtered <- neut.ab.assay %>%
  filter(
    `Participant ID` %in% filtered_ids,
    `Study Time Collected` %in% c(0, 28),
    !is.na(`Value Preferred`)
  )

# Grab the average of the antibody response
nab_composite <- nab_filtered %>%
  group_by(`Participant ID`, `Study Time Collected`) %>%
  summarize(response = mean(`Value Preferred`, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(`Participant ID`) %>%
  filter(n_distinct(`Study Time Collected`) == 2) %>%
  ungroup()

# Match expression days with response days
## We match d0 gene expression levels with d0 Ab levels, and d1 gene expression levels with d28 Ab levels
nab_final <- nab_composite %>%
  mutate(
    StudyTimeMatch = case_when(
      `Study Time Collected` == 0 ~ 0,
      `Study Time Collected` %in% c(28) ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(StudyTimeMatch))

# Merge dataframes
eligible.participants <- nab_final %>%
  inner_join(
    rawdata %>% dplyr::select("participant_id", "study_accession", "uid", "expression_time"),
    by = c("Participant ID" = "participant_id", "StudyTimeMatch" = "expression_time")
  )

# Create a final data frame with the "treatment status"
filtered.assay.data <- eligible.participants %>%
  dplyr::select(`Participant ID`,
                uid,
                `Study Time Collected`,
                response,
                study_accession) %>%
  rename(`Antibody Value` = response) %>%
  mutate(Treatment = ifelse(`Study Time Collected` == 0, 0, 1))

# Change UID from columns to rows
ge.data.transposed <- ge.data %>%
  rownames_to_column(var = "uid") %>%  # Make row names a column
  pivot_longer(-uid, names_to = "gene", values_to = "expression") %>% # Convert wide to long
  pivot_wider(names_from = uid, values_from = expression) %>% # Convert back to wide with uid as row names
  rename(uid = gene)
rm(ge.data)

# Merge gene expression data and assay data
processed.data <- filtered.assay.data %>%
  left_join(ge.data.transposed, by = "uid")

uids = processed.data %>% 
  pull(uid)


all_norm_expr = all_norm_expr %>% 
  mutate(uid = rownames(all_norm_expr)) %>% 
  relocate(uid, .before = "participant_id") %>% 
  filter(uid %in% uids)

check_df1 = all_norm_expr %>% 
  select(uid, zzz3)

rownames(check_df1) = NULL

check_df2 = processed.data %>% 
  select(uid, ZZZ3)

sum(check_df1[,2] != check_df2[,2])
