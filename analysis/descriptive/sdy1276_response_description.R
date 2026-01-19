# R script to describe the antibody responses to TIV in SDY1276

# Libraries
library(ggplot2)
library(tidyverse)
library(fs)
library(forcats)

processed_data_folder = "data"
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Use fs::path() to specify the data path robustly
p_load <- fs::path(processed_data_folder, "raw_response_influenzain.rds")

# Save dataframe
response = readRDS(p_load)

# Filter the dataframe by study
response = response %>% 
  filter(study_accession == "SDY1276") %>% 
  filter(!is.na(value_preferred))

# Figure 1 - Sample sizes for neutralising antibody and HAI assays across time
sample_counts <- response %>%
  # collapse multiple strain measurements for the same participant/assay/timepoint
  group_by(assay, study_time_collected, participant_id) %>%
  summarise(has_measure = n() >= 1, .groups = "drop") %>%
  # count participants with at least one measurement per assay x timepoint
  group_by(assay, study_time_collected) %>%
  summarise(sample_size = sum(has_measure), .groups = "drop")



# Tile plot
p1 <- ggplot(sample_counts, aes(x = as.character(study_time_collected), y = assay)) +
  geom_tile(fill = "#2ca02c", colour = "grey", size = 0.3) +
  geom_text(aes(label = sample_size), size = 6) +
  labs(x = "Timepoint",
       y = "Assay",
       title = "Availability of immune repsonse measurements") +
  theme_minimal(base_size = 25) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(),
        plot.title = element_text(face = "bold"))

# Print the plot
print(p1)

ggsave(
  filename = "sdy1276_assay_availability.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 30,
  height = 15,
  units = "cm"
)

# Figure 2 : individual antibody trajectories according to assay and strain

# Minimal preprocessing (assume value_preferred exists)
df_plot <- response %>%
  # normalize assay names to two canonical levels: HAI (row 1) and nAb (row 2)
  mutate(assay = case_when(
    tolower(assay) == "hai" ~ "HAI",
    tolower(assay)  == "nab" ~ "nAb"
  ),
  log_value_preferred = log2(value_preferred)) %>%
  # force row order HAI then nAb (if both present)
  mutate(assay = factor(assay, levels = c("HAI", "nAb"))) 

# Plot: spaghetti plots facetted with assays as rows and strains as columns
p2 = ggplot(df_plot, aes(x = as.character(study_time_collected), y = log_value_preferred, group = participant_id)) +
  geom_line(alpha = 0.2, size = 0.4) +
  facet_grid(rows = vars(assay), cols = vars(response_strain_analyte)) +
  labs(
    x = "Timepoint",
    y = "log2(Titre)",
    title = "Individual antibody trajectories"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 20)
  )

print(p2)

ggsave(
  filename = "sdy1276_antibody_trajectories.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 30,
  height = 15,
  units = "cm"
)

# Plot 3 - MFC versus pre-existing antibody titres

## Derive MFC
response_MFC <- response %>%
  filter(study_time_collected %in% c(0, 28)) %>%
  # # collapse duplicate rows per participant/assay/strain/time (use mean)
  group_by(participant_id, assay, response_strain_analyte, study_time_collected) %>%
  # make day0 and day28 columns
  pivot_wider(names_from = study_time_collected, values_from = value_preferred, names_prefix = "t") %>%
  # compute fold-change per strain (day0 / day28) ; set NA if denominator missing or <= 0
  mutate(fold_change = if_else(is.na(t0) | is.na(t28) | t0 <= 0,
                               NA_real_,
                               t28 / t0)) %>%
  # for each participant_id x assay take the maximum fold_change across strains
  group_by(participant_id, assay) %>%
  summarise(
    MFC = if (all(is.na(fold_change))) NA_real_ else max(fold_change, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  filter(!is.na(MFC)) %>% 
  mutate(log2_MFC = log2(MFC))

p3 = ggplot(response_MFC, aes(x = assay, y = log2_MFC)) +
  geom_violin(width = 0.6) +
  geom_jitter(alpha = 0.3) +
  labs(
    x = "Assay",
    y = "log2(MFC) at day 28",
    title = "Distribution of maximum fold-change by assay"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p3

ggsave(
  filename = "sdy1276_MFC_violin.pdf",
  path = descriptive_figures_folder,
  plot = p3,
  width = 20,
  height = 10,
  units = "cm"
)

# Distributions of fold-change per assay and strain

# derive fold-change per participant × assay × strain
response_FC_strain <- response %>%
  filter(study_time_collected %in% c(0, 28)) %>%
  # collapse duplicates (mean) for safety
  group_by(participant_id, assay, response_strain_analyte, study_time_collected) %>%
  summarise(value = mean(value_preferred, na.rm = TRUE), .groups = "drop") %>%
  # widen to t0 / t28
  pivot_wider(names_from = study_time_collected, values_from = value, names_prefix = "t") %>%
  # compute fold-change as t28 / t0; set NA for missing or non-positive denominators
  mutate(
    fold_change = if_else(is.na(t0) | is.na(t28) | t0 <= 0, NA_real_, t28 / t0),
    log2_fold_change = if_else(is.na(fold_change), NA_real_, log2(fold_change))
  ) %>%
  # keep relevant columns
  select(participant_id, assay, response_strain_analyte, fold_change, log2_fold_change)

# ---- Plot: violins in a grid (rows = assays, cols = strains), identical scales ----
# create ordered assay × strain factor
response_FC_strain <- response_FC_strain %>%
  mutate(
    assay = factor(assay, levels = c("hai", "nAb")),
    assay_strain = interaction(assay, response_strain_analyte, sep = " : ", lex.order = TRUE)
  )

p4 <- ggplot(response_FC_strain,
                    aes(x = assay_strain, y = log2_fold_change)) +
  geom_boxplot(width = 0.8, na.rm = TRUE, outliers = F) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, na.rm = TRUE) +
  labs(
    x = "Assay × Strain",
    y = "log2(fold-change) (Day28 / Day0)",
    title = "Per-strain fold-change distributions"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15)
  )

p4

ggsave(
  filename = "sdy1276_MFC_strain_boxplot.pdf",
  path = descriptive_figures_folder,
  plot = p4,
  width = 29,
  height = 18,
  units = "cm"
)

