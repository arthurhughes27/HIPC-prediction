# R script to interpret results from the prediction task
library(dplyr)
library(ggplot2)
library(tidyr)

# Path to save the results
results_folder = fs::path("output", "results")
figures_folder = fs::path("output", "figures", "comparison")
p_load <- fs::path(results_folder, "metrics_transformations_ranger_allData.rds")
# Save the results
metrics_df = readRDS(p_load)

levels_order <- c("d0","d1","d1fc","d3","d3fc")
row_order <- c("mean", "median", "max", "iqr", "cv", "pc1", "ssgsea")

baseline_R2 <- metrics_df %>%
  filter(data.selection == "baseline") %>%
  pull(R2) %>% .[1]

df <- metrics_df %>%
  filter(data.selection %in% levels_order) %>%
  mutate(
    data.selection = factor(data.selection, levels = levels_order),
    feature.engineering.row = factor(feature.engineering.row, levels = row_order)
  )

y_lim <- c(0,1)

p1 <- ggplot(
  df %>% filter(include.covariates == TRUE),
  aes(x = data.selection, y = R2, fill = feature.engineering.row)
) +
  geom_col(width = 0.6, position = position_dodge2(width = 0.9, padding = 0.2)) +
  geom_hline(aes(yintercept = baseline_R2, linetype = "Baseline"), color = "black") +
  scale_linetype_manual(name = "", values = c("Baseline" = "dashed")) +
  scale_fill_brewer(palette = "Set2", drop = FALSE) +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Random Forest – R² by data selection and row transformation (covariates included)",
    x = "Predictor Set",
    y = "R²",
    fill = "Row transformation"
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(size = 17))

p2 <- ggplot(
  df %>% filter(include.covariates == FALSE),
  aes(x = data.selection, y = R2, fill = feature.engineering.row)
) +
  geom_col(width = 0.6, position = position_dodge2(width = 0.9, padding = 0.2)) +
  geom_hline(aes(yintercept = baseline_R2, linetype = "Baseline"), color = "black") +
  scale_linetype_manual(name = "", values = c("Baseline" = "dashed")) +
  scale_fill_brewer(palette = "Set2", drop = FALSE) +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Random Forest – R² by data selection and row transformation (covariates excluded)",
    x = "Predictor Set",
    y = "R²",
    fill = "Row transformation"
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(size = 17))

print(p1)
print(p2)


ggsave(
  filename = "comparison_transformations_ranger_allData_CovTRUE.pdf",
  path = figures_folder,
  plot = p1,
  width = 25,
  height = 15,
  units = "cm"
)

ggsave(
  filename = "comparison_transformations_ranger_allData_CovFALSE.pdf",
  path = figures_folder,
  plot = p2,
  width = 25,
  height = 15,
  units = "cm"
)

p3 <- ggplot(
  df %>% 
    filter(include.covariates == TRUE, feature.engineering.row == "mean"),
  aes(x = data.selection, y = R2)
) +
  geom_col(width = 0.6) +
  geom_hline(
    aes(yintercept = baseline_R2, linetype = "Baseline"),
    color = "black"
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Baseline" = "dashed")
  ) +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Random Forest – R² across predictor sets (mean transformation, covariates included)",
    x = "Predictor Set",
    y = "R²"
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(size = 13))

print(p3)

ggsave(
  filename = "comparison_mean_ranger_allData_CovTRUE.pdf",
  path = figures_folder,
  plot = p3,
  width = 20,
  height = 15,
  units = "cm"
)
