# R script to interpret results of the univariate feature selection algorithm

library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidytext)   # reorder_within(), scale_y_reordered()
library(stringr)    # str_wrap()

output_folder = fs::path("output", "results")
vi_figures_folder = fs::path("output", "figures", "feature selection")

dataset_of_interest = "all_noNorm"

vaccine_of_interest = "Influenza (IN)"

study_of_interest = "SDY1276"

assay_of_interest = "nAb"

gender_of_interest = "none"

model_of_interest = "lm"

metric_of_interest = "R.spearman"


file_name = paste0(
  "feature_selection_",
  gsub("[[:space:]()]", "", tolower(vaccine_of_interest)),
  "_",
  dataset_of_interest,
  "_",
  study_of_interest,
  "_",
  assay_of_interest,
  "_",
  gender_of_interest,
  "_",
  model_of_interest,
  ".rds"
)

p_load <- fs::path(output_folder, file_name)

# Load the results
res = readRDS(p_load)

extract_pred_res <- function(res,
                             data.sel = NULL,
                             feat.eng.col = NULL,
                             feat.eng.row = NULL) {
  
  res %>%
    keep(~ (is.null(data.sel)     || .x$data.sel %in% data.sel) &
           (is.null(feat.eng.col) || .x$feat.eng.col %in% feat.eng.col) &
           (is.null(feat.eng.row) || .x$feat.eng.row %in% feat.eng.row)) %>%
    
    map_dfr(~ .x$pred.res %>%
              mutate(
                data.sel = .x$data.sel,
                feat.eng.col = .x$feat.eng.col,
                feat.eng.row = .x$feat.eng.row
              ))
}

df <- extract_pred_res(
  res,
  data.sel = c("d0", "d1", "d3"),
  feat.eng.col = "none",
  feat.eng.row = "mean"
)

# =========================
# Parameters
# =========================
top.n <- 15
plot_title <- "Top predictors by timepoint"
wrap_width <- 25    # number of characters per line for pred names

# =========================
# Prepare data
# =========================
res_top <- df %>%
  group_by(data.sel) %>%
  slice_max(order_by = metric, n = top.n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # map fold-change labels to the same day-group for colouring
    day_group = case_when(
      data.sel == "d0"               ~ "d0",
      data.sel %in% c("d1", "d1fc")  ~ "d1",
      data.sel %in% c("d3", "d3fc")  ~ "d3",
      TRUE                                 ~ as.character(data.sel)
    ),
    # wrap pred names for readability
    predictor_label = str_wrap(pred, width = wrap_width),
    # reorder predictors within each facet by metric
    predictor_reordered = reorder_within(
      predictor_label, metric, data.sel
    )
  )

# =========================
# Fixed x-axis limits across facets
# =========================
x_min <- min(res_top$metric, na.rm = TRUE)
x_max <- max(res_top$metric, na.rm = TRUE)
pad <- max(0.02 * (x_max - x_min), 0.1)

# =========================
# Colour palette (red → orange → pale yellow)
# =========================
day_cols <- c(
  d0 = "#b30000",   # deep red
  d1 = "#ff7f00",   # orange (used for d1 and d1fc)
  d3 = "#ffd966"    # pale yellow
)

# =========================
# Plot
# =========================
p <- ggplot(res_top, aes(x = metric, y = predictor_reordered, fill = day_group)) +
  geom_col(
    color = "grey20",
    size = 0.2,
    width = 0.75
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.6,
    colour = "black"
  ) +
  facet_wrap(
    ~ data.sel,
    scales = "free_y",
    nrow = 1
  ) +
  scale_y_reordered() +
  scale_fill_manual(values = day_cols, name = "Day") +
  coord_cartesian(xlim = c(x_min - pad, x_max + pad)) +
  labs(
    title = plot_title,
    x = paste0("Relative Percent Gain in ", metric_of_interest),
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 25, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 30),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# =========================
# Draw
# =========================
print(p)

file_name_figure = paste0(
  "feature_selection_",
  gsub("[[:space:]()]", "", tolower(vaccine_of_interest)),
  "_",
  dataset_of_interest,
  "_",
  study_of_interest,
  "_",
  assay_of_interest,
  "_",
  gender_of_interest,
  "_",
  model_of_interest,
  ".pdf"
)

ggsave(
  filename = file_name_figure,
  path = vi_figures_folder,
  plot = p,
  width = 45,
  height = 30,
  units = "cm"
)

