# R script to interpret results of the univariate feature selection algorithm

library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidytext)   # reorder_within(), scale_y_reordered()
library(stringr)    # str_wrap()

output_folder = fs::path("output", "results")
vi_figures_folder = fs::path("output", "figures", "variable selection")
p_load <- fs::path(output_folder, "feature.selection.lm.mean.all.rds")
# Save the results
res = readRDS(p_load)

# =========================
# Parameters
# =========================
top.n <- 25
plot_title <- "Top predictors by data selection"
wrap_width <- 25    # number of characters per line for predictor names

# =========================
# Prepare data
# =========================
res_top <- res %>%
  group_by(data.selection) %>%
  slice_max(order_by = metric, n = top.n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # map fold-change labels to the same day-group for colouring
    day_group = case_when(
      data.selection == "d0"               ~ "d0",
      data.selection %in% c("d1", "d1fc")  ~ "d1",
      data.selection %in% c("d3", "d3fc")  ~ "d3",
      TRUE                                 ~ as.character(data.selection)
    ),
    # wrap predictor names for readability
    predictor_label = str_wrap(predictor, width = wrap_width),
    # reorder predictors within each facet by metric
    predictor_reordered = reorder_within(
      predictor_label, metric, data.selection
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
    ~ data.selection,
    scales = "free_y",
    nrow = 1
  ) +
  scale_y_reordered() +
  scale_fill_manual(values = day_cols, name = "Day") +
  coord_cartesian(xlim = c(x_min - pad, x_max + pad)) +
  labs(
    title = plot_title,
    x = "Relative Percent Gain in sRMSE",
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

ggsave(
  filename = "varImp_univariate.pdf",
  path = vi_figures_folder,
  plot = p,
  width = 45,
  height = 40,
  units = "cm"
)

