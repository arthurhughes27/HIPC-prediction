library(tidyverse)
library(RColorBrewer)   # for Brewer palettes
library(fs)

# Path for results 
p_results = fs::path("output", "results", "metrics_transformations_ranger.rds")
figures_folder = fs::path("output", "figures", "comparison")
metrics_df = readRDS(p_results)

# ---- factor order for x-axis ----
order_sel <- c("baseline","d0", "d1", "d1fc", "d3", "d3fc")

metrics_df <- metrics_df %>%
  mutate(data.selection = factor(data.selection, levels = order_sel))

metric = "sRMSE"
feat.eng.col = "z"

# ---- Baseline R2 ----------------------------------------------------------
baseline_metric <- metrics_df %>%
  filter(data.selection == "baseline") %>%
  pull(metric) %>%
  unique()

plot_df <- metrics_df %>%
  filter(data.selection != "baseline",
         feature.engineering.col == feat.eng.col) %>%
  mutate(
    # keep x-order consistent
    data.selection = factor(data.selection, levels = order_sel[order_sel != "baseline"]),
    # use row transform only for fill/legend
    row_trans = factor(feature.engineering.row, levels = unique(feature.engineering.row))
  )

# ---- Plot -----------------------------------------------------------------

p1 <- ggplot(plot_df, aes(x = data.selection, y = .data[[metric]], fill = row_trans)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, colour = "grey30") +
  geom_hline(yintercept = baseline_metric, linetype = "dashed", size = 0.8, colour = "black") +
  # baseline annotation placed on the right, make sure it's not clipped
  annotate(
    "text",
    x = Inf, y = baseline_metric,
    label = paste0("Baseline ", metric," = ", format(round(baseline_metric, 3), nsmall = 3)),
    hjust = 1.05, vjust = -0.5, size = 4.5
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "GS Aggregation") +
  labs(
    title = paste0("Cross-validated ",metric ," by predictor set (column transform: ", feat.eng.col, ")"),
    x = "Predictor set",
    y = metric
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    # allow annotation to be drawn outside plotting area on the right
    plot.margin = margin(t = 8, r = 40, b = 8, l = 8)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  coord_cartesian(ylim = c(0, 1), clip = "off")   # ensures annotation and dashed line remain visible

# Print plot
p1

# Optionally save
ggsave(
  filename = paste0("metrics_comparison_ranger_col",feat.eng.col,"_", metric,".pdf"),
  path = figures_folder,
  plot = p1,
  width = 27,
  height = 13,
  units = "cm"
)
