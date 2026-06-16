############################################################
## supfig_s2_model_performance_heatmap.R
## Output: Sup2.pdf, Sup2.png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "model_metric.csv"

## -------------------------------------------------------------------------
## 1. Read and prepare data
## -------------------------------------------------------------------------
df <- read.csv(input_file, check.names = FALSE)

## Define metric order, mode order, and method order (bottom to top)
metric_order <- c("NDCG", "NWPC", "PCC", "RMSE", "SCC")
mode_order <- c("Cell-blind", "Drug-blind", "Mixed")

method_order <- c(
  "ML baseline",
  "BANDRP",
  "DeepAEG",
  "DeepCCDS",
  "DeepCDR",
  "DeepTTA",
  "DIPK",
  "GADRP",
  "GPDRP_GAT",
  "GPDRP_GCN",
  "GPDRP_GIN",
  "GPDRP_Trans",
  "GraphDRP-GINConvNet",
  "GraphDRP_GAT_GCN",
  "GraphDRP_GATNet",
  "GraphDRP_GCNNet",
  "NERD",
  "paccmann",
  "Precily"
)

## -------------------------------------------------------------------------
## 2. Reshape to long format, average replicates
## -------------------------------------------------------------------------
plot_df <- df %>%
  select(Method, mode, all_of(metric_order)) %>%
  pivot_longer(
    cols = all_of(metric_order),
    names_to = "Metric",
    values_to = "value"
  ) %>%
  group_by(Method, mode, Metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

## -------------------------------------------------------------------------
## 3. Min-max normalization (RMSE reversed)
## -------------------------------------------------------------------------
plot_df <- plot_df %>%
  group_by(Metric) %>%
  mutate(
    norm_score = case_when(
      max(value, na.rm = TRUE) == min(value, na.rm = TRUE) ~ 1,
      Metric == "RMSE" ~ (max(value, na.rm = TRUE) - value) /
        (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)),
      TRUE ~ (value - min(value, na.rm = TRUE)) /
        (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
    )
  ) %>%
  ungroup() %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Metric = factor(Metric, levels = metric_order),
    mode   = factor(mode, levels = mode_order),
    label  = sprintf("%.2f", norm_score)
  )

## Optionally remove ML baseline (uncomment if needed)
# plot_df <- plot_df %>% filter(Method != "ML baseline")

## -------------------------------------------------------------------------
## 4. Generate heatmap
## -------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = Metric, y = Method, fill = norm_score)) +
  geom_tile(
    color = "white",
    linewidth = 1.15,
    width = 0.92,
    height = 0.92
  ) +
  geom_text(
    aes(label = label),
    color = "white",
    fontface = "bold",
    size = 3.4
  ) +
  facet_grid(
    . ~ mode,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = number_format(accuracy = 0.01),
    name = "norm_score"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 18,
      color = "grey15",
      margin = margin(b = 8)
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 12,
      color = "grey30"
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 10.5,
      color = "grey30"
    ),
    panel.spacing.x = unit(0.45, "cm"),
    legend.title = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    legend.text = element_text(size = 10, color = "black"),
    legend.key.height = unit(2.0, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p)

## -------------------------------------------------------------------------
## 5. Save output
## -------------------------------------------------------------------------
ggsave(
  filename = "Sup2.pdf",
  plot = p,
  width = 12,
  height = 9,
  units = "in"
)

ggsave(
  filename = "Sup2.png",
  plot = p,
  width = 12,
  height = 8.8,
  units = "in",
  dpi = 300
)

message("Supplementary Fig. S2 generated successfully: Sup2.pdf and Sup2.png")