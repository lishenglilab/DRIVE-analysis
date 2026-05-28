############################################################
## fig3c_performance_rank_score.R
##
## Purpose:
##   Generate Fig. 3C performance RANK_SCORE plot from one fixed
##   performance result table.
##
## Important:
##   This script assumes one final metric result per method.
##   There are no folds, no repeated runs, and no grouped settings.
##
## Rank-score rule:
##   For each metric, all methods are ranked once.
##   - RMSE: lower value is better, rank 1 is best.
##   - R2, PCC, SCC, NDCG, NWPC: higher value is better, rank 1 is best.
##
##   RANK_SCORE = number_of_methods - Rank + 1
##
##   Therefore:
##   - the best method for a metric gets the highest RANK_SCORE
##   - the worst method for a metric gets RANK_SCORE = 1
##
## Input:
##   final_metrics_evaluation_ELITE_INTERSECTION.csv
##
## Required columns:
##   Method, RMSE, PCC, SCC, R2, NDCG, NWPC
##
## Output:
##   figure3_outputs/Fig3C_performance_rank_score.pdf
##   figure3_outputs/Fig3C_performance_rank_score.png
##   figure3_outputs/Fig3C_rankscore_calculation_details.csv
##   figure3_outputs/Fig3C_method_rankscore_summary.csv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

############################################################
## 1. User settings
############################################################

## Example path. Replace this with your own project directory.
work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "final_metrics_evaluation_ELITE_INTERSECTION.csv"

out_dir <- "figure3_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metrics_cols <- c("RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC")
lower_is_better <- c("RMSE")

############################################################
## 2. Load input data
############################################################

data <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
## 3. Check and clean input
############################################################

required_cols <- c("Method", metrics_cols)
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

data <- data %>%
  select(Method, all_of(metrics_cols)) %>%
  mutate(across(all_of(metrics_cols), as.numeric))

## Enforce one fixed result per method.
## If duplicated methods are present, their metric values are averaged first.
## For the final figure, the input should ideally already contain one row per method.
data_single <- data %>%
  group_by(Method) %>%
  summarise(
    across(all_of(metrics_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

############################################################
## 4. Calculate metric-wise rank and RANK_SCORE
############################################################

data_long <- data_single %>%
  pivot_longer(
    cols = all_of(metrics_cols),
    names_to = "Metric",
    values_to = "Metric_Value"
  )

rank_details <- data_long %>%
  filter(!is.na(Metric_Value)) %>%
  group_by(Metric) %>%
  mutate(
    Better_Direction = ifelse(
      Metric %in% lower_is_better,
      "Lower is better",
      "Higher is better"
    ),
    Rank = ifelse(
      Metric %in% lower_is_better,
      rank(Metric_Value, ties.method = "min", na.last = "keep"),
      rank(-Metric_Value, ties.method = "min", na.last = "keep")
    )
  ) %>%
  ungroup()

number_of_methods <- length(unique(rank_details$Method))

rank_details <- rank_details %>%
  mutate(
    Number_of_Methods = number_of_methods,
    RANK_SCORE = Number_of_Methods - Rank + 1
  )

############################################################
## 5. Summarise method-level rank score
############################################################

method_summary <- rank_details %>%
  group_by(Method) %>%
  summarise(
    Median_RANK_SCORE = median(RANK_SCORE, na.rm = TRUE),
    Mean_RANK_SCORE = mean(RANK_SCORE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Median_RANK_SCORE), desc(Mean_RANK_SCORE))

method_order <- method_summary$Method

rank_details <- rank_details %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Metric = factor(Metric, levels = c("NDCG", "NWPC", "PCC", "R2", "RMSE", "SCC"))
  )

method_summary <- method_summary %>%
  mutate(Method = factor(Method, levels = method_order))

write.csv(
  rank_details,
  file.path(out_dir, "Fig3C_rankscore_calculation_details.csv"),
  row.names = FALSE
)

write.csv(
  method_summary,
  file.path(out_dir, "Fig3C_method_rankscore_summary.csv"),
  row.names = FALSE
)

############################################################
## 6. Metric colors
############################################################

metric_colors <- c(
  "NDCG" = "#66A61E",
  "NWPC" = "#E6AB02",
  "PCC"  = "#D95F02",
  "R2"   = "#E7298A",
  "RMSE" = "#1B9E77",
  "SCC"  = "#7570B3"
)

############################################################
## 7. Plot Fig. 3C
############################################################

p <- ggplot() +

  ## Grey summary stems
  geom_segment(
    data = method_summary,
    aes(
      x = Method,
      xend = Method,
      y = 0,
      yend = Median_RANK_SCORE
    ),
    color = "#C8C8C8",
    linewidth = 3.2,
    lineend = "butt"
  ) +

  ## Grey median summary circles
  geom_point(
    data = method_summary,
    aes(x = Method, y = Median_RANK_SCORE),
    size = 8,
    shape = 16,
    color = "#C8C8C8",
    alpha = 0.95
  ) +

  ## Colored metric points
  geom_point(
    data = rank_details,
    aes(
      x = Method,
      y = RANK_SCORE,
      color = Metric
    ),
    shape = 18,
    size = 2.8,
    alpha = 0.95,
    position = position_jitter(width = 0.14, height = 0)
  ) +

  scale_color_manual(values = metric_colors, name = "Metric") +

  scale_y_continuous(
    limits = c(0, number_of_methods + 1),
    breaks = seq(0, number_of_methods, by = 5)
  ) +

  labs(
    x = "Methods",
    y = "Performance RANK_SCORE"
  ) +

  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "#F2F2F2", color = NA),
    panel.background = element_rect(fill = "#F2F2F2", color = NA),
    legend.background = element_rect(fill = "#F2F2F2", color = NA),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "#333333",
      size = 9
    ),
    axis.text.y = element_text(color = "#333333"),
    axis.title = element_text(color = "#333333", size = 14),
    axis.line = element_line(color = "#333333", linewidth = 0.6),
    axis.ticks = element_line(color = "#333333"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

############################################################
## 8. Export
############################################################

ggsave(
  filename = file.path(out_dir, "Fig3C_performance_rank_score.pdf"),
  plot = p,
  width = 12.5,
  height = 6.8,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Fig3C_performance_rank_score.png"),
  plot = p,
  width = 12.5,
  height = 6.8,
  dpi = 300
)

############################################################
## Finished
############################################################

message("Fig. 3C rank-score plot finished.")
message("Output directory: ", normalizePath(out_dir))
message("Output files:")
message("  - Fig3C_performance_rank_score.pdf")
message("  - Fig3C_performance_rank_score.png")
message("  - Fig3C_rankscore_calculation_details.csv")
message("  - Fig3C_method_rankscore_summary.csv")
