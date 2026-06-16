############################################################
## fig4c_drug_phase_donut.R
## Output: concentric_donut_polyline_final.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "fig4c_comparison_data.csv"
out_dir <- "fig4c_drug_phase_donut"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Required columns: drug_name, Max_Phase, actual_log_ic50, predicted_ic50
comparison_data <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("drug_name", "Max_Phase", "actual_log_ic50", "predicted_ic50")
missing_cols <- setdiff(required_cols, colnames(comparison_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

comparison_data <- comparison_data %>%
  mutate(
    actual_log_ic50 = as.numeric(actual_log_ic50),
    predicted_ic50 = as.numeric(predicted_ic50),
    Max_Phase = ifelse(is.na(Max_Phase) | Max_Phase == "", "Unknown", Max_Phase)
  )

actual_summary <- comparison_data %>%
  filter(!is.na(actual_log_ic50), actual_log_ic50 < 0) %>%
  count(Max_Phase, name = "N") %>%
  mutate(Source = "DepMap measured")

predicted_summary <- comparison_data %>%
  filter(!is.na(predicted_ic50), predicted_ic50 < 0) %>%
  count(Max_Phase, name = "N") %>%
  mutate(Source = "DRIVE predicted")

summary_df <- bind_rows(actual_summary, predicted_summary) %>%
  group_by(Source) %>%
  mutate(Proportion = N / sum(N)) %>%
  ungroup()

write.csv(summary_df, file.path(out_dir, "fig4c_drug_phase_summary.csv"), row.names = FALSE)

custom_colors <- c(
  "Approved" = "#4CAF50",
  "Phase 3" = "#2196F3",
  "Phase 2" = "#FFC107",
  "Phase 1" = "#FF9800",
  "Early Phase 1" = "#1f3b73",
  "Preclinical" = "#F44336",
  "Unknown" = "#9E9E9E"
)

summary_df <- summary_df %>%
  arrange(Source, desc(Max_Phase)) %>%
  group_by(Source) %>%
  mutate(
    ymax = cumsum(Proportion),
    ymin = ymax - Proportion,
    LabelY = (ymax + ymin) / 2,
    x_position = ifelse(Source == "DepMap measured", 2.0, 3.0),
    x_label = ifelse(Source == "DepMap measured", 1.2, 3.9),
    hjust = ifelse(Source == "DepMap measured", 1, 0),
    x_segment_start = ifelse(Source == "DepMap measured", 1.85, 3.15),
    x_segment_end = ifelse(Source == "DepMap measured", 1.35, 3.75)
  ) %>%
  ungroup()

p <- ggplot(summary_df, aes(x = x_position, y = Proportion, fill = Max_Phase)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.4) +
  geom_segment(
    aes(x = x_segment_start, y = LabelY, xend = x_segment_end, yend = LabelY),
    color = "grey50", linewidth = 0.35, inherit.aes = FALSE
  ) +
  geom_text(
    aes(x = x_label, y = LabelY,
        label = paste0(Max_Phase, "\nN=", format(N, big.mark = ",")),
        hjust = hjust),
    size = 3.2, lineheight = 0.9, inherit.aes = FALSE
  ) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = custom_colors, name = "Drug phase", na.value = "grey80") +
  xlim(c(0, 4.6)) +
  labs(title = "Proportion of Drug Phases (log(IC50) < 0)") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

ggsave(file.path(out_dir, "concentric_donut_polyline_final.pdf"),
       p, width = 14, height = 12, dpi = 300)

message("Fig. 4C donut plot finished.")
