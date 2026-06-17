############################################################
## Read one Excel file and draw three Fig. 2-style plots
## Rules:
## 1) Each metric is ordered from high to low
## 2) No method or model names are displayed in the figure
## 3) English only
############################################################

## =========================
## 0. Packages
## =========================

pkgs <- c("readxl", "tidyverse", "patchwork")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(readxl)
library(tidyverse)
library(patchwork)

## =========================
## 1. Input / output
## =========================

excel_candidates <- list.files(
  pattern = "\\.(xls|xlsx)$",
  ignore.case = TRUE,
  full.names = TRUE
)

if (length(excel_candidates) == 0) {
  stop("No Excel file was found in the working directory.")
}

input_file <- excel_candidates[1]
out_dir <- "Fig2_three_settings"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2. Read data
## =========================

raw_df <- read_excel(input_file)

colnames(raw_df) <- colnames(raw_df) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.", "_")

required_cols <- c("Method", "RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC", "Setting")

missing_cols <- setdiff(required_cols, colnames(raw_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df_long <- raw_df %>%
  mutate(
    Method = as.character(Method),
    Setting = as.character(Setting),
    Method = str_replace_all(Method, "鈭?", "-")
  ) %>%
  pivot_longer(
    cols = c(RMSE, PCC, SCC, R2, NDCG, NWPC),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Setting = case_when(
      str_detect(tolower(Setting), "mixed") ~ "Mixed",
      str_detect(tolower(Setting), "drug") ~ "Drug-blind",
      str_detect(tolower(Setting), "cell") ~ "Cell-blind",
      TRUE ~ Setting
    ),
    Metric = factor(
      Metric,
      levels = c("NDCG", "NWPC", "PCC", "R2", "RMSE", "SCC")
    ),
    Setting = factor(
      Setting,
      levels = c("Mixed", "Drug-blind", "Cell-blind")
    )
  ) %>%
  filter(!is.na(Value))

## =========================
## 3. Summary statistics
## =========================

summary_df <- df_long %>%
  group_by(Setting, Method, Metric) %>%
  summarise(
    n = n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    SEM = ifelse(is.na(SEM), 0, SEM)
  )

write_csv(df_long, file.path(out_dir, "fig2_three_setting_long.csv"))
write_csv(summary_df, file.path(out_dir, "fig2_three_setting_summary.csv"))

## =========================
## 4. Plot settings
## =========================

bar_fill <- "#6FA8DC"
errorbar_color <- "gray30"

theme_fig2_panel <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = base_size),
      axis.title.y = element_text(color = "black", size = base_size + 3),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
      strip.text = element_text(color = "black", size = base_size + 1),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      legend.position = "none",
      plot.title = element_text(size = base_size + 5, face = "bold", hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
}

## =========================
## 5. Plot function
## =========================

plot_metric_panel <- function(setting_name, metric_name, data_summary, show_y_title = FALSE) {
  y_lab <- if (show_y_title) "Mean Value" else ""

  plot_df <- data_summary %>%
    filter(
      Setting == setting_name,
      Metric == metric_name
    ) %>%
    arrange(desc(Mean)) %>%
    mutate(
      DisplayOrder = factor(seq_len(n()), levels = seq_len(n()))
    )

  ggplot(plot_df, aes(x = DisplayOrder, y = Mean)) +
    geom_col(
      width = 0.82,
      fill = bar_fill,
      color = NA,
      alpha = 0.95
    ) +
    geom_errorbar(
      aes(ymin = Mean - SEM, ymax = Mean + SEM),
      width = 0.22,
      linewidth = 0.45,
      color = errorbar_color
    ) +
    labs(
      x = NULL,
      y = y_lab
    ) +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_fig2_panel()
}

plot_one_setting <- function(setting_name, data_summary) {
  p1 <- plot_metric_panel(setting_name, "NDCG", data_summary, show_y_title = TRUE)
  p2 <- plot_metric_panel(setting_name, "NWPC", data_summary, show_y_title = FALSE)
  p3 <- plot_metric_panel(setting_name, "PCC",  data_summary, show_y_title = FALSE)

  p4 <- plot_metric_panel(setting_name, "R2",   data_summary, show_y_title = TRUE)
  p5 <- plot_metric_panel(setting_name, "RMSE", data_summary, show_y_title = FALSE)
  p6 <- plot_metric_panel(setting_name, "SCC",  data_summary, show_y_title = FALSE)

  (p1 + p2 + p3) / (p4 + p5 + p6) +
    plot_layout(guides = "collect") +
    plot_annotation(title = setting_name) &
    theme(
      legend.position = "none",
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
}

## =========================
## 6. Generate three plots
## =========================

settings_to_plot <- c("Mixed", "Drug-blind", "Cell-blind")

for (s in settings_to_plot) {
  p <- plot_one_setting(s, summary_df)
  file_stub <- str_replace_all(s, "-", "_")

  ggsave(
    filename = file.path(out_dir, paste0("Fig2_", file_stub, "_benchmark.pdf")),
    plot = p,
    width = 14,
    height = 9,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(out_dir, paste0("Fig2_", file_stub, "_benchmark.png")),
    plot = p,
    width = 14,
    height = 9,
    units = "in",
    dpi = 600,
    limitsize = FALSE
  )
}

## =========================
## 7. Optional previews
## =========================

p_mixed <- plot_one_setting("Mixed", summary_df)
p_drug  <- plot_one_setting("Drug-blind", summary_df)
p_cell  <- plot_one_setting("Cell-blind", summary_df)

ggsave(
  filename = file.path(out_dir, "Fig2_Mixed_benchmark_preview.pdf"),
  plot = p_mixed,
  width = 14,
  height = 9,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "Fig2_Drug_blind_benchmark_preview.pdf"),
  plot = p_drug,
  width = 14,
  height = 9,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "Fig2_Cell_blind_benchmark_preview.pdf"),
  plot = p_cell,
  width = 14,
  height = 9,
  units = "in",
  limitsize = FALSE
)
