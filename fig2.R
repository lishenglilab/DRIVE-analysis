############################################################
## fig2_three_setting_benchmark.R
##
## Purpose:
##   Generate the Fig. 2 benchmark plots for Mixed, Drug-blind,
##   and Cell-blind settings from one external Excel table.
##
## Input:
##   synthetic_benchmark_settings.xlsx
##
## Required columns:
##   1. Method
##   2. RMSE
##   3. PCC
##   4. SCC
##   5. R2
##   6. NDCG
##   7. NWPC
##   8. Setting
##
## Output:
##   1. Fig2_three_settings/fig2_three_setting_long.csv
##   2. Fig2_three_settings/fig2_three_setting_summary.csv
##   3. Fig2_three_settings/Fig2_model_palette.csv
##   4. Fig2_three_settings/Fig2_*_benchmark.pdf
##   5. Fig2_three_settings/Fig2_*_benchmark.png
############################################################

## =========================
## 0. Packages
## =========================

pkgs <- c("readxl", "tidyverse", "patchwork")

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running this script."
  )
}

library(readxl)
library(tidyverse)
library(patchwork)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

work_dir <- get_script_dir()
setwd(work_dir)

## =========================
## 1. Input / output
## =========================

input_file <- "synthetic_benchmark_settings.xlsx"
out_dir <- "Fig2_three_settings"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

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
    Method = trimws(as.character(Method)),
    Setting = as.character(Setting),
    Method = str_replace_all(Method, "[\u2010-\u2015\u2212]", "-")
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

model_palette <- c(
  "BANDRP" = "#DB9D85",
  "DeepAEG" = "#CFA373",
  "DeepCDR" = "#ABB065",
  "DeepTTA" = "#93B66E",
  "DIPK" = "#79BA7E",
  "GraphDRP" = "#ACA4E2",
  "GPDRP" = "#41BEA7",
  "NERD" = "#E093C3",
  "PaccMann" = "#E2969A",
  "Precily" = "#E494AF",
  "GADRP" = "#5CBB92",
  "DeepCCDS" = "#BFA967"
)

errorbar_color <- "gray30"

write_csv(
  tibble(Method = names(model_palette), Color = unname(model_palette)),
  file.path(out_dir, "Fig2_model_palette.csv")
)

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
      legend.position = "right",
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
      DisplayOrder = factor(seq_len(n()), levels = seq_len(n())),
      Method = factor(Method, levels = names(model_palette))
    )

  ggplot(plot_df, aes(x = DisplayOrder, y = Mean, fill = Method)) +
    geom_col(
      width = 0.82,
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
      y = y_lab,
      fill = "Model"
    ) +
    scale_fill_manual(values = model_palette, drop = FALSE) +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_fig2_panel()
}

plot_one_setting <- function(setting_name, data_summary) {
  p1 <- plot_metric_panel(setting_name, "NDCG", data_summary, show_y_title = TRUE)
  p2 <- plot_metric_panel(setting_name, "NWPC", data_summary, show_y_title = FALSE)
  p3 <- plot_metric_panel(setting_name, "PCC", data_summary, show_y_title = FALSE)

  p4 <- plot_metric_panel(setting_name, "R2", data_summary, show_y_title = TRUE)
  p5 <- plot_metric_panel(setting_name, "RMSE", data_summary, show_y_title = FALSE)
  p6 <- plot_metric_panel(setting_name, "SCC", data_summary, show_y_title = FALSE)

  (p1 + p2 + p3) / (p4 + p5 + p6) +
    plot_layout(guides = "collect") +
    plot_annotation(title = setting_name) &
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    ) &
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
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
p_drug <- plot_one_setting("Drug-blind", summary_df)
p_cell <- plot_one_setting("Cell-blind", summary_df)

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
