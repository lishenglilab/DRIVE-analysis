############################################################
## fig2_three_setting_benchmark_final.R
##
## Read model_metric.csv, remove ML baseline, and generate
## three 12:9 benchmark PDFs with raw-value scatter points.
############################################################

pkgs <- c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "patchwork")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

setwd(get_script_dir())
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "model_metric.csv"
out_dir <- "Fig2_three_settings"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

raw_df <- read_csv(input_file, show_col_types = FALSE)
names(raw_df) <- names(raw_df) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.", "_")

required_cols <- c("Method", "RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC", "mode")
missing_cols <- setdiff(required_cols, names(raw_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

metric_levels <- c("NDCG", "NWPC", "PCC", "R2", "RMSE", "SCC")
setting_levels <- c("Mixed", "Drug-blind", "Cell-blind")

df_long <- raw_df %>%
  transmute(
    Method = trimws(as.character(Method)),
    Setting = trimws(as.character(mode)),
    across(all_of(metric_levels), as.numeric)
  ) %>%
  filter(!str_detect(str_to_lower(Method), "^ml(?:\\s+baseline)?$")) %>%
  pivot_longer(all_of(metric_levels), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Setting = case_when(
      str_detect(str_to_lower(Setting), "mixed") ~ "Mixed",
      str_detect(str_to_lower(Setting), "drug") ~ "Drug-blind",
      str_detect(str_to_lower(Setting), "cell") ~ "Cell-blind",
      TRUE ~ Setting
    ),
    Metric = factor(Metric, levels = metric_levels),
    Setting = factor(Setting, levels = setting_levels)
  ) %>%
  filter(!is.na(Value), !is.na(Setting), nzchar(Method))

summary_df <- df_long %>%
  group_by(Setting, Method, Metric) %>%
  summarise(
    n = n(), Mean = mean(Value), SD = sd(Value), SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(SEM = if_else(is.na(SEM), 0, SEM))

write_csv(df_long, file.path(out_dir, "fig2_three_setting_long.csv"))
write_csv(summary_df, file.path(out_dir, "fig2_three_setting_summary.csv"))

## Colors remain attached to model names, even when bars are sorted by value.
model_palette <- c(
  "GraphDRP-GINConvNet" = "#8EACE0",
  "GraphDRP_GAT_GCN" = "#C49CDD",
  "GraphDRP_GATNet" = "#ACA4E2",
  "GraphDRP_GCNNet" = "#D596D2",
  "DeepTTA" = "#79BA7E",
  "BANDRP" = "#DB9D85",
  "NERD" = "#E093C3",
  "GADRP" = "#5CBD92",
  "GPDRP_GAT" = "#41BEA7",
  "GPDRP_GIN" = "#4CB9CC",
  "GPDRP_GCN" = "#38BDBB",
  "GPDRP_Trans" = "#6CB4D9",
  "DIPK" = "#CFA373",
  "Precily" = "#E494AF",
  "paccmann" = "#E2979A",
  "DeepAEG" = "#BFAA67",
  "DeepCDR" = "#93B66E",
  "DeepCCDS" = "#ABB065"
)

legend_order <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP_GAT_GCN", "GraphDRP_GATNet", "GraphDRP_GCNNet", "GraphDRP-GINConvNet",
  "GPDRP_GAT", "GPDRP_GCN", "GPDRP_GIN", "GPDRP_Trans",
  "NERD", "paccmann", "Precily", "GADRP", "DeepCCDS"
)
legend_order <- legend_order[legend_order %in% unique(df_long$Method)]
write_csv(tibble(Method = legend_order, Color = unname(model_palette[legend_order])),
          file.path(out_dir, "Fig2_model_palette.csv"))

theme_fig2_panel <- function(base_size = 13) {
  panel_aspect <- (50.354 + 3.223) / 31.175
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(color = "black", size = base_size + 1),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1.1),
      strip.text = element_text(color = "black", size = base_size + 1, face = "plain"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.1),
      aspect.ratio = panel_aspect,
      legend.position = "right", plot.margin = margin(4, 4, 4, 4)
    )
}

plot_metric_panel <- function(setting_name, metric_name, data_summary, data_long,
                              show_y_title = FALSE) {
  plot_df <- data_summary %>%
    filter(Setting == setting_name, Metric == metric_name) %>%
    { if (metric_name == "RMSE") arrange(., Mean, Method) else arrange(., desc(Mean), Method) } %>%
    mutate(
      DisplayOrder = factor(seq_len(n()), levels = seq_len(n())),
      Method = factor(as.character(Method), levels = legend_order)
    )

  point_df <- data_long %>%
    filter(Setting == setting_name, Metric == metric_name) %>%
    mutate(Method = factor(as.character(Method), levels = legend_order)) %>%
    inner_join(plot_df %>% select(Method, DisplayOrder), by = "Method")

  lower_value <- min(c(point_df$Value, plot_df$Mean - plot_df$SEM), na.rm = TRUE)
  upper_value <- max(c(point_df$Value, plot_df$Mean + plot_df$SEM), na.rm = TRUE)
  lower_limit <- if (lower_value < 0) lower_value - 0.2 else 0
  upper_limit <- upper_value + 0.2
  tick_step <- if (metric_name == "RMSE") {
    if (setting_name == "Drug-blind") 1 else 0.5
  } else if (metric_name == "R2" && setting_name == "Drug-blind") {
    0.5
  } else {
    0.2
  }

  ggplot(plot_df, aes(DisplayOrder, Mean, fill = Method)) +
    geom_col(width = 0.9, color = NA, linewidth = 0) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                  width = 0.16, linewidth = 0.8, color = "black") +
    geom_jitter(data = point_df, aes(DisplayOrder, Value), inherit.aes = FALSE,
                width = 0.08, height = 0, size = 1.5, color = "#6B6B6B", alpha = 0.9) +
    labs(x = NULL, y = if (show_y_title) "Mean Value" else "", fill = "Model") +
    scale_fill_manual(values = model_palette, breaks = legend_order, drop = FALSE) +
    scale_y_continuous(
      limits = c(lower_limit, upper_limit),
      breaks = scales::breaks_width(tick_step),
      expand = c(0, 0)
    ) +
    facet_wrap(~Metric, scales = "fixed") +
    theme_fig2_panel()
}

plot_one_setting <- function(setting_name) {
  plots <- lapply(seq_along(metric_levels), function(i) {
    plot_metric_panel(setting_name, metric_levels[i], summary_df, df_long,
                      show_y_title = i %in% c(1, 4))
  })
  wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 8.5)) &
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}

for (setting_name in setting_levels) {
  p <- plot_one_setting(setting_name)
  stub <- str_replace_all(setting_name, "-", "_")
  ggsave(file.path(out_dir, paste0("Fig2_", stub, "_benchmark.pdf")), p,
         width = 12, height = 9, units = "in", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0("Fig2_", stub, "_benchmark.png")), p,
         width = 12, height = 9, units = "in", dpi = 600, limitsize = FALSE)
}
