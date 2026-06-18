############################################################
## fig1c_model_algorithm_bubble_plot.R
##
## Purpose:
##   Generate the Fig. 1C model-by-algorithm bubble matrix.
##
## Input:
##   figureC_panel.csv
##
## Required columns:
##   1. method     : model name
##   2. algorithm  : algorithm label shown on the x-axis
##   3. present    : 1 indicates a bubble should be drawn; 0 indicates empty
##
## Optional columns:
##   1. fill_color : fill color for bubbles when present == 1
##
## Output:
##   1. figureC_panel.pdf
##   2. figureC_panel.svg
##   3. figureC_panel.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

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

input_file <- "figureC_panel.csv"

method_order <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP", "GPDRP", "NERD", "PaccMann", "Precily",
  "GADRP", "DeepCCDS"
)

algo_order <- c("DNN", "GCN", "AE", "Transformer", "MCA", "KNN", "MLP", "GAT", "GIN", "GAT_GCN")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

required_cols <- c("method", "algorithm", "present")
df <- read.csv(input_file, stringsAsFactors = FALSE)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    method = factor(method, levels = rev(method_order)),
    algorithm = factor(algorithm, levels = algo_order)
  )

grid_df <- expand_grid(
  method = factor(rev(method_order), levels = rev(method_order)),
  algorithm = factor(algo_order, levels = algo_order)
)

plot_df <- grid_df %>%
  left_join(df, by = c("method", "algorithm")) %>%
  mutate(
    present = replace_na(present, 0),
    fill_color = replace_na(fill_color, "#FFFFFF")
  )

p <- ggplot(plot_df, aes(x = algorithm, y = method)) +
  geom_tile(fill = "white", color = "#D9D9D9", linewidth = 1.0) +
  geom_point(
    data = plot_df %>% filter(present == 1),
    aes(fill = fill_color),
    shape = 21,
    size = 11,
    stroke = 1.0,
    color = "white"
  ) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "#D9D9D9", fill = NA, linewidth = 1.2)
  )

ggsave("figureC_panel.pdf", p, width = 8, height = 10)
ggsave("figureC_panel.svg", p, width = 8, height = 10)
ggsave("figureC_panel.png", p, width = 8, height = 10, dpi = 600, bg = "white")
