############################################################
## fig1b_drug_encoding_tile_plot.R
##
## Purpose:
##   Generate the Fig. 1B feature-encoding tile plot.
##
## Input:
##   figureB_panel.csv
##
## Required columns:
##   1. method  : model name
##   2. feature : feature type shown on the x-axis
##   3. present : 1 indicates the fixed code-side color is used; 0 indicates gray
##
## Output:
##   1. figureB_panel.pdf
##   2. figureB_panel.svg
##   3. figureB_panel.png
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

input_file <- "figureB_panel.csv"

method_order <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP", "GPDRP", "NERD", "PaccMann", "Precily",
  "GADRP", "DeepCCDS"
)

feature_order <- c("Molecular graph", "Fingerprint", "ESPF", "ECFP", "Morgan", "SMILESVec")

feature_colors <- c(
  "Molecular graph" = "#F58B8F",
  "Fingerprint" = "#A45AED",
  "ESPF" = "#5BE5B1",
  "ECFP" = "#EC7CB3",
  "Morgan" = "#F8C867",
  "SMILESVec" = "#FBF65C"
)

background_color <- "#D9D9D9"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

required_cols <- c("method", "feature", "present")
df <- read.csv(input_file, stringsAsFactors = FALSE)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    method = factor(method, levels = rev(method_order)),
    feature = factor(feature, levels = feature_order),
    present = as.integer(present)
  )

grid_df <- expand_grid(
  method = factor(rev(method_order), levels = rev(method_order)),
  feature = factor(feature_order, levels = feature_order)
)

plot_df <- grid_df %>%
  left_join(df, by = c("method", "feature")) %>%
  mutate(
    present = replace_na(present, 0),
    fill_color = case_when(
      present == 1 ~ unname(feature_colors[as.character(feature)]),
      TRUE ~ background_color
    )
  )

p <- ggplot(plot_df, aes(x = feature, y = method)) +
  geom_tile(aes(fill = fill_color), color = "white", linewidth = 1.0) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

ggsave("figureB_panel.pdf", p, width = 6, height = 10)
ggsave("figureB_panel.svg", p, width = 6, height = 10)
ggsave("figureB_panel.png", p, width = 6, height = 10, dpi = 600, bg = "white")
