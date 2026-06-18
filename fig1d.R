############################################################
## fig1d_tissue_subtype_polar_barplot.R
##
## Purpose:
##   Generate the Fig. 1D tissue-subtype polar bar plot.
##
## Input:
##   figureD_panel.csv
##
## Required columns:
##   1. system_group : major tissue/system group
##   2. subtype      : subtype label displayed inside each bar
##   3. count        : subtype count used as bar height
##   4. fill_color   : fill color for the subtype bar
##
## Output:
##   1. figureD_panel.pdf
##   2. figureD_panel.svg
##   3. figureD_panel.png
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

input_file <- "figureD_panel.csv"

group_order <- c(
  "Lung", "Blood", "Urogenital system", "Digestive system", "Aero digestive tract",
  "Nervous system", "Breast", "Skin", "Pancreas", "Kidney", "Bone",
  "Soft tissue", "Thyroid", "Large intestine"
)

outer_colors <- c(
  "Lung" = "#A7CEE5",
  "Blood" = "#3D8CC6",
  "Urogenital system" = "#7CC495",
  "Digestive system" = "#64B84A",
  "Aero digestive tract" = "#86A850",
  "Nervous system" = "#F28A8A",
  "Breast" = "#E61E1E",
  "Skin" = "#F2B26E",
  "Pancreas" = "#FF8418",
  "Kidney" = "#D99979",
  "Bone" = "#9A82D0",
  "Soft tissue" = "#9A7AA7",
  "Thyroid" = "#EED96B",
  "Large intestine" = "#B85E28"
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

required_cols <- c("system_group", "subtype", "count", "fill_color")
df <- read.csv(input_file, stringsAsFactors = FALSE)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(system_group = factor(system_group, levels = group_order)) %>%
  group_by(system_group) %>%
  mutate(group_total = n()) %>%
  ungroup() %>%
  arrange(system_group)

df$bar_id <- seq_len(nrow(df))
inner_radius <- 48
df$bar_height <- df$count
df$bar_top <- inner_radius + df$bar_height
df$label_y <- inner_radius + pmax(df$bar_height * 0.72, df$bar_height - 3)
max_count <- max(df$bar_height)
ring_breaks <- seq(10, ceiling(max_count / 10) * 10, by = 10)

group_ranges <- df %>%
  group_by(system_group) %>%
  summarise(
    start = min(bar_id) - 0.5,
    end = max(bar_id) + 0.5,
    mid = mean(bar_id),
    .groups = "drop"
  ) %>%
  mutate(
    outer_radius = inner_radius + max(df$bar_height) + 8,
    label_radius = inner_radius + max(df$bar_height) + 12,
    outer_color = outer_colors[as.character(system_group)]
  )

p <- ggplot(df, aes(x = factor(bar_id), y = bar_height, fill = fill_color)) +
  geom_hline(
    yintercept = ring_breaks,
    linetype = "dashed",
    color = "#C9C9C9",
    linewidth = 0.6
  ) +
  geom_col(color = "black", linewidth = 0.7, width = 0.95) +
  ylim(-inner_radius, inner_radius + max(df$bar_height) + 14) +
  geom_text(aes(label = count, y = label_y), size = 4.3) +
  geom_text(
    aes(label = subtype, y = inner_radius + pmin(pmax(bar_height * 0.24, 2.8), 8)),
    angle = 90,
    size = 4
  ) +
  scale_fill_identity() +
  coord_polar(clip = "off") +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    plot.margin = margin(18, 18, 18, 18)
  )

for (i in seq_len(nrow(group_ranges))) {
  start_id <- group_ranges$start[i]
  end_id <- group_ranges$end[i]
  mid_id <- group_ranges$mid[i]
  group_name <- as.character(group_ranges$system_group[i])
  p <- p +
    annotate(
      "rect",
      xmin = start_id,
      xmax = end_id,
      ymin = inner_radius + max(df$bar_height) + 6.5,
      ymax = inner_radius + max(df$bar_height) + 7.8,
      fill = group_ranges$outer_color[i],
      color = NA
    ) +
    annotate(
      "text",
      x = mid_id,
      y = inner_radius + max(df$bar_height) - 2,
      label = group_name,
      angle = 90 - 360 * (mid_id / nrow(df)),
      size = 4.4
    )
}

ggsave("figureD_panel.pdf", p, width = 10, height = 10)
ggsave("figureD_panel.svg", p, width = 10, height = 10)
ggsave("figureD_panel.png", p, width = 10, height = 10, dpi = 600, bg = "white")
