############################################################
## fig5a_crc_ic50_distribution_panels.R
##
## Purpose:
##   Generate the separated Fig. 5A top stacked-percentage panel
##   and bottom IC50 gradient panel.
##
## Input:
##   1. figure5_top_panel.csv
##   2. figure5_bottom_panel.csv
##
## Required columns for figure5_top_panel.csv:
##   1. cell_line
##   2. ic50_lt_neg1
##   3. ic50_neg1_to_0
##   4. ic50_0_to_1
##   5. ic50_gt_1
##
## Required columns for figure5_bottom_panel.csv:
##   1. cell_line
##   2. rank_index
##   3. ic50_value
##
## Output:
##   1. figure5_top_panel.pdf
##   2. figure5_top_panel.svg
##   3. figure5_top_panel.png
##   4. figure5_bottom_panel.pdf
##   5. figure5_bottom_panel.svg
##   6. figure5_bottom_panel.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
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

top_panel_file <- "figure5_top_panel.csv"
bottom_panel_file <- "figure5_bottom_panel.csv"

if (!file.exists(top_panel_file)) {
  stop("Input file not found: ", top_panel_file)
}

if (!file.exists(bottom_panel_file)) {
  stop("Input file not found: ", bottom_panel_file)
}

top_df <- read.csv(top_panel_file, stringsAsFactors = FALSE)
bottom_df <- read.csv(bottom_panel_file, stringsAsFactors = FALSE)

required_top <- c(
  "cell_line",
  "ic50_lt_neg1",
  "ic50_neg1_to_0",
  "ic50_0_to_1",
  "ic50_gt_1"
)

required_bottom <- c("cell_line", "rank_index", "ic50_value")

missing_top <- setdiff(required_top, colnames(top_df))
missing_bottom <- setdiff(required_bottom, colnames(bottom_df))

if (length(missing_top) > 0) {
  stop("Missing top-panel columns: ", paste(missing_top, collapse = ", "))
}

if (length(missing_bottom) > 0) {
  stop("Missing bottom-panel columns: ", paste(missing_bottom, collapse = ", "))
}

top_df <- top_df %>%
  distinct(cell_line, .keep_all = TRUE)

# The cell-line order in the top-panel file is treated as the final display order.
cell_order <- unique(top_df$cell_line)
shared_cell_order <- intersect(cell_order, unique(bottom_df$cell_line))

top_df <- top_df %>%
  filter(cell_line %in% shared_cell_order)

cell_order <- shared_cell_order

top_long <- top_df %>%
  pivot_longer(
    cols = -cell_line,
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(
    cell_line = factor(cell_line, levels = cell_order),
    category = factor(
      category,
      levels = c(
        "ic50_gt_1",
        "ic50_0_to_1",
        "ic50_neg1_to_0",
        "ic50_lt_neg1"
      ),
      labels = c(
        "IC50 > 1",
        "0 < IC50 \u2264 1",
        "-1 \u2264 IC50 < 0",
        "IC50 < -1"
      )
    )
  )

bottom_df <- bottom_df %>%
  filter(cell_line %in% cell_order) %>%
  mutate(
    cell_line = factor(cell_line, levels = cell_order),
    rank_index = as.numeric(rank_index),
    ic50_value = as.numeric(ic50_value)
  ) %>%
  arrange(cell_line, rank_index)

make_bottom_raster_df <- function(bottom_df, cell_order, n_rows = 1200, strip_px = 26) {
  dense_y <- seq(0, 1, length.out = n_rows)
  strip_list <- vector("list", length(cell_order))

  for (i in seq_along(cell_order)) {
    cell_name <- cell_order[i]
    cell_dat <- bottom_df %>%
      filter(cell_line == cell_name) %>%
      arrange(rank_index)

    x <- scales::rescale(cell_dat$rank_index, to = c(0, 1))
    y <- cell_dat$ic50_value
    dense_values <- approx(x = x, y = y, xout = dense_y, rule = 2)$y
    x_positions <- seq(i - 0.39, i + 0.39, length.out = strip_px)

    strip_list[[i]] <- expand_grid(
      x = x_positions,
      y = dense_y
    ) %>%
      mutate(ic50_value = rep(dense_values, times = strip_px))
  }

  bind_rows(strip_list)
}

bottom_raster_df <- make_bottom_raster_df(bottom_df, cell_order)

top_colors <- c(
  "IC50 < -1" = "#0B0E8A",
  "-1 \u2264 IC50 < 0" = "#6C8CD5",
  "0 < IC50 \u2264 1" = "#FFAA8A",
  "IC50 > 1" = "#980000"
)

heat_colors <- c("#8B78D0", "#F7F7F7", "#FFB299", "#FF5A36")

p_top <- ggplot(top_long, aes(x = cell_line, y = percentage, fill = category)) +
  geom_col(width = 0.78, color = NA) +
  scale_fill_manual(values = top_colors) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Percentage (%)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line.x = element_blank(),
    plot.margin = margin(10, 10, 2, 45)
  )

p_bottom <- ggplot() +
  geom_raster(
    data = bottom_raster_df,
    aes(x = x, y = y, fill = ic50_value),
    interpolate = TRUE
  ) +
  scale_fill_gradientn(
    colors = heat_colors,
    values = rescale(c(-5, 0, 4, 8)),
    limits = c(-5, 8),
    breaks = c(-5, -2.5, 0, 2.5, 5, 7.5),
    name = expression(IC[50] ~ "value")
  ) +
  scale_x_continuous(
    breaks = seq_along(cell_order),
    labels = cell_order,
    expand = c(0, 0)
  ) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = expression(IC[50] ~ "value")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position.inside = c(0.01, 0.52),
    legend.justification = c(0, 0.5),
    legend.direction = "vertical",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.key.height = unit(2.2, "cm"),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "black",
      size = 11
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line = element_blank(),
    plot.margin = margin(2, 10, 10, 45)
  ) +
  guides(
    fill = guide_colorbar(
      frame.colour = NA,
      ticks.colour = "black",
      barwidth = unit(0.5, "cm"),
      barheight = unit(4.0, "cm")
    )
  )

ggsave(
  "figure5_top_panel.pdf",
  p_top,
  width = 12,
  height = 2.8,
  units = "in"
)

ggsave(
  "figure5_top_panel.svg",
  p_top,
  width = 12,
  height = 2.8,
  units = "in"
)

ggsave(
  "figure5_top_panel.png",
  p_top,
  width = 12,
  height = 2.8,
  units = "in",
  dpi = 600
)

ggsave(
  "figure5_bottom_panel.pdf",
  p_bottom,
  width = 12,
  height = 6.2,
  units = "in"
)

ggsave(
  "figure5_bottom_panel.svg",
  p_bottom,
  width = 12,
  height = 6.2,
  units = "in"
)

ggsave(
  "figure5_bottom_panel.png",
  p_bottom,
  width = 12,
  height = 6.2,
  units = "in",
  dpi = 600
)
