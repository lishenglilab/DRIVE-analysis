############################################################
## fig4f_tissue_class_checkerboard.R
##
## Purpose:
##   Generate the separated Fig. 4F top stacked bar, bottom
##   checkerboard bubble matrix, and combined layout.
##
## Input:
##   plot_data_for_R.csv
##
## Required columns:
##   1. tissue  : tissue label
##   2. pathway : compound class label
##   3. count   : number of compounds in the tissue-class pair
##
## Output:
##   1. Fig4F_top_stacked_bar.pdf
##   2. Fig4F_top_stacked_bar.png
##   3. Fig4F_bottom_bubble_matrix.pdf
##   4. Fig4F_bottom_bubble_matrix.png
##   5. Fig4F_combined.pdf
##   6. Fig4F_combined.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
  library(gridExtra)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

extract_legend <- function(plot_obj) {
  grob_obj <- ggplotGrob(plot_obj)
  guide_idx <- which(vapply(grob_obj$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_idx) == 0) {
    return(NULL)
  }
  grob_obj$grobs[[guide_idx[1]]]
}

save_grob <- function(grob_obj, filename, width, height, dpi = 300) {
  ext <- tools::file_ext(filename)
  if (tolower(ext) == "pdf") {
    pdf(filename, width = width, height = height, useDingbats = FALSE)
    grid::grid.draw(grob_obj)
    dev.off()
  } else {
    png(filename, width = width, height = height, units = "in", res = dpi, bg = "white")
    grid::grid.draw(grob_obj)
    dev.off()
  }
}

work_dir <- get_script_dir()
setwd(work_dir)

input_file <- "plot_data_for_R.csv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

final_plot_data <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("tissue", "pathway", "count")
missing_cols <- setdiff(required_cols, colnames(final_plot_data))

if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

tissue_order_full <- c(
  "Biliary tract", "Bladder", "Bone", "Brain", "Breast", "Cervix", "Esophagus",
  "Eye", "Haematopoietic and lymphoid", "Head and neck", "Kidney",
  "Large intestine", "Liver", "Lung", "Ovary", "Pancreas",
  "Peripheral nervous system", "Prostate", "Skin", "Small intestine",
  "Soft tissue", "Stomach", "Testis", "Thyroid", "Uterus", "Vulva"
)

stack_order <- c(
  "Terpenoids",
  "Shikimates and Phenylpropanoids",
  "Polyketides",
  "Fatty acids",
  "Carbohydrates",
  "Amino acids and Peptides",
  "Alkaloids"
)

matrix_order <- c(
  "Alkaloids",
  "Amino acids and Peptides",
  "Carbohydrates",
  "Fatty acids",
  "Polyketides",
  "Shikimates and Phenylpropanoids",
  "Terpenoids"
)

palette_values <- c(
  "Alkaloids" = "#C5B7E3",
  "Amino acids and Peptides" = "#C56D1E",
  "Carbohydrates" = "#6F6F6F",
  "Fatty acids" = "#FFF07C",
  "Polyketides" = "#4B79B8",
  "Shikimates and Phenylpropanoids" = "#F4BE86",
  "Terpenoids" = "#7FC57A"
)

plot_data <- final_plot_data %>%
  mutate(
    tissue = str_trim(as.character(tissue)),
    pathway = str_trim(as.character(pathway)),
    count = as.numeric(count),
    bar_weight = if ("bar_weight" %in% colnames(final_plot_data)) as.numeric(bar_weight) else as.numeric(count)
  ) %>%
  filter(
    !is.na(tissue), tissue != "",
    !is.na(pathway), pathway != "",
    !is.na(count), count > 0,
    !is.na(bar_weight), bar_weight > 0,
    pathway %in% stack_order
  ) %>%
  group_by(tissue, pathway) %>%
  summarise(
    count = sum(count),
    bar_weight = sum(bar_weight),
    .groups = "drop"
  )

tissue_levels <- tissue_order_full[tissue_order_full %in% unique(plot_data$tissue)]
if (length(tissue_levels) == 0) {
  tissue_levels <- unique(plot_data$tissue)
}

plot_data <- plot_data %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    pathway_stack = factor(pathway, levels = stack_order),
    pathway_matrix = factor(pathway, levels = rev(matrix_order)),
    log2_count = log2(count)
  ) %>%
  group_by(tissue) %>%
  mutate(
    total_bar_weight = sum(bar_weight),
    proportion = bar_weight / total_bar_weight
  ) %>%
  ungroup()

background_df <- expand.grid(
  tissue = tissue_levels,
  pathway_matrix = rev(matrix_order),
  stringsAsFactors = FALSE
) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    pathway_matrix = factor(pathway_matrix, levels = rev(matrix_order)),
    tissue_id = as.integer(tissue),
    pathway_id = as.integer(pathway_matrix),
    checker_fill = ifelse((tissue_id + pathway_id) %% 2 == 0, "#F5F5F5", "#FFFFFF")
  )

top_plot <- ggplot(plot_data, aes(x = tissue, y = proportion, fill = pathway_stack)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  scale_fill_manual(
    values = palette_values,
    breaks = matrix_order,
    name = "Compound Class"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", size = 15),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 10.5),
    plot.margin = margin(8, 8, 4, 8)
  )

bottom_plot <- ggplot(plot_data, aes(x = tissue, y = pathway_matrix)) +
  geom_tile(
    data = background_df,
    aes(x = tissue, y = pathway_matrix),
    fill = background_df$checker_fill,
    color = "#CFCFCF",
    linewidth = 0.45,
    width = 1,
    height = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(fill = log2_count),
    shape = 21,
    size = 7.2,
    stroke = 0.45,
    color = "black"
  ) +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#DCEEFF", "#AACCF2", "#74AAE8", "#4C6ACB", "#0909A9"),
    limits = c(9, 19),
    breaks = c(9, 11, 13, 15, 17, 19),
    name = "Log2(Count)",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(4.2, "cm"),
      barheight = unit(0.28, "cm")
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 10.5),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 10.5),
    plot.margin = margin(4, 8, 4, 8)
  )

legend_class <- extract_legend(
  top_plot +
    theme(
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

legend_count <- extract_legend(
  bottom_plot +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

top_panel <- ggplotGrob(top_plot + theme(legend.position = "none"))
bottom_panel <- ggplotGrob(bottom_plot + theme(legend.position = "none"))

main_panel <- arrangeGrob(
  top_panel,
  bottom_panel,
  ncol = 1,
  heights = c(0.88, 1.34)
)

legend_row <- arrangeGrob(
  legend_count,
  legend_class,
  ncol = 2,
  widths = c(0.43, 0.57)
)

combined_grob <- arrangeGrob(
  main_panel,
  legend_row,
  ncol = 1,
  heights = c(5.1, 1.9)
)

ggsave("Fig4F_top_stacked_bar.pdf", top_plot, width = 13.0, height = 3.1, dpi = 300)
ggsave("Fig4F_top_stacked_bar.png", top_plot, width = 13.0, height = 3.1, dpi = 600, bg = "white")

ggsave("Fig4F_bottom_bubble_matrix.pdf", bottom_plot, width = 13.0, height = 4.8, dpi = 300)
ggsave("Fig4F_bottom_bubble_matrix.png", bottom_plot, width = 13.0, height = 4.8, dpi = 600, bg = "white")

save_grob(combined_grob, "Fig4F_combined.pdf", width = 14.2, height = 9.6)
save_grob(combined_grob, "Fig4F_combined.png", width = 14.2, height = 9.6, dpi = 600)

message("Fig. 4F top panel, bottom panel, and combined figure finished.")
