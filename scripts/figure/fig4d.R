############################################################
## fig4d_prediction_vs_observation_bubble_plot.R
##
## Purpose:
##   Generate the Fig. 4D tissue-level measured-vs-predicted
##   sensitivity bubble plot.
##
## Input:
##   1. ml_ensemble_predictions.csv
##   2. depmap.csv
##
## Output:
##   1. Fig4D_outputs/Fig4D_plot_data.csv
##   2. Fig4D_outputs/Fig4D_ratio_summary.csv
##   3. Fig4D_outputs/Fig4D_tissue_bubble_plot.pdf
##   4. Fig4D_outputs/Fig4D_tissue_bubble_plot.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
})

############################################################
## 1. User settings
############################################################

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

pred_file <- "ml_ensemble_predictions.csv"
true_file <- "depmap.csv"

out_dir <- "Fig4D_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sensitive_cutoff <- 0

############################################################
## 2. Read data
############################################################

if (!file.exists(pred_file)) {
  stop("Input file not found: ", pred_file)
}

if (!file.exists(true_file)) {
  stop("Input file not found: ", true_file)
}

pred_raw <- read_csv(pred_file, show_col_types = FALSE)

true_matrix <- read.csv(
  true_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################################################
## 3. Standardize prediction columns
############################################################

colnames(pred_raw) <- gsub("\\s+", "_", colnames(pred_raw))
colnames(pred_raw) <- gsub("^X\\.", "", colnames(pred_raw))
colnames(pred_raw) <- gsub("\\.+", "_", colnames(pred_raw))

cat("Prediction file columns:\n")
print(colnames(pred_raw))

if (all(c("level_0", "level_1", "Ensemble_Score") %in% colnames(pred_raw))) {
  
  pred_raw <- pred_raw %>%
    dplyr::rename(
      cell_line = dplyr::all_of("level_0"),
      drug_name = dplyr::all_of("level_1")
    )
  
} else if (all(c("cell_line", "drug_name", "Ensemble_Score") %in% colnames(pred_raw))) {
  
  pred_raw <- pred_raw
  
} else if (all(c("Cell_line", "Drug_name", "Ensemble_Score") %in% colnames(pred_raw))) {
  
  pred_raw <- pred_raw %>%
    dplyr::rename(
      cell_line = dplyr::all_of("Cell_line"),
      drug_name = dplyr::all_of("Drug_name")
    )
  
} else if (all(c("CellLine", "Drug", "Ensemble_Score") %in% colnames(pred_raw))) {
  
  pred_raw <- pred_raw %>%
    dplyr::rename(
      cell_line = dplyr::all_of("CellLine"),
      drug_name = dplyr::all_of("Drug")
    )
  
} else {
  
  stop(
    "Prediction file column names are not recognized.\n",
    "Required formats include one of:\n",
    "1) level_0, level_1, Ensemble_Score\n",
    "2) cell_line, drug_name, Ensemble_Score\n",
    "3) Cell_line, Drug_name, Ensemble_Score\n",
    "4) CellLine, Drug, Ensemble_Score\n",
    "Current columns are:\n",
    paste(colnames(pred_raw), collapse = ", ")
  )
}

pred_raw <- pred_raw %>%
  dplyr::mutate(
    cell_line = as.character(cell_line),
    drug_name = as.character(drug_name),
    Ensemble_Score = as.numeric(Ensemble_Score)
  ) %>%
  dplyr::filter(
    !is.na(cell_line),
    !is.na(drug_name),
    !is.na(Ensemble_Score)
  )

############################################################
## 4. Drug-name mapping
############################################################

drug_name_mapping <- c(
  "CPP" = "MCPP",
  "BF2.649" = "PITOLISANT",
  "L-798,106" = "L-798106",
  "MLN-8054" = "MLN8054",
  "ONC201" = "TIC10",
  "PLX-4720" = "PLX4720",
  "XL-647" = "XL647",
  "blebbistatin-(-)" = "BLEBBISTATIN-(+/-)",
  "cyclosporine" = "CYCLOSPORIN-A",
  "lorlatinib" = "PF-06463922"
)

clean_drug_name <- function(x) {
  x <- as.character(x)
  mapped <- drug_name_mapping[x]
  x <- ifelse(is.na(mapped), x, mapped)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

############################################################
## 5. Experimental drug-level statistics
############################################################

true_long <- true_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_line_true") %>%
  tidyr::pivot_longer(
    cols = -cell_line_true,
    names_to = "true_drug_name",
    values_to = "Actual_IC50"
  ) %>%
  dplyr::mutate(
    true_drug_name = as.character(true_drug_name),
    Actual_IC50 = as.numeric(Actual_IC50)
  ) %>%
  dplyr::filter(!is.na(Actual_IC50))

true_drug_stats_raw <- true_long %>%
  dplyr::group_by(true_drug_name) %>%
  dplyr::summarise(
    true_total = dplyr::n_distinct(cell_line_true),
    true_negative_count = dplyr::n_distinct(cell_line_true[Actual_IC50 < sensitive_cutoff]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    true_clean_name = clean_drug_name(true_drug_name)
  )

true_drug_stats <- true_drug_stats_raw %>%
  dplyr::group_by(true_clean_name) %>%
  dplyr::summarise(
    true_drug_name = paste(unique(true_drug_name), collapse = ";"),
    true_total = sum(true_total, na.rm = TRUE),
    true_negative_count = sum(true_negative_count, na.rm = TRUE),
    true_negative_ratio = true_negative_count / true_total,
    .groups = "drop"
  )

############################################################
## 6. Predicted drug-level statistics
############################################################

pred_drug_stats_raw <- pred_raw %>%
  dplyr::group_by(drug_name) %>%
  dplyr::summarise(
    pred_total = dplyr::n_distinct(cell_line),
    pred_negative_count = dplyr::n_distinct(cell_line[Ensemble_Score < sensitive_cutoff]),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pred_clean_name = clean_drug_name(drug_name)
  )

pred_drug_stats <- pred_drug_stats_raw %>%
  dplyr::group_by(pred_clean_name) %>%
  dplyr::summarise(
    drug_name = paste(unique(drug_name), collapse = ";"),
    pred_total = sum(pred_total, na.rm = TRUE),
    pred_negative_count = sum(pred_negative_count, na.rm = TRUE),
    pred_negative_ratio = pred_negative_count / pred_total,
    .groups = "drop"
  )

############################################################
## 7. Merge true and predicted statistics
############################################################

combined_stats_final <- pred_drug_stats %>%
  dplyr::inner_join(
    true_drug_stats,
    by = c("pred_clean_name" = "true_clean_name")
  ) %>%
  dplyr::mutate(
    count_ratio = true_total / pred_total
  ) %>%
  dplyr::select(
    drug_name,
    true_drug_name,
    true_negative_ratio,
    pred_negative_ratio,
    true_total,
    pred_total,
    count_ratio,
    true_negative_count,
    pred_negative_count
  ) %>%
  dplyr::filter(
    !is.na(true_negative_ratio),
    !is.na(pred_negative_ratio),
    !is.na(count_ratio),
    pred_total > 0
  ) %>%
  dplyr::arrange(desc(count_ratio))

write_csv(
  combined_stats_final,
  file.path(out_dir, "Fig4D_combined_stats_final.csv")
)

############################################################
## 8. Save unmatched predicted drugs
############################################################

unmatched_predicted_drugs <- pred_drug_stats %>%
  dplyr::anti_join(
    true_drug_stats,
    by = c("pred_clean_name" = "true_clean_name")
  ) %>%
  dplyr::arrange(drug_name)

write_csv(
  unmatched_predicted_drugs,
  file.path(out_dir, "Fig4D_unmatched_predicted_drugs.csv")
)

cat("Predicted drugs:", nrow(pred_drug_stats), "\n")
cat("Experimental drugs:", nrow(true_drug_stats), "\n")
cat("Matched drugs for Fig.4D:", nrow(combined_stats_final), "\n")
cat("Unmatched predicted drugs:", nrow(unmatched_predicted_drugs), "\n")

cat("\ntrue_negative_ratio summary:\n")
print(summary(combined_stats_final$true_negative_ratio))

cat("\npred_negative_ratio summary:\n")
print(summary(combined_stats_final$pred_negative_ratio))

cat("\ncount_ratio summary:\n")
print(summary(combined_stats_final$count_ratio))

############################################################
## 9. Prepare plot data
############################################################

plot_df <- combined_stats_final %>%
  dplyr::mutate(
    true_negative_ratio = as.numeric(true_negative_ratio),
    pred_negative_ratio = as.numeric(pred_negative_ratio),
    count_ratio = as.numeric(count_ratio),
    
    ## Clip only for visual scaling.
    ## Values > 0.3 are displayed with the same maximum size as 0.3.
    ## The original count_ratio is still saved in combined_stats_final.
    count_ratio_plot = pmin(pmax(count_ratio, 0), 0.3)
  ) %>%
  dplyr::filter(
    !is.na(true_negative_ratio),
    !is.na(pred_negative_ratio),
    !is.na(count_ratio_plot)
  )


############################################################
## 10. Plot Fig.4D
############################################################

p <- ggplot(
  plot_df,
  aes(
    x = true_negative_ratio,
    y = pred_negative_ratio
  )
) +
  geom_point(
    aes(
      size = count_ratio_plot,
      fill = count_ratio_plot
    ),
    shape = 21,
    color = "white",
    stroke = 0.25,
    alpha = 0.82
  ) +
  
  ## size: continuous scale; legend shown as hollow circles
  scale_size_continuous(
    name = "Count ratio",
    range = c(0.15, 11),
    limits = c(0, 0.3),
    breaks = c(0.1, 0.2, 0.3),
    labels = c("0.1", "0.2", "0.3")
  ) +
  
  ## Fill: gradient for color bar
  scale_fill_gradient(
    name = "Count ratio",
    low = "#FFB3B3",
    high = "red",
    limits = c(0, 0.3),
    breaks = c(0.1, 0.2, 0.3),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(4.0, "cm"),
      barheight = unit(0.35, "cm")
    )
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  labs(
    x = "Ratio of cell lines with experimental data",
    y = "Ratio of cell lines with DRIVE prediction"
  ) +
  
  guides(
    size = guide_legend(
      title = "Count ratio",
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      nrow = 1,
      order = 1,
      override.aes = list(
        shape = 21,
        fill = NA,          # hollow circles in size legend
        color = "black",
        alpha = 1,
        stroke = 0.4
      )
    ),
    fill = guide_colorbar(
      title = "Count ratio",
      title.position = "top",
      title.hjust = 0.5,
      order = 2,
      barwidth = grid::unit(3.8, "cm"),
      barheight = grid::unit(0.32, "cm")
    )
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.65),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    
    ## legends in one row at top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    
    plot.margin = margin(16, 16, 10, 12)
  )

print(p)

ggsave(
  filename = file.path(out_dir, "Fig4D_prediction_validation_scatter.pdf"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Fig4D_prediction_validation_scatter.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)
