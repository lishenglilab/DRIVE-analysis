############################################################
## fig4b_density_heatmap.R
##
## Purpose:
##   Generate the tissue-level IC50 density heatmap used for Fig. 4B.
##
## Input:
##   high_quality_data.csv
##
## Required columns:
##   1. OncotreeLineage : tissue or cancer lineage
##   2. cell_line       : cell-line name or ID
##   3. drug            : drug name or drug ID
##   4. ic50            : IC50 or log-transformed IC50 value
##
## Output:
##   1. Fig4B_Tissue_Summary_Density_Heatmap.pdf
##   2. Fig4B_Tissue_Summary_Density_Heatmap.png
##   3. Fig4B_tissue_summary_matrix.csv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
})

############################################################
## 1. Set working directory and input file
############################################################

## Example path. Replace this with your own project directory.
work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "high_quality_data.csv"

out_dir <- "fig4b_density_heatmap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Load input data
############################################################

high_quality_data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("OncotreeLineage", "cell_line", "drug", "ic50")
missing_cols <- setdiff(required_cols, colnames(high_quality_data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

high_quality_data <- high_quality_data %>%
  filter(
    !is.na(OncotreeLineage),
    !is.na(cell_line),
    !is.na(drug),
    !is.na(ic50)
  )

high_quality_data$ic50 <- as.numeric(high_quality_data$ic50)

############################################################
## 3. Summarize IC50 values by tissue and drug
############################################################

tissue_summary <- high_quality_data %>%
  group_by(OncotreeLineage, drug) %>%
  summarise(
    mean_ic50 = mean(ic50, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 4. Convert long-format table to drug-by-tissue matrix
############################################################

tissue_matrix <- tissue_summary %>%
  select(OncotreeLineage, drug, mean_ic50) %>%
  pivot_wider(
    names_from = OncotreeLineage,
    values_from = mean_ic50
  ) %>%
  column_to_rownames("drug") %>%
  as.matrix()

write.csv(
  tissue_matrix,
  file = file.path(out_dir, "Fig4B_tissue_summary_matrix.csv"),
  row.names = TRUE
)

############################################################
## 5. Plot tissue-level density heatmap
############################################################

plot_width <- max(12, ncol(tissue_matrix) * 0.5)
plot_height <- 10

pdf(
  file = file.path(out_dir, "Fig4B_Tissue_Summary_Density_Heatmap.pdf"),
  width = plot_width,
  height = plot_height
)

ht <- densityHeatmap(
  tissue_matrix,
  show_quantiles = TRUE,
  ylab = "IC50 Values"
)

draw(
  ht,
  column_title = "Density Heatmap by Tissue"
)

dev.off()

png(
  file = file.path(out_dir, "Fig4B_Tissue_Summary_Density_Heatmap.png"),
  width = plot_width * 100,
  height = plot_height * 100,
  res = 300
)

ht <- densityHeatmap(
  tissue_matrix,
  show_quantiles = TRUE,
  ylab = "IC50 Values"
)

draw(
  ht,
  column_title = "Density Heatmap by Tissue"
)

dev.off()

############################################################
## Finished
############################################################

message("Fig. 4B density heatmap finished.")
message("Output directory: ", normalizePath(out_dir))
message("Output files:")
message("  - Fig4B_Tissue_Summary_Density_Heatmap.pdf")
message("  - Fig4B_Tissue_Summary_Density_Heatmap.png")
message("  - Fig4B_tissue_summary_matrix.csv")
