############################################################
## supfig_s15_drug_response_comparison.R
##
## Purpose:
##   Generate Supplementary Fig. S15 comparing experimental and
##   reference IC50 values across candidate compounds.
##
## Input:
##   drug_response_ic50.csv
##
## Required columns:
##   1. CellLine
##   2. Drug_Group
##   3. Type
##   4. Compound_Name
##   5. IC50
##
## Output:
##   1. Sup3.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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

input_file <- "drug_response_ic50.csv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

## -------------------------------------------------------------------------
## 1. Read data
## Expected columns:
##   CellLine      : character (e.g., HCT116, RKO, Lovo)
##   Drug_Group    : character (grouping factor, e.g., compound names)
##   Type          : character ("Test" or "Ref")
##   Compound_Name : character (specific compound name)
##   IC50          : numeric (raw IC50 values)
## -------------------------------------------------------------------------
raw_data <- read.csv(input_file, stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 2. Data processing
## -------------------------------------------------------------------------
plot_data <- raw_data %>%
  mutate(LNIC50 = log(IC50)) %>%
  mutate(
    Fill_Label = ifelse(Type == "Test", Compound_Name, "Reference (5-Fu)"),
    Type = factor(Type, levels = c("Test", "Ref")),
    Drug_Group = factor(Drug_Group, levels = c(
      "Alisporivir", "Josamycin", "NIM811", "Rifabutin",
      "Cephalomannine", "Belotecan"
    ))
  )

## -------------------------------------------------------------------------
## 3. Color mapping
## -------------------------------------------------------------------------
custom_colors <- c(
  "Alisporivir"    = "#E74C3C",
  "Josamycin"      = "#1ABC9C",
  "NIM811"         = "#F39C12",
  "Rifabutin"      = "#9B59B6",
  "Cephalomannine" = "#D35400",
  "Belotecan"      = "#C2185B",
  "Reference (5-Fu)" = "#377EB8"
)

## -------------------------------------------------------------------------
## 4. Define plotting function
## -------------------------------------------------------------------------
create_cell_plot <- function(df, cell_name) {
  ggplot(df, aes(x = Drug_Group, y = LNIC50, fill = Fill_Label, group = Type)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      size = 0.2
    ) +
    scale_fill_manual(values = custom_colors, name = "Compound") +
    labs(
      title = paste(cell_name, "Cells"),
      y = "LN IC50 Value",
      x = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      legend.position = "right"
    ) +
    geom_hline(yintercept = 0, color = "black", size = 0.5)
}

## -------------------------------------------------------------------------
## 5. Generate plots for each cell line
## -------------------------------------------------------------------------
p1 <- create_cell_plot(filter(plot_data, CellLine == "HCT116"), "HCT116")
p2 <- create_cell_plot(filter(plot_data, CellLine == "RKO"), "RKO")
p3 <- create_cell_plot(filter(plot_data, CellLine == "Lovo"), "Lovo")

## -------------------------------------------------------------------------
## 6. Save output (arrange plots in one PDF)
## -------------------------------------------------------------------------
pdf("Sup3.pdf", width = 12, height = 10)
print(p1)
print(p2)
print(p3)
dev.off()

message("Supplementary Fig. S3 generated: Sup3.pdf")
