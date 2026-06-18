############################################################
## fig1e_pathway_drug_count_barplot.R
##
## Purpose:
##   Generate the Fig. 1E pathway-level drug count bar plot.
##
## Input:
##   figureE_panel.csv
##
## Required columns:
##   1. pathway    : pathway name shown on the y-axis
##   2. drug_type  : one of Chemotherapy, Hormonal, or Targeted
##   3. drug_count : number of drugs for the pathway
##
## Output:
##   1. figureE_panel.pdf
##   2. figureE_panel.svg
##   3. figureE_panel.png
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

input_file <- "figureE_panel.csv"

type_colors <- c(
  "Chemotherapy" = "#F4A3A8",
  "Hormonal" = "#95D095",
  "Targeted" = "#9BB8E8"
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

required_cols <- c("pathway", "drug_type", "drug_count")
df <- read.csv(input_file, stringsAsFactors = FALSE)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    pathway = factor(pathway, levels = rev(pathway)),
    drug_type = factor(drug_type, levels = c("Chemotherapy", "Hormonal", "Targeted"))
  )

p <- ggplot(df, aes(x = drug_count, y = pathway, fill = drug_type)) +
  geom_col(color = "black", linewidth = 0.7, width = 0.85) +
  geom_text(aes(label = drug_count, x = drug_count + 0.2), hjust = 0, size = 5) +
  scale_fill_manual(values = type_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of drugs", y = NULL, fill = NULL) +
  theme_classic(base_size = 16) +
  theme(
    legend.position.inside = c(0.80, 0.48),
    legend.position = "inside",
    legend.text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

ggsave("figureE_panel.pdf", p, width = 11, height = 9)
ggsave("figureE_panel.svg", p, width = 11, height = 9)
ggsave("figureE_panel.png", p, width = 11, height = 9, dpi = 600, bg = "white")
