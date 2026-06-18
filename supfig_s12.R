############################################################
## supfig_s12_predicted_ic50_violin_by_tissue.R
##
## Purpose:
##   Generate Supplementary Fig. S12 showing the distribution
##   of predicted IC50 values across tissues.
##
## Input:
##   predicted_ic50_by_tissue.csv
##
## Required columns:
##   1. Tissue
##   2. IC50
##
## Output:
##   1. Predicted_IC50_by_Tissue.pdf
##   2. Predicted_IC50_by_Tissue.png
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
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

############################################################
## 1. Load data
############################################################

input_file <- "predicted_ic50_by_tissue.csv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

ic50_df <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("Tissue", "IC50")
missing_cols <- setdiff(required_cols, colnames(ic50_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

############################################################
## 2. Data preparation
############################################################

ic50_df <- ic50_df %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = unique(Tissue)
    )
  )

############################################################
## 3. Visualization
############################################################

p <- ggplot(
  ic50_df,
  aes(
    x = Tissue,
    y = IC50
  )
) +
  geom_violin(
    fill = "#73A6D1",
    color = "black",
    linewidth = 0.5,
    trim = FALSE,
    scale = "width"
  ) +
  labs(
    x = NULL,
    y = expression(
      "Predicted IC"[50] ~ "values"
    )
  ) +
  scale_y_continuous(
    limits = c(-7, 9),
    breaks = c(-5, 0, 5)
  ) +
  theme(
    panel.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 1
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 12
    ),
    
    axis.title.y = element_text(
      size = 16
    ),
    
    axis.title.x = element_blank(),
    
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    )
  )

p

ggsave(
  filename = "Predicted_IC50_by_Tissue.pdf",
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = "Predicted_IC50_by_Tissue.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
