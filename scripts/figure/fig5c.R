############################################################
## fig5c_crc_tcga_independent_top30_lollipop.R
##
## Purpose:
##   Generate Fig. 5C, a mirror lollipop plot comparing the top
##   frequently sensitive drugs between CRC cell lines and TCGA
##   patient samples.
##
## Core logic:
##   1. Read prediction files from two folders:
##        - full_prediction_results_cellline/
##        - full_prediction_results_tcga/
##   2. For each group independently:
##        - keep Ensemble_Score < -1
##        - count sensitive occurrences for each drug
##        - calculate frequency = count / total number of samples
##        - calculate mean sensitivity = mean(abs(Ensemble_Score))
##        - select top 30 drugs by frequency
##   3. Merge the two independent top-30 lists by rank.
##   4. Plot a mirror lollipop chart:
##        - left: CRC cell lines
##        - right: TCGA patients
##
## Input:
##   Folder 1:
##     full_prediction_results_cellline/*_full_predictions.csv
##
##   Folder 2:
##     full_prediction_results_tcga/*_full_predictions.csv
##
## Required columns in each prediction file:
##   - drug_name
##   - Ensemble_Score
##
## Output:
##   - Fig5C_independent_top30_crc_vs_tcga.csv
##   - Fig5C_independent_top30_crc_vs_tcga.pdf
##   - Fig5C_independent_top30_crc_vs_tcga.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

############################################################
## 1. Set working directory and parameters
############################################################

## Example path. Replace this with your own project directory.
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

cellline_prediction_dir <- "full_prediction_results_cellline"
tcga_prediction_dir <- "full_prediction_results_tcga"

out_dir <- "fig5c_crc_tcga_lollipop"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Sensitive response threshold.
## The final Fig. 5C used Ensemble_Score < -1.
sensitive_threshold <- -1

top_n_drugs <- 30

############################################################
## 2. Function: read and summarize one prediction folder
############################################################

process_drug_folder <- function(folder_path,
                                threshold = -1,
                                top_n = 30) {

  files <- list.files(
    path = folder_path,
    pattern = "_full_predictions\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No *_full_predictions.csv files found in folder: ", folder_path)
  }

  total_samples <- length(files)
  all_data_list <- list()

  for (i in seq_along(files)) {

    df <- tryCatch({
      ## First try tab-delimited format.
      d <- read.table(
        files[i],
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "\""
      )

      ## If tab-delimited reading fails to split columns, try comma-delimited format.
      if (ncol(d) < 2) {
        d <- read.table(
          files[i],
          header = TRUE,
          sep = ",",
          check.names = FALSE,
          stringsAsFactors = FALSE,
          quote = "\""
        )
      }

      d
    }, error = function(e) {
      warning("Failed to read file: ", files[i])
      return(NULL)
    })

    if (is.null(df)) {
      next
    }

    colnames(df) <- trimws(colnames(df))

    idx_name <- grep("drug_name", colnames(df), ignore.case = TRUE)
    idx_score <- grep("Ensemble_Score", colnames(df), ignore.case = TRUE)

    if (length(idx_name) > 0 && length(idx_score) > 0) {
      temp_df <- df[, c(idx_name[1], idx_score[1])]
      colnames(temp_df) <- c("drug_name", "Ensemble_Score")
      all_data_list[[length(all_data_list) + 1]] <- temp_df
    } else {
      warning("Required columns were not found in file: ", files[i])
    }
  }

  if (length(all_data_list) == 0) {
    stop("No valid prediction files were loaded from folder: ", folder_path)
  }

  full_df <- bind_rows(all_data_list)

  result <- full_df %>%
    mutate(
      drug_name = as.character(drug_name),
      Ensemble_Score = as.numeric(Ensemble_Score)
    ) %>%
    filter(
      !is.na(drug_name),
      drug_name != "",
      !is.na(Ensemble_Score),
      Ensemble_Score < threshold
    ) %>%
    group_by(drug_name) %>%
    summarise(
      count = n(),
      freq = count / total_samples,
      mean_sens = mean(abs(Ensemble_Score), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(freq), desc(mean_sens)) %>%
    slice_head(n = top_n) %>%
    mutate(rank = row_number())

  return(result)
}

############################################################
## 3. Process CRC cell-line and TCGA prediction folders
############################################################

df_cl <- process_drug_folder(
  folder_path = cellline_prediction_dir,
  threshold = sensitive_threshold,
  top_n = top_n_drugs
)

df_tcga <- process_drug_folder(
  folder_path = tcga_prediction_dir,
  threshold = sensitive_threshold,
  top_n = top_n_drugs
)

############################################################
## 4. Merge independent top-30 lists by rank
############################################################

plot_df <- full_join(
  df_cl %>%
    select(
      rank,
      drug_cl = drug_name,
      count_cl = count,
      freq_cl = freq,
      sens_cl = mean_sens
    ),
  df_tcga %>%
    select(
      rank,
      drug_tcga = drug_name,
      count_tcga = count,
      freq_tcga = freq,
      sens_tcga = mean_sens
    ),
  by = "rank"
) %>%
  arrange(rank)

write.csv(
  plot_df,
  file = file.path(out_dir, "Fig5C_independent_top30_crc_vs_tcga.csv"),
  row.names = FALSE
)

############################################################
## 5. Plot mirror lollipop chart
############################################################

center_gap <- 0.6
x_left_0 <- -center_gap
x_right_0 <- center_gap

p <- ggplot(plot_df) +

  ## Background center dashed lines
  geom_vline(
    xintercept = x_left_0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = x_right_0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +

  ## Left panel: CRC cell lines
  geom_segment(
    aes(x = x_left_0,
        xend = x_left_0 - freq_cl,
        y = -rank,
        yend = -rank),
    color = "grey85",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    aes(x = x_left_0 - freq_cl,
        y = -rank,
        size = sens_cl,
        color = "CRC cell lines"),
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_text(
    aes(x = x_left_0 + 0.03,
        y = -rank,
        label = drug_cl),
    hjust = 0,
    size = 3,
    fontface = "bold",
    na.rm = TRUE
  ) +

  ## Right panel: TCGA patients
  geom_segment(
    aes(x = x_right_0,
        xend = x_right_0 + freq_tcga,
        y = -rank,
        yend = -rank),
    color = "grey85",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    aes(x = x_right_0 + freq_tcga,
        y = -rank,
        size = sens_tcga,
        color = "TCGA patients"),
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_text(
    aes(x = x_right_0 - 0.03,
        y = -rank,
        label = drug_tcga),
    hjust = 1,
    size = 3,
    fontface = "bold",
    na.rm = TRUE
  ) +

  ## X-axis
  scale_x_continuous(
    limits = c(-1.8, 1.8),
    breaks = c(
      x_left_0 - 1.0,
      x_left_0 - 0.5,
      x_left_0,
      x_right_0,
      x_right_0 + 0.5,
      x_right_0 + 1.0
    ),
    labels = c("1.0", "0.5", "0.0", "0.0", "0.5", "1.0")
  ) +

  ## Colors and sizes
  scale_color_manual(
    values = c(
      "CRC cell lines" = "#FF8C00",
      "TCGA patients" = "#87CEEB"
    ),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(2.0, 2.5, 3.0, 3.5),
    name = "Mean Sensitivity (|Score|)"
  ) +

  ## Labels
  labs(
    title = "Top 30 Most Frequent Drugs",
    subtitle = paste0("Independent Top 30 per sample group; Ensemble_Score < ", sensitive_threshold),
    x = "Ratio (%)",
    y = NULL
  ) +

  ## Theme
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 5)),
    size = guide_legend(
      title.position = "left",
      direction = "horizontal",
      order = 2
    )
  ) +

  ## Group labels
  annotate(
    "text",
    x = -1.1,
    y = 1.5,
    label = "CRC cell lines",
    color = "#FF8C00",
    size = 5.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 1.5,
    label = "TCGA patients",
    color = "#4682B4",
    size = 5.5,
    fontface = "bold"
  )

############################################################
## 6. Export figure
############################################################

ggsave(
  filename = file.path(out_dir, "Fig5C_independent_top30_crc_vs_tcga.pdf"),
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Fig5C_independent_top30_crc_vs_tcga.png"),
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

############################################################
## Finished
############################################################

message("Fig. 5C mirror lollipop plot finished.")
message("Output directory: ", normalizePath(out_dir))
message("Output files:")
message("  - Fig5C_independent_top30_crc_vs_tcga.csv")
message("  - Fig5C_independent_top30_crc_vs_tcga.pdf")
message("  - Fig5C_independent_top30_crc_vs_tcga.png")
