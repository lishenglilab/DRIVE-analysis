############################################################
## supfig_s3_violin_performance.R
## Output: Sup3.pdf
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

input_file <- "model_metric.csv"

## -------------------------------------------------------------------------
## 1. Read data
## Expected columns: Method, mode, RMSE, PCC, SCC, R2, NDCG, NWPC
## -------------------------------------------------------------------------
df <- read.csv(input_file, check.names = FALSE)

## -------------------------------------------------------------------------
## 2. Reshape to long format
## -------------------------------------------------------------------------
df_long <- df %>%
  pivot_longer(
    cols = c("RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    mode = dplyr::case_when(
      tolower(mode) == "mixed" ~ "Mixed",
      TRUE ~ mode
    )
  )

## -------------------------------------------------------------------------
## 3. Set factor levels for ordering
## -------------------------------------------------------------------------
df_long <- df_long %>%
  mutate(
    mode = factor(mode, levels = c("Mixed", "Cell-blind", "Drug-blind")),
    Metric = factor(Metric, levels = c("RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC"))
  )

## -------------------------------------------------------------------------
## 4. Generate violin plot
## -------------------------------------------------------------------------
p <- ggplot(df_long, aes(x = Metric, y = Value, fill = mode)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    position = position_dodge(0.8),
    color = "black",
    alpha = 0.8,
    linewidth = 0.3
  ) +
  facet_wrap(~ Method, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(
      "Mixed"       = "#E69F00",
      "Cell-blind"  = "#56B4E9",
      "Drug-blind"  = "#009E73"
    ),
    name = "Mode"
  ) +
  labs(
    y = "Performance Value",
    x = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, color = "black", size = 9, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11)
  )

## -------------------------------------------------------------------------
## 5. Save plot
## -------------------------------------------------------------------------
pdf("Sup3.pdf", width = 12, height = 12)
print(p)
dev.off()

message("Supplementary Fig. S3 generated: Sup3.pdf")
