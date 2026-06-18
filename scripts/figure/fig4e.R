############################################################
## fig4e_drug_count_matrix.R
## Output: drug_count_matrix.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
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

input_file <- "fig4e_drug_count_data.csv"
out_dir <- "fig4e_drug_count_matrix"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Required columns: OncotreeLineage, Drug_Category, Count
## Optional column: Count_ratio
df <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("OncotreeLineage", "Drug_Category", "Count")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!("Count_ratio" %in% colnames(df))) {
  df <- df %>%
    group_by(OncotreeLineage) %>%
    mutate(Count_ratio = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup()
}

df <- df %>%
  mutate(Count = as.numeric(Count), Count_ratio = as.numeric(Count_ratio))

matrix_df <- df %>%
  select(OncotreeLineage, Drug_Category, Count_ratio) %>%
  pivot_wider(names_from = Drug_Category, values_from = Count_ratio, values_fill = 0)

write.csv(matrix_df, file.path(out_dir, "fig4e_drug_count_matrix.csv"), row.names = FALSE)

tissue_order <- df %>%
  group_by(OncotreeLineage) %>%
  summarise(total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  arrange(total) %>%
  pull(OncotreeLineage)

category_order <- df %>%
  group_by(Drug_Category) %>%
  summarise(total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Drug_Category)

df <- df %>%
  mutate(
    OncotreeLineage = factor(OncotreeLineage, levels = tissue_order),
    Drug_Category = factor(Drug_Category, levels = category_order)
  )

p <- ggplot(df, aes(x = Drug_Category, y = OncotreeLineage, fill = Count_ratio)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "Count ratio") +
  labs(x = "Drug category", y = "Tissue lineage") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "drug_count_matrix.pdf"),
       p,
       width = max(8, 0.45 * length(category_order)),
       height = max(6, 0.25 * length(tissue_order)),
       dpi = 300,
       limitsize = FALSE)

message("Fig. 4E drug count matrix finished.")
