############################################################
## fig4d_count_size_dotplot.R
## Output: count_size_dotplot.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "fig4d_coverage_ratio.csv"
out_dir <- "fig4d_count_size_dotplot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Required columns: Tissue, Data_Type, Count_ratio
## Optional column: Count
df <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("Tissue", "Data_Type", "Count_ratio")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!("Count" %in% colnames(df))) {
  df$Count <- df$Count_ratio
}

df <- df %>%
  mutate(
    Count_ratio = as.numeric(Count_ratio),
    Count = as.numeric(Count),
    Tissue = fct_reorder(Tissue, Count_ratio, .fun = median, .desc = FALSE)
  )

p <- ggplot(df, aes(x = Data_Type, y = Tissue)) +
  geom_point(aes(size = Count_ratio, fill = Count),
             shape = 21, color = "grey35", alpha = 0.85) +
  scale_size_continuous(name = "Count ratio", range = c(2, 10)) +
  scale_fill_gradient(name = "Count", low = "#deebf7", high = "#08519c") +
  labs(x = "", y = "") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "count_size_dotplot.pdf"),
       p, width = 8, height = max(6, 0.25 * length(unique(df$Tissue))), dpi = 300)

message("Fig. 4D count-size dotplot finished.")
