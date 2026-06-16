############################################################
## supfig_s1_model_feature_summary.R
## Output: Sup1.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "supfig_s1_model_summary.csv"

## Required columns: Model, Omics_Feature, Drug_Representation, Algorithm_Type
df <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("Model", "Omics_Feature", "Drug_Representation", "Algorithm_Type")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

omics_df <- df %>%
  separate_rows(Omics_Feature, sep = "[;|,]") %>%
  mutate(Omics_Feature = trimws(Omics_Feature)) %>%
  filter(Omics_Feature != "") %>%
  distinct(Model, Omics_Feature) %>%
  count(Omics_Feature, name = "N") %>%
  arrange(desc(N))

drug_df <- df %>%
  separate_rows(Drug_Representation, sep = "[;|,]") %>%
  mutate(Drug_Representation = trimws(Drug_Representation)) %>%
  filter(Drug_Representation != "") %>%
  distinct(Model, Drug_Representation) %>%
  count(Drug_Representation, name = "N") %>%
  arrange(desc(N))

algo_df <- df %>%
  separate_rows(Algorithm_Type, sep = "[;|,]") %>%
  mutate(Algorithm_Type = trimws(Algorithm_Type)) %>%
  filter(Algorithm_Type != "") %>%
  distinct(Model, Algorithm_Type) %>%
  count(Algorithm_Type, name = "N") %>%
  arrange(desc(N))

p1 <- ggplot(omics_df, aes(x = reorder(Omics_Feature, N), y = N)) +
  geom_col(fill = "#4C78A8") +
  coord_flip() +
  labs(title = "A. Omics features", x = "", y = "Number of models") +
  theme_bw(base_size = 12)

p2 <- ggplot(drug_df, aes(x = reorder(Drug_Representation, N), y = N)) +
  geom_col(fill = "#F58518") +
  coord_flip() +
  labs(title = "B. Drug representations", x = "", y = "Number of models") +
  theme_bw(base_size = 12)

p3 <- ggplot(algo_df, aes(x = reorder(Algorithm_Type, N), y = N)) +
  geom_col(fill = "#54A24B") +
  coord_flip() +
  labs(title = "C. Algorithmic architectures", x = "", y = "Number of models") +
  theme_bw(base_size = 12)

pdf("Sup1.pdf", width = 10, height = 12)
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
dev.off()

message("Supplementary Fig. S1 summary finished.")
