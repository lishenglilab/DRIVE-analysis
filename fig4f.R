library(tidyverse)
library(scales)

input_file <- "plot_data_for_R.csv"
output_pdf <- "Fig4F_stacked_bar.pdf"
output_png <- "Fig4F_stacked_bar.png"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

final_plot_data <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("tissue", "pathway", "count")
missing_cols <- setdiff(required_cols, colnames(final_plot_data))

if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

plot_data <- final_plot_data %>%
  mutate(
    tissue = str_trim(as.character(tissue)),
    pathway = str_trim(as.character(pathway)),
    count = as.numeric(count)
  ) %>%
  filter(
    !is.na(tissue), tissue != "",
    !is.na(pathway), pathway != "",
    !is.na(count), count > 0
  ) %>%
  group_by(tissue, pathway) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(
    total_count = sum(count),
    proportion = count / total_count
  ) %>%
  ungroup()

pathway_levels <- plot_data %>%
  group_by(pathway) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  pull(pathway)

tissue_levels <- plot_data %>%
  group_by(tissue) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  pull(tissue)

plot_data <- plot_data %>%
  mutate(
    pathway = factor(pathway, levels = pathway_levels),
    tissue = factor(tissue, levels = tissue_levels)
  )

palette_values <- c(
  "Amino acids and Peptides" = "#A65628",
  "Carbohydrates" = "#636363",
  "Unclassified" = "#E7298A",
  "Polyketides" = "#377EB8",
  "Fatty acids" = "#FFFF99",
  "Shikimates and Phenylpropanoids" = "#FDB462",
  "Alkaloids" = "#BEBADA",
  "Terpenoids" = "#A1D99B"
)

missing_pathways <- setdiff(levels(plot_data$pathway), names(palette_values))
if (length(missing_pathways) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_pathways))
  names(extra_colors) <- missing_pathways
  palette_values <- c(palette_values, extra_colors)
}

p <- ggplot(plot_data, aes(x = tissue, y = proportion, fill = pathway)) +
  geom_col(width = 0.72, color = NA) +
  scale_fill_manual(values = palette_values, name = "Compound Class") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Tissue",
    y = "Proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.18, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )

print(p)

ggsave(output_pdf, p, width = 10, height = 8)
ggsave(output_png, p, width = 10, height = 8, dpi = 600)
