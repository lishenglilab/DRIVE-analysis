library(tidyverse)

input_file <- "figureE_panel.csv"

type_colors <- c(
  "Chemotherapy" = "#F4A3A8",
  "Hormonal" = "#95D095",
  "Targeted" = "#9BB8E8"
)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
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
ggsave("figureE_panel.png", p, width = 11, height = 9, dpi = 600)
