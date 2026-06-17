library(tidyverse)

input_file <- "figureB_panel.csv"

method_order <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP", "GPDRP", "NERD", "PaccMann", "Precily",
  "GADRP", "DeepCCDS"
)

feature_order <- c("Molecular graph", "Fingerprint", "ESPF", "ECFP", "Morgan", "SMILESVec")

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    method = factor(method, levels = rev(method_order)),
    feature = factor(feature, levels = feature_order)
  )

grid_df <- expand_grid(
  method = factor(rev(method_order), levels = rev(method_order)),
  feature = factor(feature_order, levels = feature_order)
)

plot_df <- grid_df %>%
  left_join(df, by = c("method", "feature")) %>%
  mutate(
    present = replace_na(present, 0),
    fill_color = replace_na(fill_color, "#D9D9D9")
  )

p <- ggplot(plot_df, aes(x = feature, y = method)) +
  geom_tile(aes(fill = fill_color), color = "white", linewidth = 1.0) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

ggsave("figureB_panel.pdf", p, width = 6, height = 10)
ggsave("figureB_panel.svg", p, width = 6, height = 10)
ggsave("figureB_panel.png", p, width = 6, height = 10, dpi = 600)
