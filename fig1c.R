library(tidyverse)

input_file <- "figureC_panel.csv"

method_order <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP", "GPDRP", "NERD", "PaccMann", "Precily",
  "GADRP", "DeepCCDS"
)

algo_order <- c("DNN", "GCN", "AE", "Transformer", "MCA", "KNN", "MLP", "GAT", "GIN", "GAT_GCN")

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    method = factor(method, levels = rev(method_order)),
    algorithm = factor(algorithm, levels = algo_order)
  )

grid_df <- expand_grid(
  method = factor(rev(method_order), levels = rev(method_order)),
  algorithm = factor(algo_order, levels = algo_order)
)

plot_df <- grid_df %>%
  left_join(df, by = c("method", "algorithm")) %>%
  mutate(
    present = replace_na(present, 0),
    fill_color = replace_na(fill_color, "#FFFFFF")
  )

p <- ggplot(plot_df, aes(x = algorithm, y = method)) +
  geom_tile(fill = "white", color = "#D9D9D9", linewidth = 1.0) +
  geom_point(
    data = plot_df %>% filter(present == 1),
    aes(fill = fill_color),
    shape = 21,
    size = 11,
    stroke = 1.0,
    color = "white"
  ) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "#D9D9D9", fill = NA, linewidth = 1.2)
  )

ggsave("figureC_panel.pdf", p, width = 8, height = 10)
ggsave("figureC_panel.svg", p, width = 8, height = 10)
ggsave("figureC_panel.png", p, width = 8, height = 10, dpi = 600)
