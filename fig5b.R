library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(tibble)

# ============================================================
# Load data
# ============================================================

mat <- read.csv(
  "path/to/drug_response_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

# rows   : drugs
# columns: cell lines
# values : sensitivity score / predicted IC50

# ============================================================
# Convert to long format
# ============================================================

heatmap_df <- mat %>%
  rownames_to_column("drug_name") %>%
  pivot_longer(
    -drug_name,
    names_to = "cell_line",
    values_to = "Sensitivity_Score"
  )

# ============================================================
# Preserve original ordering
# ============================================================

heatmap_df$drug_name <- factor(
  heatmap_df$drug_name,
  levels = rev(rownames(mat))
)

heatmap_df$cell_line <- factor(
  heatmap_df$cell_line,
  levels = colnames(mat)
)

# ============================================================
# Heatmap
# ============================================================

p <- ggplot(
  heatmap_df,
  aes(
    x = cell_line,
    y = drug_name,
    fill = Sensitivity_Score
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  scale_fill_gradientn(
    colours = c(
      "#F0E8D0",  # light beige
      "#F4C0AC",  # soft pink
      "#F58A73"   # coral
    ),
    values = rescale(
      c(-5, -3.5, -2)
    ),
    limits = c(-5, -2),
    breaks = c(-5, -4, -3, -2),
    name = "Sensitivity\nScore"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.ticks = element_blank(),
    
    axis.title = element_blank(),
    
    panel.border = element_blank(),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )

p

# ============================================================
# Save figure
# ============================================================

ggsave(
  "Heatmap_minus1.pdf",
  p,
  width = 12,
  height = 8,
  units = "in"
)