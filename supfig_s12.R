library(ggplot2)
library(dplyr)

# ============================================================
# Load data
# ============================================================

ic50_df <- read.csv(
  "path/to/predicted_ic50_by_tissue.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# Data preparation
# ============================================================

ic50_df <- ic50_df %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = unique(Tissue)
    )
  )

# ============================================================
# Visualization
# ============================================================

p <- ggplot(
  ic50_df,
  aes(
    x = Tissue,
    y = IC50
  )
) +
  geom_violin(
    fill = "#73A6D1",
    color = "black",
    linewidth = 0.5,
    trim = FALSE,
    scale = "width"
  ) +
  labs(
    x = NULL,
    y = expression(
      "Predicted IC"[50] ~ "values"
    )
  ) +
  scale_y_continuous(
    limits = c(-7, 9),
    breaks = c(-5, 0, 5)
  ) +
  theme(
    panel.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 1
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 12
    ),
    
    axis.title.y = element_text(
      size = 16
    ),
    
    axis.title.x = element_blank(),
    
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    )
  )

p

ggsave(
  filename = "Predicted_IC50_by_Tissue.pdf",
  plot = p,
  width = 8,
  height = 6
)