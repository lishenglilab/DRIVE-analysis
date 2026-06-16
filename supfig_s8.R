############################################################
## supfig_s4_ensemble_strategy_comparison.R
## Output: Sup4.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "C3_loo_18_to_2_curves_all_strategies.csv"

## -------------------------------------------------------------------------
## 1. Read data
## Expected columns:
##   strategy, num_models, train_RMSE, validation_RMSE
## -------------------------------------------------------------------------
df <- read.csv(input_file, check.names = FALSE)

## -------------------------------------------------------------------------
## 2. Factor ordering for strategies
## -------------------------------------------------------------------------
df <- df %>%
  mutate(
    strategy = factor(
      strategy,
      levels = c("RandomForest", "RidgeCV", "LightGBM", "KNN", "XGBoost")
    )
  ) %>%
  arrange(strategy, num_models)

## -------------------------------------------------------------------------
## 3. Define color palette (Nature NPG style)
## -------------------------------------------------------------------------
npg_cols <- pal_npg("nrc")(5)
names(npg_cols) <- levels(df$strategy)

## -------------------------------------------------------------------------
## 4. Identify best points per strategy and overall
## -------------------------------------------------------------------------
train_best_each <- df %>%
  group_by(strategy) %>%
  slice_min(order_by = train_RMSE, n = 1, with_ties = FALSE) %>%
  ungroup()

valid_best_each <- df %>%
  group_by(strategy) %>%
  slice_min(order_by = validation_RMSE, n = 1, with_ties = FALSE) %>%
  ungroup()

## -------------------------------------------------------------------------
## 5. Define shared theme
## -------------------------------------------------------------------------
theme_nat <- theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

## -------------------------------------------------------------------------
## 6. Generate train and validation plots
## -------------------------------------------------------------------------
p_train <- ggplot(df, aes(x = num_models, y = train_RMSE, color = strategy)) +
  geom_line(linewidth = 1.05, lineend = "round") +
  geom_point(size = 2.4) +
  geom_point(
    data = train_best_each,
    aes(x = num_models, y = train_RMSE, color = strategy),
    shape = 21, fill = "white", stroke = 1.0, size = 3.4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = npg_cols) +
  scale_x_continuous(
    breaks = seq(min(df$num_models), max(df$num_models), by = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    title = "Train",
    x = "Number of models in ensemble",
    y = "Train RMSE",
    color = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_nat

p_valid <- ggplot(df, aes(x = num_models, y = validation_RMSE, color = strategy)) +
  geom_line(linewidth = 1.05, lineend = "round") +
  geom_point(size = 2.4) +
  geom_point(
    data = valid_best_each,
    aes(x = num_models, y = validation_RMSE, color = strategy),
    shape = 21, fill = "white", stroke = 1.0, size = 3.4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = npg_cols) +
  scale_x_continuous(
    breaks = seq(min(df$num_models), max(df$num_models), by = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    title = "Validation",
    x = "Number of models in ensemble",
    y = "Validation RMSE",
    color = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_nat

## -------------------------------------------------------------------------
## 7. Combine plots and save
## -------------------------------------------------------------------------
p_final <- p_train + p_valid +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

pdf("Sup4.pdf", width = 10, height = 5)
print(p_final)
dev.off()

message("Supplementary Fig. S4 generated: Sup4.pdf")