library(tidyverse)

input_file <- "performance_data.csv"

data <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("Cell_N", "Time", "Peak_RAM")
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

data_long <- data %>%
  select(Cell_N, Time, Peak_RAM) %>%
  pivot_longer(
    cols = c(Time, Peak_RAM),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Cell_N = as.numeric(Cell_N),
    value = as.numeric(value)
  )

p <- ggplot(data_long, aes(x = Cell_N, y = value, color = variable)) +
  geom_line(linewidth = 2) +
  geom_point(
    aes(fill = variable),
    size = 5,
    shape = 21,
    color = "black",
    stroke = 1.2
  ) +
  facet_wrap(~ variable, scales = "free", ncol = 1) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 50, 100, 500, 1000),
    labels = c(1, 2, 5, 10, 50, 100, 500, 1000)
  ) +
  scale_color_manual(values = c("Time" = "#5494cc", "Peak_RAM" = "#e18283")) +
  scale_fill_manual(values = c("Time" = "#5494cc", "Peak_RAM" = "#e18283")) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.2),
    axis.text = element_text(color = "black", face = "bold"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(
    x = "Number of Points",
    y = "Measured Value"
  )

print(p)