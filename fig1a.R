############################################################
## fig1a_model_omics_polar_ring_formal.R
##
## Purpose:
##   Formal cleaned script for Fig. 1A-style model-by-omics
##   polar ring plot.
##
## Key points:
##   1. Built-in test data are included for direct testing.
##   2. The omics legend is placed at the side.
##   3. The center hole is preserved.
##   4. The middle can be adjusted later in Illustrator/AI.
##
## Output:
##   figure1_outputs/Fig1A_model_omics_polar_ring.pdf
##   figure1_outputs/Fig1A_model_omics_polar_ring.png
##   figure1_outputs/fig1_model_omics_matrix_test.csv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

############################################################
## 1. User settings
############################################################

## TRUE  = use built-in test data
## FALSE = read external CSV file
use_test_data <- TRUE

## External file used only when use_test_data = FALSE
input_file <- "fig1_model_omics_matrix.csv"

out_dir <- "figure1_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Built-in test dataset
############################################################

method_order <- c(
  "DeepCCDS", "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK",
  "GraphDRP", "GPDRP", "NERD", "PaccMann", "Precily", "GADRP"
)

omics_order <- c(
  "Copy numbe variation",
  "Gene mutation",
  "mRNA expression",
  "DNA methylation",
  "miRNA expression",
  "Protein-protien interaction"
)

presence_list <- list(
  DeepCCDS = c("Gene mutation", "mRNA expression"),
  BANDRP   = c("Copy numbe variation", "Gene mutation", "mRNA expression", "DNA methylation"),
  DeepAEG  = c("Copy numbe variation", "Gene mutation", "mRNA expression"),
  DeepCDR  = c("Copy numbe variation", "Gene mutation", "mRNA expression"),
  DeepTTA  = c("mRNA expression"),
  DIPK     = c("mRNA expression", "Protein-protien interaction"),
  GraphDRP = c("Gene mutation", "mRNA expression"),
  GPDRP    = c("mRNA expression"),
  NERD     = c("Copy numbe variation", "mRNA expression", "miRNA expression"),
  PaccMann = c("mRNA expression"),
  Precily  = c("mRNA expression"),
  GADRP    = c("Copy numbe variation", "mRNA expression", "DNA methylation", "miRNA expression")
)

build_test_data <- function(methods, omics, presence_list) {
  expand_grid(Method = methods, Omics = omics) %>%
    rowwise() %>%
    mutate(Present = ifelse(Omics %in% presence_list[[Method]], 1, 0)) %>%
    ungroup()
}

############################################################
## 3. Load data
############################################################

if (use_test_data) {
  message("Using built-in test dataset.")
  df <- build_test_data(method_order, omics_order, presence_list)
  write.csv(df, file.path(out_dir, "fig1_model_omics_matrix_test.csv"), row.names = FALSE)
} else {
  message("Reading external input file: ", input_file)
  df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
}

required_cols <- c("Method", "Omics", "Present")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

############################################################
## 4. Format data
############################################################

## Main geometry controls
inner_radius <- 2.4
ring_width   <- 0.78
outer_margin <- 0.75

df <- df %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Omics = factor(Omics, levels = omics_order),
    Present = as.numeric(Present),
    FillGroup = ifelse(Present > 0, as.character(Omics), NA),
    method_id = as.numeric(Method),
    omics_id = as.numeric(Omics),
    ring_y = inner_radius + omics_id
  ) %>%
  filter(!is.na(Method), !is.na(Omics))

omics_colors <- c(
  "Copy numbe variation"        = "#FF5A5F",
  "Gene mutation"               = "#5B5CF6",
  "mRNA expression"             = "#5B96C8",
  "DNA methylation"             = "#D58AE6",
  "miRNA expression"            = "#FF929B",
  "Protein-protien interaction" = "#9A9ACB"
)

n_methods <- length(method_order)
n_omics   <- length(omics_order)

sector_lines <- tibble(
  x = seq(0.5, n_methods + 0.5, by = 1),
  y_start = inner_radius + 0.5,
  y_end   = inner_radius + n_omics + 0.5
)

circle_lines <- tibble(
  y = c(inner_radius + 0.5, inner_radius + n_omics + 0.5)
)

############################################################
## 5. Plot
############################################################

p <- ggplot() +
  ## Background tiles
  geom_tile(
    data = df,
    aes(x = method_id, y = ring_y),
    fill = "#EDEDED",
    color = "white",
    linewidth = 0.8,
    width = 0.96,
    height = ring_width
  ) +
  ## Colored tiles
  geom_tile(
    data = df %>% filter(Present > 0),
    aes(x = method_id, y = ring_y, fill = FillGroup),
    color = "white",
    linewidth = 0.8,
    width = 0.96,
    height = ring_width
  ) +
  ## Sector lines
  geom_segment(
    data = sector_lines,
    aes(x = x, xend = x, y = y_start, yend = y_end),
    color = "black",
    linewidth = 0.55
  ) +
  ## Inner/outer circular borders
  geom_hline(
    data = circle_lines,
    aes(yintercept = y),
    color = "black",
    linewidth = 0.65
  ) +
  coord_polar(start = 0, clip = "off") +
  scale_x_continuous(
    breaks = seq_along(method_order),
    labels = method_order,
    limits = c(0.5, n_methods + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, inner_radius + n_omics + outer_margin),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = omics_colors,
    na.value = "#EDEDED",
    name = "Omics",
    drop = FALSE
  ) +
  labs(x = "", y = "") +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold.italic", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(12, 24, 12, 12)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(color = NA),
      ncol = 1
    )
  )

############################################################
## 6. Export
############################################################

ggsave(
  filename = file.path(out_dir, "Fig1A_model_omics_polar_ring.pdf"),
  plot = p,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Fig1A_model_omics_polar_ring.png"),
  plot = p,
  width = 9,
  height = 8,
  dpi = 300
)

message("Fig. 1A formal plot finished.")
message("Output directory: ", normalizePath(out_dir))
message("Output files:")
message("  - Fig1A_model_omics_polar_ring.pdf")
message("  - Fig1A_model_omics_polar_ring.png")
message("  - fig1_model_omics_matrix_test.csv")
