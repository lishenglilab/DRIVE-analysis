############################################################
## supfig_s7_internal_external_grouped_comparison.R
##
## Purpose:
##   Generate Supplementary Fig. S7 comparing internal and
##   external performance across random, tissue-holdout, and
##   chemical-cluster-holdout settings.
##
## Input:
##   ALL_split_level_internal_external_metrics.csv
##
## Required columns:
##   strategy, split, holdout_group, internal_RMSE, internal_PCC,
##   external_RMSE, external_PCC, n_test_cells, n_test_drugs,
##   external_n_cells, external_n_drugs
##
## Output:
##   1. DRIVE_internal_external_comparison_barplot.pdf
##   2. DRIVE_internal_external_comparison_barplot.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

work_dir <- get_script_dir()
setwd(work_dir)

metric_file <- "ALL_split_level_internal_external_metrics.csv"

if (!file.exists(metric_file)) {
  stop("Input file not found: ", metric_file)
}

df <- read_csv(metric_file, show_col_types = FALSE)

required_cols <- c(
  "strategy", "split", "holdout_group", "internal_RMSE", "internal_PCC",
  "external_RMSE", "external_PCC", "n_test_cells", "n_test_drugs",
  "external_n_cells", "external_n_drugs"
)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

tissue_groups_to_plot <- c("blood", "lung", "urogenital_system")
chemical_groups_to_plot <- c("Chem-A", "Chem-B", "Chem-C")

pretty_group <- function(x) {
  case_when(
    x == "blood" ~ "Blood",
    x == "lung" ~ "Lung",
    x == "urogenital_system" ~ "Urogenital",
    TRUE ~ str_replace_all(x, "_", " ") %>% str_to_title()
  )
}

save_pdf_png <- function(plot, file_prefix, width = 12, height = 9) {
  ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

random_internal <- df %>%
  filter(strategy == "random_5fold_cv") %>%
  arrange(split) %>%
  transmute(
    group = "Random",
    group_type = "Random",
    unit = split,
    RMSE = internal_RMSE,
    PCC = internal_PCC,
    n_label = "n=5"
  )

random_external <- df %>%
  filter(strategy == "random_5fold_cv") %>%
  arrange(split) %>%
  transmute(
    group = "Random",
    group_type = "Random",
    unit = split,
    RMSE = external_RMSE,
    PCC = external_PCC,
    n_label = "n=5"
  )

tissue_rows <- df %>%
  filter(strategy == "leave_tissue_organ_out")

if (!is.null(tissue_groups_to_plot)) {
  tissue_rows <- tissue_rows %>%
    filter(holdout_group %in% tissue_groups_to_plot)
}

chemical_rows <- df %>%
  filter(
    strategy == "leave_chemical_cluster_out",
    holdout_group %in% chemical_groups_to_plot
  )

grouped_internal <- bind_rows(
  tissue_rows %>%
    transmute(
      group = pretty_group(holdout_group),
      raw_group = holdout_group,
      group_type = "Tissue",
      RMSE = internal_RMSE,
      PCC = internal_PCC,
      n_label = paste0("n=", n_test_cells)
    ),
  chemical_rows %>%
    transmute(
      group = holdout_group,
      raw_group = holdout_group,
      group_type = "Chemical",
      RMSE = internal_RMSE,
      PCC = internal_PCC,
      n_label = paste0("n=", n_test_drugs)
    )
)

grouped_external <- bind_rows(
  tissue_rows %>%
    transmute(
      group = pretty_group(holdout_group),
      raw_group = holdout_group,
      group_type = "Tissue",
      RMSE = external_RMSE,
      PCC = external_PCC,
      n_label = paste0("n=", external_n_cells)
    ),
  chemical_rows %>%
    transmute(
      group = holdout_group,
      raw_group = holdout_group,
      group_type = "Chemical",
      RMSE = external_RMSE,
      PCC = external_PCC,
      n_label = paste0("n=", external_n_drugs)
    )
)

make_bar_df <- function(random_df, grouped_df, dataset_name) {
  random_sum <- random_df %>%
    summarise(
      group = "Random",
      group_type = "Random",
      RMSE = mean(RMSE, na.rm = TRUE),
      RMSE_sd = sd(RMSE, na.rm = TRUE),
      PCC = mean(PCC, na.rm = TRUE),
      PCC_sd = sd(PCC, na.rm = TRUE),
      n_label = "n=5"
    )

  grouped_sum <- grouped_df %>%
    mutate(
      RMSE_sd = NA_real_,
      PCC_sd = NA_real_
    ) %>%
    select(group, group_type, RMSE, RMSE_sd, PCC, PCC_sd, n_label)

  bind_rows(random_sum, grouped_sum) %>%
    mutate(dataset = dataset_name)
}

internal_df <- make_bar_df(random_internal, grouped_internal, "Internal validation")
external_df <- make_bar_df(random_external, grouped_external, "External validation")
plot_df <- bind_rows(internal_df, external_df)

group_order <- c(
  "Random",
  grouped_internal %>% filter(group_type == "Tissue") %>% pull(group),
  grouped_internal %>% filter(group_type == "Chemical") %>% pull(group)
) %>%
  unique()

plot_df <- plot_df %>%
  mutate(
    group = factor(group, levels = group_order),
    x_label = paste0(group, "\n", n_label)
  )

compare_colors <- c(
  "Internal validation" = "#4C97C2",
  "External validation" = "#D98C2B"
)

nature_theme <- theme_classic(base_size = 14, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.75),
    axis.ticks = element_line(color = "black", linewidth = 0.60),
    axis.ticks.length = unit(0.18, "cm"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      face = "bold",
      color = "black",
      size = 10.2,
      lineheight = 0.9
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 12
    ),
    axis.title.y = element_text(
      face = "bold",
      color = "black",
      size = 14
    ),
    axis.title.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16,
      color = "black"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11, color = "black"),
    legend.key.size = unit(0.45, "cm"),
    plot.margin = margin(20, 12, 10, 12)
  )

add_section_annotations <- function(plot_top, n_tissue, n_chemical) {
  x_random <- 1
  tissue_start <- 2
  tissue_end <- 1 + n_tissue
  chem_start <- tissue_end + 1
  chem_end <- tissue_end + n_chemical

  ann <- list()

  if (n_tissue > 0) {
    ann <- c(
      ann,
      list(
        geom_vline(
          xintercept = 1.5,
          linetype = "dashed",
          color = "grey50",
          linewidth = 0.45
        )
      )
    )
  }

  if (n_chemical > 0) {
    ann <- c(
      ann,
      list(
        geom_vline(
          xintercept = tissue_end + 0.5,
          linetype = "dashed",
          color = "grey50",
          linewidth = 0.45
        )
      )
    )
  }

  ann <- c(
    ann,
    list(
      annotate(
        "text",
        x = x_random,
        y = plot_top * 1.015,
        label = "Random 5-fold CV",
        fontface = "bold",
        size = 4.0
      ),
      annotate(
        "segment",
        x = 0.75,
        xend = 1.25,
        y = plot_top * 0.965,
        yend = plot_top * 0.965,
        linewidth = 0.55
      )
    )
  )

  if (n_tissue > 0) {
    ann <- c(
      ann,
      list(
        annotate(
          "text",
          x = mean(c(tissue_start, tissue_end)),
          y = plot_top * 1.015,
          label = "Leave-tissue/organ-out holdout",
          fontface = "bold",
          size = 4.0
        ),
        annotate(
          "segment",
          x = tissue_start - 0.35,
          xend = tissue_end + 0.35,
          y = plot_top * 0.965,
          yend = plot_top * 0.965,
          linewidth = 0.55
        )
      )
    )
  }

  if (n_chemical > 0) {
    ann <- c(
      ann,
      list(
        annotate(
          "text",
          x = mean(c(chem_start, chem_end)),
          y = plot_top * 1.015,
          label = "Leave-chemical-cluster-out holdout",
          fontface = "bold",
          size = 4.0
        ),
        annotate(
          "segment",
          x = chem_start - 0.35,
          xend = chem_end + 0.35,
          y = plot_top * 0.965,
          yend = plot_top * 0.965,
          linewidth = 0.55
        )
      )
    )
  }

  ann
}

n_tissue <- grouped_internal %>%
  filter(group_type == "Tissue") %>%
  distinct(group) %>%
  nrow()

n_chemical <- grouped_internal %>%
  filter(group_type == "Chemical") %>%
  distinct(group) %>%
  nrow()

make_comparison_panel <- function(metric = c("RMSE", "PCC")) {
  metric <- match.arg(metric)

  dat <- plot_df %>%
    mutate(
      group_chr = as.character(group),
      x_id = as.numeric(factor(group_chr, levels = group_order)),
      x_bar = case_when(
        dataset == "Internal validation" ~ x_id - 0.18,
        dataset == "External validation" ~ x_id + 0.18,
        TRUE ~ x_id
      )
    )

  random_points <- bind_rows(
    random_internal %>% mutate(dataset = "Internal validation"),
    random_external %>% mutate(dataset = "External validation")
  ) %>%
    mutate(
      group_chr = "Random",
      x_id = 1,
      x_bar = case_when(
        dataset == "Internal validation" ~ x_id - 0.18,
        dataset == "External validation" ~ x_id + 0.18,
        TRUE ~ x_id
      )
    ) %>%
    group_by(dataset) %>%
    mutate(
      x_point = x_bar + seq(-0.025, 0.025, length.out = n())
    ) %>%
    ungroup()

  value_fmt <- "%.3f"

  if (metric == "RMSE") {
    dat <- dat %>%
      mutate(
        value = RMSE,
        sd_value = RMSE_sd
      )

    random_points <- random_points %>%
      mutate(value = RMSE)

    y_lab <- "RMSE"
    plot_top <- max(dat$value + replace_na(dat$sd_value, 0), na.rm = TRUE) * 1.35
  } else {
    dat <- dat %>%
      mutate(
        value = PCC,
        sd_value = PCC_sd
      )

    random_points <- random_points %>%
      mutate(value = PCC)

    y_lab <- "PCC"
    plot_top <- 1.12
  }

  dat_random <- dat %>%
    filter(group_chr == "Random") %>%
    mutate(
      sd_value = replace_na(sd_value, 0),
      value_label = sprintf(value_fmt, value),
      label_y = value + sd_value + plot_top * 0.030
    )

  dat_other <- dat %>%
    filter(group_chr != "Random") %>%
    mutate(
      value_label = sprintf(value_fmt, value),
      label_y = value + plot_top * 0.025
    )

  ggplot() +
    geom_col(
      data = dat,
      aes(x = x_bar, y = value, fill = dataset),
      width = 0.32,
      color = "black",
      linewidth = 0.42
    ) +
    geom_errorbar(
      data = dat_random,
      aes(
        x = x_bar,
        ymin = value - sd_value,
        ymax = value + sd_value
      ),
      width = 0.08,
      linewidth = 0.55,
      color = "black"
    ) +
    geom_point(
      data = random_points,
      aes(x = x_point, y = value, fill = dataset),
      shape = 21,
      size = 0.75,
      stroke = 0.22,
      color = "black",
      inherit.aes = FALSE
    ) +
    add_section_annotations(plot_top, n_tissue, n_chemical) +
    scale_fill_manual(values = compare_colors) +
    scale_x_continuous(
      breaks = seq_along(group_order),
      labels = group_order,
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    coord_cartesian(ylim = c(0, plot_top), clip = "off") +
    labs(
      y = y_lab
    ) +
    nature_theme +
    geom_text(
      data = dat_random,
      aes(x = x_bar, y = label_y, label = value_label),
      size = 3.10,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = dat_other,
      aes(x = x_bar, y = label_y, label = value_label),
      size = 3.10,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    )
}

p_comp_rmse <- make_comparison_panel("RMSE")
p_comp_pcc <- make_comparison_panel("PCC")

fig_comparison <- p_comp_rmse / p_comp_pcc +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18,
        color = "black"
      ),
      plot.tag = element_text(
        face = "bold",
        size = 18,
        color = "black"
      )
    )
  )

print(fig_comparison)

save_pdf_png(
  plot = fig_comparison,
  file_prefix = "DRIVE_internal_external_comparison_barplot",
  width = 12,
  height = 9
)
