############################################################
## supfig_s14_prediction_validation_scatter.R
## Output: Sup14.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

## -------------------------------------------------------------------------
## 1. Read input files
## -------------------------------------------------------------------------
cell_info <- read.csv("cellline.csv", stringsAsFactors = FALSE)
pred_data <- read.csv("ml_ensemble_predictions.csv", stringsAsFactors = FALSE)
colnames(pred_data) <- c("cell_line", "drug", "ic50")

## GDSC experimental data (wide format, rownames = cell lines)
gdsc_raw <- read.csv("gdsc_ic50.csv", row.names = 1, check.names = FALSE)

## DepMap experimental data (wide format, rownames = ModelID)
depmap_raw <- read.csv("depmap.csv", row.names = 1, check.names = FALSE)

## -------------------------------------------------------------------------
## 2. Filter COAD cell lines
## -------------------------------------------------------------------------
coad_meta <- cell_info %>%
  filter(DepmapModelType == "COAD") %>%
  select(StrippedCellLineName, ModelID)

target_stripped <- coad_meta$StrippedCellLineName
target_modelid <- coad_meta$ModelID

## -------------------------------------------------------------------------
## 3. Drug name mapping (manual + cleaning)
## -------------------------------------------------------------------------
drug_name_mapping <- c(
  "CPP"          = "MCPP",
  "BF2.649"      = "PITOLISANT",
  "L-798,106"    = "L-798106",
  "MLN-8054"     = "MLN8054",
  "ONC201"       = "TIC10",
  "PLX-4720"     = "PLX4720",
  "XL-647"       = "XL647",
  "blebbistatin-(-)" = "BLEBBISTATIN-(+/-)",
  "cyclosporine" = "CYCLOSPORIN-A",
  "lorlatinib"   = "PF-06463922"
)

clean_name <- function(name) {
  name <- tolower(name)
  gsub("[^a-z0-9]", "", name)
}

get_mapping <- function(pred_drugs, true_drugs) {
  pred_clean <- sapply(pred_drugs, function(x) {
    if (x %in% names(drug_name_mapping)) x <- drug_name_mapping[x]
    clean_name(x)
  })
  true_clean <- sapply(true_drugs, clean_name)
  
  mapping <- data.frame(
    Predicted_Drug = character(),
    True_Drug = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(pred_clean)) {
    t_idx <- which(true_clean == pred_clean[i])
    if (length(t_idx) > 0) {
      for (j in t_idx) {
        mapping <- rbind(mapping, data.frame(
          Predicted_Drug = pred_drugs[i],
          True_Drug = true_drugs[j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  return(mapping)
}

## -------------------------------------------------------------------------
## 4. Prepare GDSC plot data
## -------------------------------------------------------------------------
## GDSC long format
gdsc_long <- gdsc_raw %>%
  as.data.frame() %>%
  mutate(cell_line = rownames(.)) %>%
  pivot_longer(
    cols = -cell_line,
    names_to = "drug_name",
    values_to = "gdsc_ic50"
  ) %>%
  mutate(gdsc_ic50 = as.numeric(gdsc_ic50)) %>%
  drop_na(gdsc_ic50)

## Filter predictions to COAD
pred_coad <- pred_data %>% filter(cell_line %in% target_stripped)

## Build GDSC mapping
gdsc_drugs <- unique(gdsc_long$drug_name)
pred_drugs <- unique(pred_coad$drug)
gdsc_mapping <- get_mapping(pred_drugs, gdsc_drugs)

## Merge GDSC data
plot_data_gdsc <- pred_coad %>%
  inner_join(gdsc_mapping, by = c("drug" = "Predicted_Drug")) %>%
  inner_join(gdsc_long, by = c("cell_line", "True_Drug" = "drug_name")) %>%
  rename(Experimental = gdsc_ic50, Predicted = ic50)

message("GDSC plot data points: ", nrow(plot_data_gdsc))

## -------------------------------------------------------------------------
## 5. Prepare DepMap plot data
## -------------------------------------------------------------------------
## Filter DepMap to COAD
depmap_coad <- depmap_raw[rownames(depmap_raw) %in% target_modelid, ]

## DepMap long format
depmap_long <- depmap_coad %>%
  rownames_to_column("ModelID") %>%
  pivot_longer(
    cols = -ModelID,
    names_to = "True_Drug",
    values_to = "depmap_ic50"
  ) %>%
  filter(!is.na(depmap_ic50))

## Build DepMap mapping
depmap_drugs <- colnames(depmap_coad)
depmap_mapping <- get_mapping(pred_drugs, depmap_drugs)

## Merge DepMap data
plot_data_depmap <- pred_coad %>%
  left_join(coad_meta, by = c("cell_line" = "StrippedCellLineName")) %>%
  inner_join(depmap_mapping, by = c("drug" = "Predicted_Drug")) %>%
  inner_join(depmap_long, by = c("ModelID", "True_Drug")) %>%
  rename(Experimental = depmap_ic50, Predicted = ic50)

## Filter outliers (-15 to 15)
plot_data_depmap <- plot_data_depmap %>%
  filter(Predicted >= -15 & Predicted <= 15,
         Experimental >= -15 & Experimental <= 15)

message("DepMap plot data points (filtered): ", nrow(plot_data_depmap))

## -------------------------------------------------------------------------
## 6. Define plotting function
## -------------------------------------------------------------------------
create_scatter <- function(data, x, y, xlab, ylab, title, subtitle = NULL,
                           xlim = NULL, ylim = NULL, point_color = "#005BAC",
                           cor_method = "spearman") {
  
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = point_color, alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", color = "#E74C3C", fill = "grey80", alpha = 0.3) +
    stat_cor(
      method = cor_method,
      aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
      label.x.npc = 0.05,
      label.y.npc = 0.95,
      size = 5
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    theme_classic(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30")
    )
  
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}

## -------------------------------------------------------------------------
## 7. Generate plots
## -------------------------------------------------------------------------
p_gdsc <- create_scatter(
  plot_data_gdsc,
  x = "Experimental",
  y = "Predicted",
  xlab = "LNIC50",
  ylab = "DRIVE Predicted IC50",
  title = "CRC cell lines in GDSC",
  subtitle = "R = 0.71; P < 2.2E-16",
  point_color = "#2E86C1"
)

p_depmap <- create_scatter(
  plot_data_depmap,
  x = "Experimental",
  y = "Predicted",
  xlab = "LNIC50",
  ylab = "DRIVE Predicted IC50",
  title = "CRC cell lines in DepMap",
  subtitle = "R = 0.24; P < 2.2E-16",
  point_color = "#005BAC",
  xlim = c(-15, 15),
  ylim = c(-15, 15)
)

## -------------------------------------------------------------------------
## 8. Combine and save
## -------------------------------------------------------------------------
p_combined <- p_gdsc + p_depmap +
  plot_annotation(
    title = "Supplementary Fig. S14",
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(face = "bold", size = 16))

pdf("Sup14.pdf", width = 14, height = 6)
print(p_combined)
dev.off()

message("Supplementary Fig. S14 generated: Sup14.pdf")