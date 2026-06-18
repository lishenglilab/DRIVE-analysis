options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
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

set.seed(42)

overview_entries <- list()

record_overview <- function(file, obj) {
  rows <- if (is.data.frame(obj) || is.matrix(obj)) nrow(obj) else NA_integer_
  cols <- if (is.data.frame(obj) || is.matrix(obj)) ncol(obj) else NA_integer_
  overview_entries[[length(overview_entries) + 1]] <<- data.frame(
    file = file,
    rows = rows,
    columns = cols
  )
}

write_csv_record <- function(obj, file, row.names = FALSE) {
  write.csv(obj, file, row.names = row.names)
  record_overview(file, obj)
}

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

methods <- c(
  "BANDRP", "DeepAEG", "DeepCDR", "DeepTTA", "DIPK", "GraphDRP",
  "GPDRP", "NERD", "PaccMann", "Precily", "GADRP", "DeepCCDS"
)

features <- c("Molecular graph", "Fingerprint", "ESPF", "ECFP", "Morgan", "SMILESVec")
algorithms <- c("DNN", "GCN", "AE", "Transformer", "MCA", "KNN", "MLP", "GAT", "GIN", "GAT_GCN")
omics <- c(
  "Copy number variation",
  "Gene mutation",
  "mRNA expression",
  "DNA methylation",
  "miRNA expression",
  "Protein-protein interaction"
)

feature_colors <- c(
  "Molecular graph" = "#5B8FF9",
  "Fingerprint" = "#61DDAA",
  "ESPF" = "#65789B",
  "ECFP" = "#F6BD16",
  "Morgan" = "#7262FD",
  "SMILESVec" = "#78D3F8"
)

algorithm_colors <- c(
  "DNN" = "#5470C6",
  "GCN" = "#91CC75",
  "AE" = "#FAC858",
  "Transformer" = "#EE6666",
  "MCA" = "#73C0DE",
  "KNN" = "#3BA272",
  "MLP" = "#FC8452",
  "GAT" = "#9A60B4",
  "GIN" = "#EA7CCC",
  "GAT_GCN" = "#6E7074"
)

presence_map <- list(
  BANDRP = c("Copy number variation", "Gene mutation", "mRNA expression", "DNA methylation"),
  DeepAEG = c("Copy number variation", "Gene mutation", "mRNA expression"),
  DeepCDR = c("Copy number variation", "Gene mutation", "mRNA expression"),
  DeepTTA = c("mRNA expression"),
  DIPK = c("mRNA expression", "Protein-protein interaction"),
  GraphDRP = c("Gene mutation", "mRNA expression"),
  GPDRP = c("mRNA expression"),
  NERD = c("Copy number variation", "mRNA expression", "miRNA expression"),
  PaccMann = c("mRNA expression"),
  Precily = c("mRNA expression"),
  GADRP = c("Copy number variation", "mRNA expression", "DNA methylation", "miRNA expression"),
  DeepCCDS = c("Gene mutation", "mRNA expression")
)

fig1a_df <- do.call(
  rbind,
  lapply(methods, function(method) {
    data.frame(
      Method = method,
      Omics = omics,
      Present = as.integer(omics %in% presence_map[[method]])
    )
  })
)
write_csv_record(fig1a_df, "fig1_model_omics_matrix.csv")

fig1b_df <- data.frame(
  method = rep(methods, each = length(features)),
  feature = rep(features, times = length(methods)),
  present = c(
    0, 1, 1, 1, 0, 0,
    1, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0
  )
)
write_csv_record(fig1b_df, "figureB_panel.csv")

fig1c_df <- do.call(
  rbind,
  lapply(seq_along(methods), function(i) {
    active_algorithms <- algorithms[c(1, 2, 4)]
    if (i %% 2 == 0) {
      active_algorithms <- c(active_algorithms, "GAT")
    }
    if (i %% 4 == 0) {
      active_algorithms <- c(active_algorithms, "GIN")
    }
    data.frame(
      method = methods[i],
      algorithm = active_algorithms,
      present = 1,
      fill_color = unname(algorithm_colors[active_algorithms])
    )
  })
)
write_csv_record(fig1c_df, "figureC_panel.csv")

system_palette <- c(
  "Lung" = "#A7CEE5",
  "Blood" = "#3D8CC6",
  "Urogenital system" = "#7CC495",
  "Digestive system" = "#64B84A",
  "Aero digestive tract" = "#86A850",
  "Nervous system" = "#F28A8A",
  "Breast" = "#E61E1E",
  "Skin" = "#F2B26E",
  "Pancreas" = "#FF8418",
  "Kidney" = "#D99979",
  "Bone" = "#9A82D0",
  "Soft tissue" = "#9A7AA7",
  "Thyroid" = "#EED96B",
  "Large intestine" = "#B85E28"
)

subtypes <- data.frame(
  system_group = c(
    "Lung", "Lung", "Blood", "Urogenital system", "Digestive system", "Aero digestive tract",
    "Nervous system", "Breast", "Skin", "Pancreas", "Kidney", "Bone",
    "Soft tissue", "Thyroid", "Large intestine", "Large intestine"
  ),
  subtype = c(
    "NSCLC", "SCLC", "AML", "Ovary", "Stomach", "Esophagus",
    "Glioma", "TNBC", "Melanoma", "PDAC", "RCC", "Osteosarcoma",
    "Sarcoma", "ATC", "CRC", "Rectum"
  ),
  count = c(18, 9, 7, 8, 11, 6, 5, 10, 7, 4, 5, 3, 4, 2, 12, 8)
)
subtypes$fill_color <- unname(system_palette[subtypes$system_group])
write_csv_record(subtypes, "figureD_panel.csv")

fig1e_df <- data.frame(
  pathway = c("EGFR", "MAPK", "PI3K/AKT", "DNA damage", "Apoptosis", "Cell cycle", "HDAC"),
  drug_count = c(8, 6, 7, 5, 4, 3, 2),
  drug_type = c("Targeted", "Targeted", "Targeted", "Chemotherapy", "Hormonal", "Chemotherapy", "Hormonal")
)
write_csv_record(fig1e_df, "figureE_panel.csv")

settings <- c("Mixed", "Drug-blind", "Cell-blind")
setting_offsets <- c("Mixed" = 0.00, "Drug-blind" = 0.06, "Cell-blind" = 0.09)
method_index <- seq_along(methods)

benchmark_df <- do.call(
  rbind,
  lapply(settings, function(setting) {
    do.call(
      rbind,
      lapply(1:3, function(rep_id) {
        data.frame(
          Method = methods,
          Setting = setting,
          RMSE = round(0.62 + setting_offsets[[setting]] + method_index * 0.015 + rnorm(length(methods), 0, 0.01), 3),
          PCC = round(0.82 - setting_offsets[[setting]] - method_index * 0.01 + rnorm(length(methods), 0, 0.008), 3),
          SCC = round(0.79 - setting_offsets[[setting]] - method_index * 0.008 + rnorm(length(methods), 0, 0.008), 3),
          R2 = round(0.68 - setting_offsets[[setting]] - method_index * 0.012 + rnorm(length(methods), 0, 0.008), 3),
          NDCG = round(0.86 - setting_offsets[[setting]] - method_index * 0.009 + rnorm(length(methods), 0, 0.008), 3),
          NWPC = round(0.84 - setting_offsets[[setting]] - method_index * 0.009 + rnorm(length(methods), 0, 0.008), 3)
        )
      })
    )
  })
)

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", benchmark_df)
saveWorkbook(wb, "synthetic_benchmark_settings.xlsx", overwrite = TRUE)
record_overview("synthetic_benchmark_settings.xlsx", benchmark_df)

rank_df <- benchmark_df %>%
  group_by(Method) %>%
  summarise(
    RMSE = round(mean(RMSE), 3),
    PCC = round(mean(PCC), 3),
    SCC = round(mean(SCC), 3),
    R2 = round(mean(R2), 3),
    NDCG = round(mean(NDCG), 3),
    NWPC = round(mean(NWPC), 3),
    .groups = "drop"
  )
write_csv_record(rank_df, "final_metrics_evaluation_ELITE_INTERSECTION.csv")

canonical_drugs <- c(
  "MCPP", "PITOLISANT", "L-798106", "MLN8054", "TIC10",
  "PLX4720", "XL647", "BLEBBISTATIN-(+/-)", "CYCLOSPORIN-A", "PF-06463922"
)

predicted_drug_alias <- c(
  "CPP", "BF2.649", "L-798,106", "MLN-8054", "ONC201",
  "PLX-4720", "XL-647", "blebbistatin-(-)", "cyclosporine", "lorlatinib"
)

stripped_cells <- c("HCT116", "RKO", "Lovo", "SW480", "DLD1", "HT29")
model_ids <- paste0("ACH-", sprintf("%06d", seq_along(stripped_cells)))
canonical_ids <- paste0("CL", sprintf("%03d", seq_along(stripped_cells)))
cosmic_ids <- as.character(9001:(9000 + length(stripped_cells)))

lineage_map <- data.frame(
  cell_line = stripped_cells,
  OncotreeLineage = c("Large Intestine", "Large Intestine", "Large Intestine", "Large Intestine", "Pancreas", "Stomach"),
  ModelID = model_ids,
  CanonicalID = canonical_ids,
  CosmicID = cosmic_ids
)

drug_base <- seq(-2.4, 1.2, length.out = length(canonical_drugs))
cell_effect <- c(-0.6, -0.2, 0.0, 0.2, 0.45, 0.7)

truth_grid <- expand.grid(
  cell_idx = seq_along(stripped_cells),
  drug_idx = seq_along(canonical_drugs)
)

truth_grid$gdsc_value <- round(
  drug_base[truth_grid$drug_idx] + cell_effect[truth_grid$cell_idx] + rnorm(nrow(truth_grid), 0, 0.15),
  3
)
truth_grid$depmap_value <- round(truth_grid$gdsc_value + rnorm(nrow(truth_grid), 0, 0.22), 3)
truth_grid$canonical_drug <- canonical_drugs[truth_grid$drug_idx]
truth_grid$pred_drug <- predicted_drug_alias[truth_grid$drug_idx]
truth_grid$cell_line <- stripped_cells[truth_grid$cell_idx]
truth_grid$model_id <- model_ids[truth_grid$cell_idx]
truth_grid$canonical_id <- canonical_ids[truth_grid$cell_idx]
truth_grid$cosmic_id <- cosmic_ids[truth_grid$cell_idx]
truth_grid$OncotreeLineage <- lineage_map$OncotreeLineage[truth_grid$cell_idx]

gdsc_matrix <- xtabs(gdsc_value ~ cell_line + canonical_drug, data = truth_grid)
depmap_matrix <- xtabs(depmap_value ~ model_id + canonical_drug, data = truth_grid)

write_csv_record(as.data.frame.matrix(gdsc_matrix), "gdsc_ic50.csv", row.names = TRUE)
write_csv_record(as.data.frame.matrix(depmap_matrix), "depmap.csv", row.names = TRUE)

cellline_df <- lineage_map %>%
  transmute(
    DepmapModelType = "COAD",
    StrippedCellLineName = cell_line,
    ModelID = ModelID
  )
write_csv_record(cellline_df, "cellline.csv")

cellline2_df <- data.frame(
  ID = canonical_ids,
  `cell.names` = stripped_cells,
  `cosmic.id` = cosmic_ids,
  `cell.names_nerd` = tolower(stripped_cells),
  `cell.names.1` = model_ids
)
write_csv_record(cellline2_df, "cellline2.csv", row.names = TRUE)

ic50_training <- data.frame(
  CellLineID = canonical_ids,
  as.data.frame.matrix(xtabs(gdsc_value ~ canonical_id + canonical_drug, data = truth_grid)),
  check.names = FALSE
)
write_csv_record(ic50_training, "ic50.csv")

ml_ensemble_predictions <- truth_grid %>%
  transmute(
    cell_line = cell_line,
    drug_name = pred_drug,
    Ensemble_Score = round(gdsc_value + rnorm(n(), 0, 0.18), 3)
  )
write_csv_record(ml_ensemble_predictions, "ml_ensemble_predictions.csv")

high_quality_data <- truth_grid %>%
  transmute(
    OncotreeLineage = OncotreeLineage,
    cell_line = cell_line,
    drug = canonical_drug,
    ic50 = gdsc_value
  )
write_csv_record(high_quality_data, "high_quality_data.csv")

fig4a_resource_overview <- data.frame(
  id = c(
    "DRIVE", "Data", "Models", "Validation",
    "GDSC", "DepMap", "Omics", "Benchmark",
    "Base models", "Ensemble", "Internal", "External"
  ),
  parent = c(
    "", "DRIVE", "DRIVE", "DRIVE",
    "Data", "Data", "Data", "Data",
    "Models", "Models", "Validation", "Validation"
  ),
  label = c(
    "DRIVE", "Data", "Models", "Validation",
    "GDSC", "DepMap", "Omics", "Benchmark",
    "Base models", "Ensemble", "Internal", "External"
  ),
  value = c(100, 40, 30, 30, 12, 10, 8, 10, 18, 12, 15, 15)
)
write_csv_record(fig4a_resource_overview, "fig4a_resource_overview.csv")

phase_levels <- c("Approved", "Phase 3", "Phase 2", "Phase 1", "Preclinical")
fig4c_comparison <- truth_grid %>%
  group_by(canonical_drug) %>%
  summarise(
    Max_Phase = phase_levels[(cur_group_id() %% length(phase_levels)) + 1],
    actual_log_ic50 = round(mean(gdsc_value), 3),
    predicted_ic50 = round(mean(gdsc_value) + rnorm(1, 0, 0.18), 3),
    .groups = "drop"
  ) %>%
  rename(drug_name = canonical_drug)
write_csv_record(fig4c_comparison, "fig4c_comparison_data.csv")

fig4e_drug_count <- expand.grid(
  OncotreeLineage = c("Large Intestine", "Pancreas", "Stomach", "Lung", "Breast"),
  Drug_Category = c("Targeted", "Cytotoxic", "Hormonal", "Epigenetic", "Metabolic")
) %>%
  mutate(
    Count = c(12, 8, 3, 4, 5, 6, 7, 2, 3, 4, 5, 2, 2, 3, 1, 4, 2, 1, 2, 1, 9, 4, 2, 3, 2),
    Count_ratio = Count / ave(Count, OncotreeLineage, FUN = sum)
  )
write_csv_record(fig4e_drug_count, "fig4e_drug_count_data.csv")

fig4f_tissues <- c(
  "Biliary tract", "Bladder", "Bone", "Brain", "Breast", "Cervix", "Esophagus",
  "Eye", "Haematopoietic and lymphoid", "Head and neck", "Kidney",
  "Large intestine", "Liver", "Lung", "Ovary", "Pancreas",
  "Peripheral nervous system", "Prostate", "Skin", "Small intestine",
  "Soft tissue", "Stomach", "Testis", "Thyroid", "Uterus", "Vulva"
)

fig4f_log2_template <- data.frame(
  tissue = fig4f_tissues,
  Terpenoids = c(18.0, 17.8, 17.1, 17.2, 16.8, 18.2, 17.7, 17.8, 16.0, 17.5, 17.6, 17.4, 18.8, 17.1, 17.3, 17.2, 17.0, 18.0, 17.5, 17.6, 17.4, 17.0, 18.9, 18.4, 18.1, 19.0),
  `Shikimates and Phenylpropanoids` = c(12.8, 13.1, 13.3, 13.1, 13.5, 12.6, 13.0, 13.0, 15.1, 12.9, 13.0, 13.1, 11.0, 13.2, 13.7, 13.8, 13.9, 12.7, 13.1, 13.2, 13.2, 13.7, 12.4, 12.3, 12.5, 12.2),
  Polyketides = c(14.0, 14.2, 14.3, 14.2, 14.1, 13.8, 14.1, 14.1, 13.0, 14.0, 14.1, 14.1, 14.0, 14.2, 14.2, 14.3, 14.2, 13.9, 14.1, 14.2, 14.1, 14.1, 13.8, 13.9, 13.8, 13.7),
  `Fatty acids` = c(10.2, 10.3, 10.4, 10.3, 10.5, 10.2, 10.3, 10.3, 12.1, 10.3, 10.2, 10.2, 10.1, 10.4, 10.5, 10.5, 10.4, 10.2, 10.3, 10.4, 10.3, 10.3, 10.1, 10.1, 10.1, 10.0),
  Carbohydrates = c(9.4, 9.4, 9.5, 9.5, 9.5, 9.3, 9.4, 9.4, 10.7, 9.4, 9.4, 9.4, 9.3, 9.5, 9.5, 9.6, 9.5, 9.4, 9.4, 9.5, 9.4, 9.4, 9.3, 9.3, 9.3, 9.3),
  `Amino acids and Peptides` = c(10.0, 10.1, 10.1, 10.1, 10.2, 9.9, 10.0, 10.0, 11.3, 10.0, 10.0, 10.0, 9.8, 10.1, 10.1, 10.2, 10.2, 10.0, 10.0, 10.1, 10.0, 10.0, 9.9, 9.9, 9.9, 9.9),
  Alkaloids = c(11.2, 11.2, 11.3, 11.2, 11.3, 11.1, 11.2, 11.2, 12.5, 11.2, 11.2, 11.2, 11.0, 11.3, 11.3, 11.4, 11.3, 11.1, 11.2, 11.2, 11.2, 11.2, 11.0, 11.0, 11.0, 11.0),
  check.names = FALSE
)

fig4f_bar_template <- data.frame(
  tissue = fig4f_tissues,
  Terpenoids = c(59, 54, 49, 50, 46, 61, 54, 54, 37, 52, 52, 51, 64, 49, 51, 48, 49, 58, 52, 54, 53, 49, 64, 60, 58, 65),
  `Shikimates and Phenylpropanoids` = c(7, 8, 9, 9, 10, 5, 8, 8, 16, 7, 8, 8, 4, 9, 9, 10, 9, 6, 8, 9, 9, 9, 4, 5, 6, 5),
  Polyketides = c(10, 13, 14, 14, 15, 12, 13, 13, 20, 13, 14, 14, 15, 13, 13, 14, 15, 12, 14, 14, 14, 14, 13, 12, 12, 12),
  `Fatty acids` = c(4, 5, 5, 5, 6, 4, 5, 5, 8, 4, 4, 4, 3, 5, 5, 5, 5, 4, 4, 5, 5, 5, 4, 4, 4, 4),
  Carbohydrates = c(8, 9, 10, 10, 11, 8, 9, 9, 8, 9, 9, 9, 6, 10, 10, 10, 10, 8, 9, 9, 9, 9, 8, 8, 8, 7),
  `Amino acids and Peptides` = c(3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2),
  Alkaloids = c(9, 8, 10, 9, 9, 8, 8, 8, 8, 12, 10, 11, 6, 11, 9, 10, 9, 9, 10, 6, 7, 11, 5, 9, 10, 5),
  check.names = FALSE
)

plot_data_for_R <- fig4f_log2_template %>%
  tidyr::pivot_longer(
    cols = -tissue,
    names_to = "pathway",
    values_to = "log2_count"
  ) %>%
  left_join(
    fig4f_bar_template %>%
      tidyr::pivot_longer(
        cols = -tissue,
        names_to = "pathway",
        values_to = "bar_weight"
      ),
    by = c("tissue", "pathway")
  ) %>%
  mutate(count = as.integer(round(2 ^ log2_count))) %>%
  select(tissue, pathway, count, bar_weight)
write_csv_record(plot_data_for_R, "plot_data_for_R.csv")

figure5_top_panel <- data.frame(
  cell_line = stripped_cells,
  ic50_lt_neg1 = c(55, 48, 43, 38, 29, 24),
  ic50_neg1_to_0 = c(20, 23, 25, 28, 31, 29),
  ic50_0_to_1 = c(15, 18, 20, 21, 24, 28),
  ic50_gt_1 = c(10, 11, 12, 13, 16, 19)
)
write_csv_record(figure5_top_panel, "figure5_top_panel.csv")

figure5_bottom_panel <- do.call(
  rbind,
  lapply(seq_along(stripped_cells), function(i) {
    x <- seq_len(200)
    data.frame(
      cell_line = stripped_cells[i],
      rank_index = x,
      ic50_value = round(sort(seq(-4.5, 7.2, length.out = 200) + rnorm(200, 0, 0.35)), 3)
    )
  })
)
write_csv_record(figure5_bottom_panel, "figure5_bottom_panel.csv")

drug_response_matrix <- round(
  matrix(
    seq(-4.8, -2.1, length.out = length(canonical_drugs) * length(stripped_cells)) +
      rnorm(length(canonical_drugs) * length(stripped_cells), 0, 0.12),
    nrow = length(canonical_drugs),
    ncol = length(stripped_cells),
    dimnames = list(canonical_drugs, stripped_cells)
  ),
  3
)
write_csv_record(as.data.frame.matrix(drug_response_matrix), "drug_response_matrix.csv", row.names = TRUE)

ensure_dir("full_prediction_results_cellline")
ensure_dir("full_prediction_results_tcga")

top_freq_drugs <- c(
  canonical_drugs,
  "TRAMETINIB", "DASATINIB", "BORTEZOMIB", "ERLOTINIB", "OXALIPLATIN",
  "PALBOCICLIB", "VINCRISTINE", "GEMCITABINE", "PACLITAXEL", "IRINOTECAN"
)

for (i in seq_len(8)) {
  cl_df <- data.frame(
    drug_name = top_freq_drugs,
    Ensemble_Score = round(seq(-2.8, 0.6, length.out = length(top_freq_drugs)) + rnorm(length(top_freq_drugs), 0, 0.22), 3)
  )
  tcga_df <- data.frame(
    drug_name = rev(top_freq_drugs),
    Ensemble_Score = round(seq(-2.6, 0.8, length.out = length(top_freq_drugs)) + rnorm(length(top_freq_drugs), 0, 0.25), 3)
  )
  write.csv(cl_df, file.path("full_prediction_results_cellline", sprintf("sample_%02d_full_predictions.csv", i)), row.names = FALSE)
  write.csv(tcga_df, file.path("full_prediction_results_tcga", sprintf("sample_%02d_full_predictions.csv", i)), row.names = FALSE)
}
record_overview("full_prediction_results_cellline/*.csv", data.frame(rows = 8, columns = 2))
record_overview("full_prediction_results_tcga/*.csv", data.frame(rows = 8, columns = 2))

gene_expression <- matrix(
  round(rexp(length(stripped_cells) * 18, rate = 0.08), 3),
  nrow = length(stripped_cells),
  dimnames = list(stripped_cells, paste0("Gene", seq_len(18)))
)
cnv_matrix <- matrix(
  round(rnorm(length(stripped_cells) * 12, 0, 1.2), 3),
  nrow = length(stripped_cells),
  dimnames = list(stripped_cells, paste0("CNV", seq_len(12)))
)
mutation_matrix <- matrix(
  sample(0:1, length(stripped_cells) * 14, replace = TRUE, prob = c(0.7, 0.3)),
  nrow = length(stripped_cells),
  dimnames = list(stripped_cells, paste0("MUT", seq_len(14)))
)
write_csv_record(as.data.frame.matrix(gene_expression), "gene_expression.csv", row.names = TRUE)
write_csv_record(as.data.frame.matrix(cnv_matrix), "cnv.csv", row.names = TRUE)
write_csv_record(as.data.frame.matrix(mutation_matrix), "mutation.csv", row.names = TRUE)

gene_raw <- as.data.frame(t(gene_expression), check.names = FALSE)
write_csv_record(gene_raw, "gene_raw.csv", row.names = TRUE)

gmt_lines <- c(
  paste("KEGG_CELL_CYCLE", "Synthetic", paste(colnames(gene_expression)[1:6], collapse = "\t"), sep = "\t"),
  paste("KEGG_DNA_REPAIR", "Synthetic", paste(colnames(gene_expression)[7:12], collapse = "\t"), sep = "\t"),
  paste("KEGG_MAPK_SIGNALING", "Synthetic", paste(colnames(gene_expression)[13:18], collapse = "\t"), sep = "\t")
)
writeLines(gmt_lines, "KEGG_human_latest_Symbol.gmt")
record_overview("KEGG_human_latest_Symbol.gmt", data.frame(rows = 3, columns = 8))

response_df <- data.frame(matrix("", nrow = nrow(truth_grid), ncol = 17))
colnames(response_df) <- paste0("col", seq_len(17))
response_df$col2 <- truth_grid$cell_line
response_df$col11 <- round(exp(truth_grid$gdsc_value + 2.2), 4)
response_df$col12 <- paste0("DRUG", sprintf("%03d", truth_grid$drug_idx))
response_df$col17 <- truth_grid$canonical_drug
write_csv_record(response_df, "secondary-screen-dose-response-curve-parameters.csv")

cellnames_df <- data.frame(ID = stripped_cells[1:4])
write_csv_record(cellnames_df, "cellnames.csv")

supfig_s1_model_summary <- data.frame(
  Model = methods,
  Omics_Feature = c(
    "mRNA expression;Gene mutation",
    "mRNA expression;Copy number variation",
    "mRNA expression;Gene mutation;Copy number variation",
    "mRNA expression",
    "mRNA expression;Protein-protein interaction",
    "mRNA expression;Gene mutation",
    "mRNA expression",
    "mRNA expression;miRNA expression",
    "mRNA expression",
    "mRNA expression",
    "mRNA expression;DNA methylation",
    "mRNA expression;Gene mutation"
  ),
  Drug_Representation = c(
    "Molecular graph;Fingerprint",
    "Molecular graph;ESPF",
    "Molecular graph;Morgan",
    "ESPF;SMILESVec",
    "Molecular graph;Fingerprint",
    "Molecular graph",
    "Molecular graph;ECFP",
    "SMILESVec;Fingerprint",
    "ESPF",
    "Molecular graph",
    "Molecular graph;Morgan",
    "Fingerprint;ECFP"
  ),
  Algorithm_Type = c(
    "GCN;Transformer",
    "AE;DNN",
    "GCN;MLP",
    "Transformer",
    "GAT;DNN",
    "GIN",
    "GAT",
    "DNN",
    "Transformer",
    "GCN",
    "MCA;DNN",
    "AE;MLP"
  )
)
write_csv_record(supfig_s1_model_summary, "supfig_s1_model_summary.csv")

performance_data <- data.frame(
  Cell_N = c(1, 2, 5, 10, 50, 100, 500, 1000),
  Time = c(0.4, 0.7, 1.3, 2.4, 8.8, 16.2, 64.5, 122.0),
  Peak_RAM = c(0.15, 0.20, 0.28, 0.41, 0.95, 1.55, 4.80, 8.90)
)
write_csv_record(performance_data, "performance_data.csv")

predicted_ic50_by_tissue <- do.call(
  rbind,
  lapply(c("Colon", "Rectum", "Pancreas", "Stomach", "Lung", "Breast"), function(tissue) {
    data.frame(
      Tissue = tissue,
      IC50 = round(rnorm(40, mean = runif(1, -1.5, 2.5), sd = runif(1, 0.6, 1.2)), 3)
    )
  })
)
write_csv_record(predicted_ic50_by_tissue, "predicted_ic50_by_tissue.csv")

method_order_full <- c(
  "ML baseline",
  "BANDRP",
  "DeepAEG",
  "DeepCCDS",
  "DeepCDR",
  "DeepTTA",
  "DIPK",
  "GADRP",
  "GPDRP_GAT",
  "GPDRP_GCN",
  "GPDRP_GIN",
  "GPDRP_Trans",
  "GraphDRP-GINConvNet",
  "GraphDRP_GAT_GCN",
  "GraphDRP_GATNet",
  "GraphDRP_GCNNet",
  "NERD",
  "paccmann",
  "Precily"
)

mode_offsets <- c("Mixed" = 0.00, "Cell-blind" = 0.08, "Drug-blind" = 0.05)
model_metric <- do.call(
  rbind,
  lapply(seq_along(method_order_full), function(i) {
    do.call(
      rbind,
      lapply(names(mode_offsets), function(mode_name) {
        do.call(
          rbind,
          lapply(1:5, function(rep_id) {
            base_shift <- (i - 1) * 0.008 + mode_offsets[[mode_name]]
            data.frame(
              Method = method_order_full[i],
              mode = mode_name,
              RMSE = round(0.68 + base_shift + rnorm(1, 0, 0.012), 3),
              PCC = round(0.83 - base_shift + rnorm(1, 0, 0.010), 3),
              SCC = round(0.79 - base_shift + rnorm(1, 0, 0.010), 3),
              R2 = round(0.70 - base_shift + rnorm(1, 0, 0.012), 3),
              NDCG = round(0.85 - base_shift / 1.2 + rnorm(1, 0, 0.008), 3),
              NWPC = round(0.82 - base_shift / 1.2 + rnorm(1, 0, 0.008), 3)
            )
          })
        )
      })
    )
  })
)
write_csv_record(model_metric, "model_metric.csv")

all_split_metrics <- bind_rows(
  data.frame(
    strategy = "random_5fold_cv",
    split = paste0("Fold", 1:5),
    holdout_group = NA_character_,
    internal_RMSE = c(0.71, 0.70, 0.69, 0.72, 0.70),
    internal_PCC = c(0.82, 0.83, 0.84, 0.81, 0.83),
    external_RMSE = c(0.79, 0.80, 0.78, 0.81, 0.80),
    external_PCC = c(0.74, 0.75, 0.76, 0.73, 0.75),
    n_test_cells = 120,
    external_n_cells = 95,
    n_test_drugs = 30,
    external_n_drugs = 24
  ),
  data.frame(
    strategy = "leave_tissue_organ_out",
    split = c("blood", "lung", "urogenital_system"),
    holdout_group = c("blood", "lung", "urogenital_system"),
    internal_RMSE = c(0.73, 0.76, 0.78),
    internal_PCC = c(0.79, 0.76, 0.74),
    external_RMSE = c(0.82, 0.86, 0.88),
    external_PCC = c(0.69, 0.65, 0.62),
    n_test_cells = c(48, 61, 37),
    external_n_cells = c(32, 40, 25),
    n_test_drugs = c(NA, NA, NA),
    external_n_drugs = c(NA, NA, NA)
  ),
  data.frame(
    strategy = "leave_chemical_cluster_out",
    split = c("Chem-A", "Chem-B", "Chem-C"),
    holdout_group = c("Chem-A", "Chem-B", "Chem-C"),
    internal_RMSE = c(0.75, 0.77, 0.79),
    internal_PCC = c(0.77, 0.75, 0.73),
    external_RMSE = c(0.84, 0.87, 0.90),
    external_PCC = c(0.67, 0.63, 0.60),
    n_test_cells = c(NA, NA, NA),
    external_n_cells = c(NA, NA, NA),
    n_test_drugs = c(14, 16, 18),
    external_n_drugs = c(10, 11, 12)
  )
)
write_csv_record(all_split_metrics, "ALL_split_level_internal_external_metrics.csv")

c3_curves <- expand.grid(
  strategy = c("RandomForest", "RidgeCV", "LightGBM", "KNN", "XGBoost"),
  num_models = 2:10
) %>%
  mutate(
    strategy_offset = c(0.00, 0.03, 0.02, 0.05, 0.01)[match(strategy, c("RandomForest", "RidgeCV", "LightGBM", "KNN", "XGBoost"))],
    train_RMSE = round(0.76 - 0.018 * num_models + strategy_offset + rnorm(n(), 0, 0.01), 3),
    validation_RMSE = round(0.88 - 0.013 * num_models + strategy_offset + rnorm(n(), 0, 0.012), 3)
  ) %>%
  select(strategy, num_models, train_RMSE, validation_RMSE)
write_csv_record(c3_curves, "C3_loo_18_to_2_curves_all_strategies.csv")

drug_groups <- c("Alisporivir", "Josamycin", "NIM811", "Rifabutin", "Cephalomannine", "Belotecan")
drug_response_ic50 <- do.call(
  rbind,
  lapply(c("HCT116", "RKO", "Lovo"), function(cell_name) {
    do.call(
      rbind,
      lapply(seq_along(drug_groups), function(i) {
        ref_value <- round(exp(0.9 + i * 0.08 + runif(1, 0, 0.12)), 3)
        test_value <- round(exp(0.1 + i * 0.12 + runif(1, 0, 0.18)), 3)
        rbind(
          data.frame(CellLine = cell_name, Drug_Group = drug_groups[i], Type = "Test", Compound_Name = drug_groups[i], IC50 = test_value),
          data.frame(CellLine = cell_name, Drug_Group = drug_groups[i], Type = "Ref", Compound_Name = "5-Fu", IC50 = ref_value)
        )
      })
    )
  })
)
write_csv_record(drug_response_ic50, "drug_response_ic50.csv")

dataset_overview <- bind_rows(overview_entries) %>%
  distinct(file, .keep_all = TRUE)
write.csv(dataset_overview, "synthetic_test_data_overview.csv", row.names = FALSE)

message("Synthetic test data generated successfully.")
