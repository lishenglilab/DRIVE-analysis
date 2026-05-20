############################################################
## ic50_matrix_process.R
##
## Purpose:
##   Convert long-format drug response data into an IC50 matrix.
##
## Input:
##   1. secondary-screen-dose-response-curve-parameters.csv
##      A long-format response table containing cell-line IDs,
##      IC50 values, drug IDs, and drug names.
##
## Optional input:
##   1. cellnames.csv
##      A cell-line list used to filter the IC50 matrix.
##
## Output:
##   1. predict_drug.csv
##      Unique drug ID and drug name table.
##   2. celllines_all.csv
##      Unique cell-line list.
##   3. ic50_matrix_raw.csv
##      Raw IC50 matrix, with cell lines as rows and drugs as columns.
##   4. ic50_matrix_log.csv
##      Natural log-transformed IC50 matrix.
##   5. ic50_matrix_log_selected_cells.csv
##      Optional output if cellnames.csv is provided.
############################################################

############################################################
## 1. Set working directory and input files
############################################################

## Example path. Replace this with your own project directory.
work_dir <- "path/to/your/project"
setwd(work_dir)

response_file <- "secondary-screen-dose-response-curve-parameters.csv"

## Optional file. If this file exists, the script will filter the IC50 matrix
## to the selected cell lines.
selected_cell_file <- "cellnames.csv"

out_dir <- "ic50_matrix_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Set column indexes
############################################################

## These indexes follow the original processing script:
## response <- response[, c(2, 11, 12, 17)]
##
## After extraction:
##   column 1: cell-line identifier
##   column 2: IC50 value
##   column 3: drug identifier
##   column 4: drug name or drug annotation

cell_col_index <- 2
ic50_col_index <- 11
drug_id_col_index <- 12
drug_name_col_index <- 17

############################################################
## 3. Read response data
############################################################

response_raw <- read.csv(
  response_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

response <- response_raw[, c(
  cell_col_index,
  ic50_col_index,
  drug_id_col_index,
  drug_name_col_index
)]

colnames(response) <- c(
  "cell_line",
  "IC50",
  "Drug.Id",
  "Drug.Name"
)

response$IC50 <- as.numeric(response$IC50)

############################################################
## 4. Remove invalid records
############################################################

response <- response[
  !is.na(response$cell_line) &
    response$cell_line != "" &
    !is.na(response$Drug.Id) &
    response$Drug.Id != "" &
    !is.na(response$IC50),
]

response <- response[response$IC50 > 0, ]

message("Valid response records: ", nrow(response))

############################################################
## 5. Export unique drug and cell-line information
############################################################

drug_info <- unique(response[, c("Drug.Id", "Drug.Name")])
rownames(drug_info) <- NULL

cell_line_info <- data.frame(
  cell_line = unique(response$cell_line),
  stringsAsFactors = FALSE
)

write.csv(
  drug_info,
  file = file.path(out_dir, "predict_drug.csv"),
  row.names = FALSE
)

write.csv(
  cell_line_info,
  file = file.path(out_dir, "celllines_all.csv"),
  row.names = FALSE
)

############################################################
## 6. Aggregate duplicated cell-line and drug pairs
############################################################

## If the same cell-line/drug pair appears more than once,
## the mean IC50 value is used.

response_agg <- aggregate(
  IC50 ~ cell_line + Drug.Id,
  data = response,
  FUN = mean,
  na.rm = TRUE
)

############################################################
## 7. Convert long-format response data to IC50 matrix
############################################################

unique_cell_lines <- unique(response_agg$cell_line)
unique_drugs <- unique(response_agg$Drug.Id)

ic50_matrix <- matrix(
  NA,
  nrow = length(unique_cell_lines),
  ncol = length(unique_drugs),
  dimnames = list(unique_cell_lines, unique_drugs)
)

for (i in seq_len(nrow(response_agg))) {
  ic50_matrix[
    response_agg$cell_line[i],
    response_agg$Drug.Id[i]
  ] <- response_agg$IC50[i]
}

############################################################
## 8. Log-transform IC50 matrix
############################################################

## Natural log transformation, consistent with log(IC50)
## used in the original processing script.

log_ic50_matrix <- log(ic50_matrix)

############################################################
## 9. Export IC50 matrices
############################################################

write.csv(
  ic50_matrix,
  file = file.path(out_dir, "ic50_matrix_raw.csv"),
  row.names = TRUE
)

write.csv(
  log_ic50_matrix,
  file = file.path(out_dir, "ic50_matrix_log.csv"),
  row.names = TRUE
)

############################################################
## 10. Optional: filter selected cell lines
############################################################

if (file.exists(selected_cell_file)) {
  selected_cells <- read.csv(
    selected_cell_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  ## The selected cell-line file is expected to contain one of these columns.
  selected_col <- intersect(
    c("ID", "cell_line", "DepMap_ID", "Cell_line"),
    colnames(selected_cells)
  )[1]

  if (!is.na(selected_col)) {
    selected_cell_ids <- selected_cells[[selected_col]]

    log_ic50_selected <- log_ic50_matrix[
      rownames(log_ic50_matrix) %in% selected_cell_ids,
      ,
      drop = FALSE
    ]

    write.csv(
      log_ic50_selected,
      file = file.path(out_dir, "ic50_matrix_log_selected_cells.csv"),
      row.names = TRUE
    )

    message("Selected cell-line matrix exported: ",
            nrow(log_ic50_selected), " cell lines x ",
            ncol(log_ic50_selected), " drugs")
  } else {
    warning("cellnames.csv was found, but no valid cell-line ID column was detected.")
  }
}

############################################################
## Finished
############################################################

message("IC50 matrix processing finished.")
message("Raw IC50 matrix dimension: ",
        nrow(ic50_matrix), " cell lines x ", ncol(ic50_matrix), " drugs")
message("Output directory: ", normalizePath(out_dir))
