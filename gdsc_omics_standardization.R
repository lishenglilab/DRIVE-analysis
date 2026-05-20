############################################################
## gdsc_omics_standardization.R
##
## Purpose:
##   Standardize omics feature matrices used for downstream
##   drug response prediction analysis.
##
## Input:
##   1. gene_expression.csv  : cell lines x genes
##   2. cnv.csv              : cell lines x CNV features
##   3. mutation.csv         : cell lines x mutation features
##
## Output:
##   1. gene_z.csv           : z-score standardized expression matrix
##   2. cnv_z.csv            : z-score standardized CNV matrix
##   3. mutation_numeric.csv : numeric mutation matrix
############################################################

############################################################
## 1. Set working directory and input files
############################################################

## Example path. Replace this with your own project directory.
work_dir <- "path/to/your/project"
setwd(work_dir)

expr_file <- "gene_expression.csv"
cnv_file  <- "cnv.csv"
mut_file  <- "mutation.csv"

out_dir <- "omics_standardization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Helper functions
############################################################

read_omics_matrix <- function(file) {
  mat <- read.csv(
    file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  mat <- as.data.frame(mat, check.names = FALSE)
  mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
  mat <- as.matrix(mat)

  return(mat)
}

clean_matrix <- function(mat) {
  ## Remove samples or features with missing names
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- mat[, !is.na(colnames(mat)) & colnames(mat) != "", drop = FALSE]

  ## Remove duplicated samples and features
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]

  return(mat)
}

zscore_by_feature <- function(mat) {
  ## Z-score standardization is performed for each feature
  ## across all cell lines.
  mat_z <- scale(mat)

  ## Features with zero variance will generate NA values after scaling.
  ## These values are replaced with 0.
  mat_z[is.na(mat_z)] <- 0

  mat_z <- as.matrix(mat_z)
  rownames(mat_z) <- rownames(mat)
  colnames(mat_z) <- colnames(mat)

  return(mat_z)
}

############################################################
## 3. Gene expression standardization
############################################################

if (file.exists(expr_file)) {
  gene_expr <- read_omics_matrix(expr_file)
  gene_expr <- clean_matrix(gene_expr)

  ## Log2 transformation for expression data
  gene_expr <- log2(gene_expr + 1)

  ## Z-score standardization by gene
  gene_z <- zscore_by_feature(gene_expr)

  write.csv(
    gene_z,
    file = file.path(out_dir, "gene_z.csv"),
    row.names = TRUE
  )

  message("Gene expression matrix standardized: ",
          nrow(gene_z), " cell lines x ", ncol(gene_z), " genes")
}

############################################################
## 4. CNV standardization
############################################################

if (file.exists(cnv_file)) {
  cnv <- read_omics_matrix(cnv_file)
  cnv <- clean_matrix(cnv)

  ## Z-score standardization by CNV feature
  cnv_z <- zscore_by_feature(cnv)

  write.csv(
    cnv_z,
    file = file.path(out_dir, "cnv_z.csv"),
    row.names = TRUE
  )

  message("CNV matrix standardized: ",
          nrow(cnv_z), " cell lines x ", ncol(cnv_z), " features")
}

############################################################
## 5. Mutation feature processing
############################################################

if (file.exists(mut_file)) {
  mutation <- read_omics_matrix(mut_file)
  mutation <- clean_matrix(mutation)

  ## Mutation features are usually binary or discrete.
  ## Here they are converted to numeric values and kept without z-score scaling.
  mutation[is.na(mutation)] <- 0

  write.csv(
    mutation,
    file = file.path(out_dir, "mutation_numeric.csv"),
    row.names = TRUE
  )

  message("Mutation matrix processed: ",
          nrow(mutation), " cell lines x ", ncol(mutation), " features")
}

############################################################
## Finished
############################################################

message("Omics standardization finished.")
message("Output directory: ", normalizePath(out_dir))
