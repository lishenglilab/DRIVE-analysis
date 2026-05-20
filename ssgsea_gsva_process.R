############################################################
## GSVA and ssGSEA pathway activity score calculation
##
## Input:
##   1. gene expression matrix: genes x samples/cell lines
##   2. pathway gene set file: GMT format
##
## Output:
##   1. gsva_scores.csv
##   2. ssgsea_scores.csv
############################################################

library(GSVA)
library(GSEABase)

############################################################
## 1. Set working directory
############################################################

work_dir <- "path/to/your/project"
setwd(work_dir)

############################################################
## 2. Load expression matrix
############################################################

# The expression matrix should use genes as rows and cell lines/samples as columns.
# Row names should be gene symbols.
expr <- read.csv(
  "gene_raw.csv",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(expr)

# Convert all values to numeric
expr <- apply(expr, 2, as.numeric)
rownames(expr) <- rownames(read.csv(
  "gene_raw.csv",
  row.names = 1,
  check.names = FALSE
))

############################################################
## 3. Basic preprocessing
############################################################

# Remove genes with missing gene names
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[rownames(expr) != "", ]

# Remove duplicated gene symbols
expr <- expr[!duplicated(rownames(expr)), ]

# Remove genes with constant expression across all samples
non_constant_genes <- apply(expr, 1, function(x) length(unique(x)) > 1)
expr <- expr[non_constant_genes, ]

# Log2 transformation
expr <- log2(expr + 1)

############################################################
## 4. Load pathway gene sets
############################################################

gene_sets <- getGmt("KEGG_human_latest_Symbol.gmt")

# Alternative:
# gene_sets <- getGmt("c2.cp.v6.1.symbols.gmt")

############################################################
## 5. Calculate GSVA scores
############################################################

gsva_scores <- gsva(
  expr = expr,
  gset.idx.list = gene_sets,
  method = "gsva",
  kcdf = "Gaussian",
  abs.ranking = FALSE,
  min.sz = 5
)

############################################################
## 6. Calculate ssGSEA scores
############################################################

ssgsea_scores <- gsva(
  expr = expr,
  gset.idx.list = gene_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = FALSE,
  min.sz = 5
)

############################################################
## 7. Export results
############################################################

# Output format 1: pathways x samples
write.csv(
  gsva_scores,
  "gsva_scores_pathway_by_sample.csv"
)

write.csv(
  ssgsea_scores,
  "ssgsea_scores_pathway_by_sample.csv"
)

# Output format 2: samples x pathways
gsva_scores_t <- t(gsva_scores)
ssgsea_scores_t <- t(ssgsea_scores)

write.csv(
  gsva_scores_t,
  "gsva_scores_sample_by_pathway.csv"
)

write.csv(
  ssgsea_scores_t,
  "ssgsea_scores_sample_by_pathway.csv"
)

############################################################
## Finished
############################################################

message("GSVA and ssGSEA calculation finished.")