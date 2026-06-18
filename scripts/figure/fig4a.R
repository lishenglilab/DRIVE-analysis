############################################################
## fig4a_sunburst_resource_overview.R
##
## Purpose:
##   Generate the Fig. 4A DRIVE resource overview sunburst plot.
##
## Input:
##   fig4a_resource_overview.csv
##
## Required columns:
##   1. id     : unique node identifier
##   2. parent : parent node identifier; use empty value for the root
##   3. label  : display label
##   4. value  : numeric node value
##
## Output:
##   1. fig4a_sunburst/sunburst_plot_validate.html
##   2. fig4a_sunburst/sunburst_plot_validate.pdf
############################################################

suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
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

input_file <- "fig4a_resource_overview.csv"
out_dir <- "fig4a_sunburst"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

df <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("id", "parent", "label", "value")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df$parent[is.na(df$parent)] <- ""
df$value <- as.numeric(df$value)

p <- plot_ly(
  data = df,
  ids = ~id,
  labels = ~label,
  parents = ~parent,
  values = ~value,
  type = "sunburst",
  branchvalues = "total",
  textinfo = "label+value",
  insidetextorientation = "radial"
) %>%
  layout(
    title = list(text = "DRIVE resource overview", x = 0.5),
    margin = list(l = 0, r = 0, b = 0, t = 50)
  )

html_file <- file.path(out_dir, "sunburst_plot_validate.html")
pdf_file <- file.path(out_dir, "sunburst_plot_validate.pdf")

saveWidget(p, html_file, selfcontained = FALSE)

if (requireNamespace("webshot2", quietly = TRUE)) {
  webshot2::webshot(html_file, pdf_file, vwidth = 1200, vheight = 900)
} else {
  message("webshot2 is not installed. HTML was exported; install webshot2 to export PDF.")
}

message("Fig. 4A sunburst plot finished.")
