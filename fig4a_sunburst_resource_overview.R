############################################################
## fig4a_sunburst_resource_overview.R
## Output: sunburst_plot_validate.pdf
############################################################

suppressPackageStartupMessages({
  library(plotly)
  library(htmlwidgets)
})

work_dir <- "path/to/your/project"
setwd(work_dir)

input_file <- "fig4a_resource_overview.csv"
out_dir <- "fig4a_sunburst"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Required input columns:
## id,parent,label,value
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

saveWidget(p, html_file, selfcontained = TRUE)

if (requireNamespace("webshot2", quietly = TRUE)) {
  webshot2::webshot(html_file, pdf_file, vwidth = 1200, vheight = 900)
} else {
  message("webshot2 is not installed. HTML was exported; install webshot2 to export PDF.")
}

message("Fig. 4A sunburst plot finished.")
