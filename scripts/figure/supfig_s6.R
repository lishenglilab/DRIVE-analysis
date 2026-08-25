############################################################
## figure_s6_cv_base_r.R
##
## Purpose:
##   Generate Supplementary Fig. S6 showing coefficient of
##   variation (CV) across metrics and validation modes.
##
## Input:
##   fig2.csv
##
## Required columns:
##   1. Method
##   2. Setting
##   3. Metric
##   4. Value
##
## Output:
##   1. Figure_S6_CV_plotted_points.csv
##   2. Figure_S6_raw_input_values.csv
##   3. Supplementary_Figure_S6_CV_recreated_final3.pdf
############################################################

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1) args[[1]] else "C:/Users/26336/Documents/Codex/2026-08-18/new-chat/work/fig2_cv/fig2.csv"
out_dir <- if (length(args) >= 2) args[[2]] else "C:/Users/26336/Documents/Codex/2026-08-18/new-chat/output/fig2_cv"
if (!file.exists(input_csv)) {
  stop("Input file not found: ", input_csv)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## 1. Read data
## -------------------------------------------------------------------------
raw <- read.csv(input_csv, check.names = FALSE)
raw$Value <- as.numeric(raw$Value)

required_cols <- c("Method", "Setting", "Metric", "Value")
missing_cols <- setdiff(required_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

## -------------------------------------------------------------------------
## 2. Define factor ordering
## -------------------------------------------------------------------------
metric_order <- c("RMSE", "PCC", "SCC", "R2", "NDCG", "NWPC")
setting_order <- c("Mixed", "Cell-blind", "Drug-blind")
raw$Metric <- factor(raw$Metric, levels = metric_order)
raw$Setting <- factor(raw$Setting, levels = setting_order)

## -------------------------------------------------------------------------
## 3. Calculate CV for each method, setting, and metric
## -------------------------------------------------------------------------
# Each plotted point is the coefficient of variation for one
# method/setting/metric combination.
keys <- interaction(raw$Method, raw$Setting, raw$Metric, drop = TRUE, sep = "\u001f")
split_raw <- split(raw, keys)
cv_list <- lapply(split_raw, function(z) {
  m <- mean(z$Value, na.rm = TRUE)
  s <- sd(z$Value, na.rm = TRUE)
  data.frame(Method = as.character(z$Method[1]), Setting = as.character(z$Setting[1]),
             Metric = as.character(z$Metric[1]), N_replicates = sum(!is.na(z$Value)),
             Mean = m, SD = s, CV = if (is.finite(m) && m != 0) s / m else NA_real_)
})
cv <- do.call(rbind, cv_list)
cv$Metric <- factor(cv$Metric, levels = metric_order)
cv$Setting <- factor(cv$Setting, levels = setting_order)
cv <- cv[order(cv$Metric, cv$Setting, cv$Method), ]
rownames(cv) <- NULL

write.csv(cv, file.path(out_dir, "Figure_S6_CV_plotted_points.csv"), row.names = FALSE, na = "")
write.csv(raw, file.path(out_dir, "Figure_S6_raw_input_values.csv"), row.names = FALSE, na = "")

## -------------------------------------------------------------------------
## 4. Define plotting parameters
## -------------------------------------------------------------------------
cols <- c("Mixed" = "#E7B83C", "Cell-blind" = "#4C91B9", "Drug-blind" = "#37A878")

## -------------------------------------------------------------------------
## 5. Generate figure
## -------------------------------------------------------------------------
pdf(
  file.path(out_dir, "Supplementary_Figure_S6_CV_recreated_final3.pdf"),
  width = 10.5,
  height = 6.2,
  useDingbats = FALSE
)

par(
  mar = c(5.3, 7.0, 2.0, 1.2),
  mgp = c(3.0, 0.8, 0),
  family = "sans",
  las = 1,
  xaxs = "i",
  yaxs = "i"
)
plot.new()
usr <- c(0.5, 6.5, -3.1, 2.6)
plot.window(xlim = usr[1:2], ylim = usr[3:4])
line_width <- 1.15
axis(
  1,
  at = 1:6,
  labels = metric_order,
  cex.axis = 0.95,
  tck = -0.018,
  lwd = line_width,
  lwd.ticks = line_width
)
axis(
  2,
  at = -3:2,
  labels = -3:2,
  cex.axis = 0.95,
  tck = -0.018,
  lwd = line_width,
  lwd.ticks = line_width
)
box(lwd = line_width)
abline(h = 0, col = "#666666", lty = 2, lwd = line_width)
mtext("Coefficient of Variation (CV)", side = 2, line = 3.0, las = 0, cex = 1.05)

set.seed(20260825)
offsets <- c("Mixed" = -0.25, "Cell-blind" = 0, "Drug-blind" = 0.25)
for (mi in seq_along(metric_order)) {
  met <- metric_order[mi]
  vals_list <- lapply(setting_order, function(st) {
    z <- cv$CV[cv$Metric == met & cv$Setting == st]
    z[is.finite(z)]
  })
  boxplot(
    vals_list,
    at = mi + unname(offsets[setting_order]),
    add = TRUE,
    width = rep(0.01, 3),
    axes = FALSE,
    outline = TRUE,
    outpch = 16,
    outcex = 0.70,
    outcol = adjustcolor("#4A4A4A", alpha.f = 0.55),
    col = unname(adjustcolor(cols[setting_order], alpha.f = 0.88)),
    border = "#4D4D4D",
    pars = list(
      boxlwd = line_width,
      medlwd = line_width,
      whisklwd = line_width,
      staplelwd = line_width
    ),
    staplewex = 0.55
  )
  for (si in seq_along(setting_order)) {
    vals <- vals_list[[si]]
    x <- mi + offsets[[setting_order[si]]]
    if (length(vals) > 0) {
      if (length(vals) == 1) {
        jitter_x <- x
      } else {
        jitter_x <- x + seq(-0.035, 0.035, length.out = length(vals))
      }
      points(
        jitter_x,
        vals,
        pch = 16,
        cex = 0.70,
        col = adjustcolor("#4A4A4A", alpha.f = 0.55)
      )
    }
  }
}
## -------------------------------------------------------------------------
## 6. Add legend and save
## -------------------------------------------------------------------------
legend("top", inset = c(0, -0.08), xpd = NA, horiz = TRUE, bty = "n",
       legend = setting_order, fill = unname(cols[setting_order]), border = "#4D4D4D",
       pt.cex = 0.8, cex = 0.9)
dev.off()

cat("Wrote", nrow(cv), "CV points and", nrow(raw), "raw values to", out_dir, "\n")
