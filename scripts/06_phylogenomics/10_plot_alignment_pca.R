#!/usr/bin/env Rscript
# --- Parse PCA table and figure paths ---
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) stop("Missing ", flag); args[[i+1]] }
scores <- read.delim(value("--scores"), check.names = FALSE)
variance <- read.delim(value("--variance"), check.names = FALSE)
out <- value("--output-dir"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
cols <- c(Reproductive = "#0072B2", Bacteraemia = "#D55E00")
# --- Validate and prepare the three requested label styles ---
required <- c("assembly_accession", "strain", "reproductive_bacteraemia_category")
if (!all(required %in% names(scores))) {
  stop("PCA scores lack required columns: ", paste(setdiff(required, names(scores)), collapse = ", "))
}
scores$strain[is.na(scores$strain) | scores$strain == ""] <-
  scores$assembly_accession[is.na(scores$strain) | scores$strain == ""]
label_sets <- list(
  strain = scores$strain,
  assembly_accession = scores$assembly_accession
)

# --- Spread dense labels and connect them to their points ---
spread_labels <- function(x, y, labels, cex) {
  position_x <- x
  position_y <- y + 1.2 * strheight(labels, cex = cex)
  width <- strwidth(labels, cex = cex)
  height <- strheight(labels, cex = cex)
  limits <- par("usr")
  for (iteration in seq_len(500)) {
    for (i in seq_along(labels)) {
      for (j in seq_len(i - 1L)) {
        overlap_x <- abs(position_x[i] - position_x[j]) < (width[i] + width[j]) / 2
        overlap_y <- abs(position_y[i] - position_y[j]) < (height[i] + height[j]) * 0.7
        if (overlap_x && overlap_y) {
          direction <- if (position_y[i] >= position_y[j]) 1 else -1
          if (position_y[i] == position_y[j]) direction <- if (i %% 2L) 1 else -1
          push <- (height[i] + height[j]) * 0.38
          position_y[i] <- position_y[i] + direction * push
          position_y[j] <- position_y[j] - direction * push
        }
      }
    }
    position_x <- position_x + 0.003 * (x - position_x)
    position_y <- position_y + 0.003 * (y - position_y)
    position_x <- pmin(pmax(position_x, limits[1] + width / 2), limits[2] - width / 2)
    position_y <- pmin(pmax(position_y, limits[3] + height), limits[4] - height)
  }
  list(x = position_x, y = position_y)
}

# --- Draw one labelled PCA view with a bottom-right legend ---
draw_plot <- function(x, y, labels, label_type) {
  group <- scores$reproductive_bacteraemia_category
  x_label <- sprintf("%s (%.2f%%)", x, 100 * variance$explained_variance_proportion[variance$component == x])
  y_label <- sprintf("%s (%.2f%%)", y, 100 * variance$explained_variance_proportion[variance$component == y])
  par(mar = c(5, 5, 1, 1))
  x_limits <- extendrange(scores[[x]], f = 0.12)
  y_limits <- extendrange(scores[[y]], f = 0.15)
  plot(scores[[x]], scores[[y]], pch = 19, col = cols[group], xlab = x_label,
       ylab = y_label, xlim = x_limits, ylim = y_limits)
  label_cex <- c(strain = 0.42, assembly_accession = 0.32)[[label_type]]
  positions <- spread_labels(scores[[x]], scores[[y]], labels, label_cex)
  segments(scores[[x]], scores[[y]], positions$x, positions$y, col = "grey70", lwd = 0.5)
  text(positions$x, positions$y, labels = labels, cex = label_cex)
  legend("bottomright", names(cols), col = cols, pch = 19, bty = "n")
}

# --- Save PNG and PDF versions for every component and label style ---
plot_pair <- function(x, y, stem) {
  for (label_type in names(label_sets)) {
    suffix <- if (label_type == "strain") "" else paste0("_", label_type, "_labels")
    output_stem <- paste0(stem, suffix)
    png(file.path(out, paste0(output_stem, ".png")), width = 2400, height = 1800, res = 200)
    draw_plot(x, y, label_sets[[label_type]], label_type)
    dev.off()
    pdf(file.path(out, paste0(output_stem, ".pdf")), width = 12, height = 9)
    draw_plot(x, y, label_sets[[label_type]], label_type)
    dev.off()
  }
}
plot_pair("PC1", "PC2", "pca_pc1_pc2")
plot_pair("PC1", "PC3", "pca_pc1_pc3")
message("PCA plots written for ", nrow(scores), " genomes using separate strain and assembly labels")
