#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}

scores_path <- value("--scores")
variance_path <- value("--variance")
out <- value("--output-dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

scores <- read.delim(scores_path, check.names = FALSE, stringsAsFactors = FALSE)
variance <- read.delim(variance_path, check.names = FALSE, stringsAsFactors = FALSE)

required <- c("assembly_accession", "strain", "reproductive_bacteraemia_category")
if (!all(required %in% names(scores))) {
  stop("PCA scores lack required columns: ", paste(setdiff(required, names(scores)), collapse = ", "))
}

scores$strain <- ifelse(is.na(scores$strain) | scores$strain == "", scores$assembly_accession, scores$strain)

source_colours <- c(Reproductive = "#0072B2", Bacteraemia = "#D55E00")
source_shapes <- c(Reproductive = 17L, Bacteraemia = 16L)

# Keep the PCA coordinates unchanged; only the visualisation and metadata join are refined.
plot_data <- scores
plot_data$source_group <- factor(
  plot_data$reproductive_bacteraemia_category,
  levels = c("Reproductive", "Bacteraemia")
)

label_var <- function(component) {
  variance_row <- variance[variance$component == component, , drop = FALSE]
  pct <- as.numeric(variance_row$explained_variance_proportion)
  sprintf("%s (%0.1f%%)", component, 100 * pct)
}

make_plot <- function(x, y, stem) {
  x_label <- label_var(x)
  y_label <- label_var(y)

  gg <- ggplot(
    plot_data,
    aes(x = .data[[x]], y = .data[[y]], colour = .data$source_group, shape = .data$source_group)
  ) +
    geom_point(size = 2.4, alpha = 0.9) +
    scale_colour_manual(values = source_colours, name = "Source group") +
    scale_shape_manual(values = source_shapes, name = "Source group") +
    labs(
      x = x_label,
      y = y_label,
      title = "Nine-locus PCA by source group",
      subtitle = "Sample identifiers are retained in the accompanying PCA score table."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.45),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "#374151")
    )

  ggsave(file.path(out, paste0(stem, ".png")), gg, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), gg, width = 10, height = 8, bg = "white")
  gg
}

make_plot("PC1", "PC2", "pca_pc1_pc2")
make_plot("PC1", "PC3", "pca_pc1_pc3")

message("PCA plots written for ", nrow(scores), " genomes without point labels; sample identifiers remain unchanged in the input score table")
