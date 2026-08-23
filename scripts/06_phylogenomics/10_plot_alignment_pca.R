#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    library(ggrepel)
  }
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
plot_data$label <- ifelse(
  plot_data$reproductive_bacteraemia_category == "Reproductive",
  plot_data$strain,
  NA_character_
)

label_var <- function(component) {
  variance_row <- variance[variance$component == component, , drop = FALSE]
  pct <- as.numeric(variance_row$explained_variance_proportion)
  sprintf("%s (%0.1f%%)", component, 100 * pct)
}

make_plot <- function(x, y, stem, label_reproductive_only = TRUE) {
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
      subtitle = "Reproductive isolates are labelled; bacteraemia-associated isolates remain unlabelled."
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

  if (label_reproductive_only && any(!is.na(plot_data$label))) {
    label_df <- subset(plot_data, !is.na(label))
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      gg <- gg + geom_text_repel(
        data = label_df,
        aes(label = .data$label),
        color = "#111827",
        size = 3.1,
        box.padding = 0.35,
        point.padding = 0.2,
        segment.alpha = 0.7,
        max.overlaps = Inf,
        force = 1.2
      )
    } else {
      gg <- gg + geom_text(
        data = label_df,
        aes(label = .data$label),
        color = "#111827",
        size = 3.1,
        vjust = -0.6
      )
    }
  }

  ggsave(file.path(out, paste0(stem, ".png")), gg, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), gg, width = 10, height = 8, bg = "white")
  gg
}

make_plot("PC1", "PC2", "pca_pc1_pc2")
make_plot("PC1", "PC3", "pca_pc1_pc3")

# Export a concise PCA summary table with ST information where available.
project_root <- normalizePath(file.path(dirname(scores_path), "..", "..", "..", ".."), mustWork = TRUE)
integrated_path <- file.path(
  project_root,
  "results",
  "tables",
  "mlst_amr_phylogeny",
  "integrated_genome_mlst_amr_72.csv"
)
if (file.exists(integrated_path)) {
  mlst <- read.csv(integrated_path, check.names = FALSE, stringsAsFactors = FALSE)
  mlst_table <- mlst[, c("Genome", "ST")]
  names(mlst_table)[1] <- "assembly_accession"
  pca_table <- merge(
    plot_data[, c("assembly_accession", "strain", "reproductive_bacteraemia_category", "PC1", "PC2", "PC3")],
    mlst_table,
    by = "assembly_accession",
    all.x = TRUE,
    sort = FALSE
  )
  pca_table <- pca_table[, c("assembly_accession", "strain", "reproductive_bacteraemia_category", "PC1", "PC2", "PC3", "ST")]
  write.table(
    pca_table,
    file.path(dirname(scores_path), "pca_scores_with_metadata_72.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

message("PCA plots written for ", nrow(scores), " genomes with reproductive-only labels and source-group colouring")
