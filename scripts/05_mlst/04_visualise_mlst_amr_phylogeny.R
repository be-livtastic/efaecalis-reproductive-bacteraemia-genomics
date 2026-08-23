#!/usr/bin/env Rscript

# Visualise formal MLST and frozen AMR metadata alongside the existing nine-locus tree.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ape)
  library(phangorn)
  library(ggtree)
})

# Resolve repository-relative analytical inputs and figure destinations.
get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  derived <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  derived
}

# Parse the overwrite safeguard only.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(setdiff(args, "--overwrite"))) stop("Unknown command-line argument", call. = FALSE)
  list(overwrite = "--overwrite" %in% args)
}

# Save both publication and preview formats without silent replacement.
save_pair <- function(plot, stem, output_dir, overwrite, width, height) {
  paths <- file.path(output_dir, paste0(stem, c(".png", ".pdf")))
  if (any(file.exists(paths)) && !overwrite) stop("Refusing to overwrite figure(s): ", paste(paths[file.exists(paths)], collapse = ", "), call. = FALSE)
  ggsave(paths[[1]], plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paths[[2]], plot, width = width, height = height, bg = "white")
}

# Create the presentation-only MLST tree display while keeping midpoint rooting non-analytical.
main <- function() {
  args <- parse_args(); root <- get_project_root()
  integrated_path <- file.path(root, "results", "tables", "mlst_amr_phylogeny", "integrated_genome_mlst_amr_72.csv")
  tree_path <- file.path(root, "analysis", "phylogenomics", "72_genomes_9_locus_observed_indels", "iqtree", "efaecalis_72_genomes_9_locus_observed_indels.treefile")
  if (!all(file.exists(c(integrated_path, tree_path)))) stop("Integrated table or tree is missing", call. = FALSE)
  output_dir <- file.path(root, "results", "figures", "mlst_amr_phylogeny")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- read_csv(integrated_path, col_types = cols(.default = col_character(), AMR_burden = col_double()))
  tree <- read.tree(tree_path)
  if (length(tree$tip.label) != 72L || !setequal(tree$tip.label, metadata$Genome)) stop("Tree/integration accession mismatch", call. = FALSE)

  display_tree <- midpoint(tree)
  plot_data <- metadata |>
    transmute(
      label = Genome, Source,
      Proxy = if_else(str_detect(HLGR_proxy, " proxy detected$"), "Proxy detected", "Proxy not detected"),
      Tip_label = paste0(Genome, " | ", if_else(ST == "Unassigned", "ST unassigned", paste0("ST", ST)))
    )
  tree_plot <- ggtree(display_tree) %<+% plot_data +
    geom_tippoint(aes(colour = Source, shape = Proxy), size = 2.3) +
    geom_tiplab(aes(label = Tip_label), size = 2.2, align = TRUE, linesize = 0.2) +
    scale_colour_manual(values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"), name = "Dataset group") +
    scale_shape_manual(values = c("Proxy detected" = 17, "Proxy not detected" = 16), name = "HLGR-associated genotype") +
    theme_tree2() + theme(legend.position = "right") +
    labs(title = "Formal MLST and HLGR-associated genotype proxy on the nine-locus phylogeny",
         subtitle = "Midpoint rooting is used for presentation only; all analytical distances use the original unrooted tree")
  tree_plot <- tree_plot + xlim(NA, max(tree_plot$data$x, na.rm = TRUE) * 1.72)
  save_pair(tree_plot, "mlst_hlgr_proxy_midpoint_display_tree_72", output_dir, args$overwrite, 18, 13)

  # Record the exact R session and reproducible command sequence used for this analysis layer.
  log_dir <- file.path(root, "analysis", "logs", "mlst_amr_phylogeny")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  session_path <- file.path(log_dir, "mlst_amr_phylogeny_session_info_72.txt")
  if (file.exists(session_path) && !args$overwrite) stop("Refusing to overwrite session record: ", session_path, call. = FALSE)
  writeLines(c(
    paste("Recorded UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "Commands:",
    "Rscript scripts/05_mlst/01_prepare_pubmlst_reference.R --date=2026-08-16",
    "Rscript scripts/05_mlst/02_assign_mlst_72_genomes.R --reference-dir references/mlst/pubmlst/efaecalis_2026-08-16",
    "Rscript scripts/05_mlst/03_integrate_mlst_amr_phylogeny.R",
    "Rscript scripts/05_mlst/04_visualise_mlst_amr_phylogeny.R",
    "", capture.output(sessionInfo())
  ), session_path)
  message("Generated presentation-only MLST tree figure in ", output_dir)
}

# Run only during a deliberate Rscript invocation.
if (!interactive() && Sys.getenv("EFAECALIS_SKIP_MAIN", "0") != "1") main()
