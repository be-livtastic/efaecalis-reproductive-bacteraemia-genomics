#!/usr/bin/env Rscript

# Compare unrooted trees and produce descriptive nearest-neighbour/ST summaries.
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtree)
})

find_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  list(overwrite = "--overwrite" %in% args)
}

safe_write <- function(x, path, overwrite) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  temporary <- paste0(path, ".tmp")
  readr::write_csv(x, temporary, na = "NA")
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path, call. = FALSE)
}

nearest_rows <- function(tree, metadata, tree_name) {
  distances <- cophenetic.phylo(unroot(tree))
  reproductive <- metadata |> filter(Source == "Reproductive") |> pull(Genome) |> sort()
  bind_rows(lapply(reproductive, function(genome) {
    candidates <- distances[genome, setdiff(colnames(distances), genome)]
    minimum <- min(candidates)
    tied <- names(candidates)[abs(candidates - minimum) <= 1e-12]
    tibble(
      Genome = genome,
      Tree = tree_name,
      Neighbour = sort(tied),
      Patristic_distance = unname(candidates[sort(tied)])
    )
  })) |>
    left_join(metadata |> select(Genome, ST, HLGR_proxy), by = "Genome") |>
    left_join(metadata |> select(Neighbour = Genome, Neighbour_source = Source,
                                 Neighbour_ST = ST, Neighbour_HLGR_proxy = HLGR_proxy), by = "Neighbour") |>
    mutate(Same_ST = !is.na(ST) & !is.na(Neighbour_ST) & ST == Neighbour_ST,
           HLGR_proxy_concordance = HLGR_proxy == Neighbour_HLGR_proxy)
}

main <- function() {
  args <- parse_args()
  root <- find_root()
  core_path <- file.path(root, "results/core_genome/iqtree_primary/full_core_72_gtrg4_20260819T042202Z.treefile")
  nine_path <- file.path(root, "analysis/phylogenomics/72_genomes_9_locus_observed_indels/iqtree/efaecalis_72_genomes_9_locus_observed_indels.treefile")
  metadata_path <- file.path(root, "results/tables/mlst_amr_phylogeny/integrated_genome_mlst_amr_72.csv")
  core <- read.tree(core_path)
  nine <- read.tree(nine_path)
  metadata <- read_csv(metadata_path, show_col_types = FALSE, na = c("", "NA"))
  expected <- sort(metadata$Genome)
  if (length(expected) != 72L || anyDuplicated(expected)) stop("Integrated metadata must contain 72 unique genomes", call. = FALSE)
  if (!identical(sort(core$tip.label), expected) || !identical(sort(nine$tip.label), expected)) {
    stop("Core, nine-locus and metadata tip sets are not identical", call. = FALSE)
  }
  core_unrooted <- unroot(core)
  nine_unrooted <- unroot(nine)
  comparison <- tibble(
    Metric = c("taxa", "unrooted_RF_distance", "normalised_unrooted_RF_distance"),
    Value = c(72, RF.dist(core_unrooted, nine_unrooted, normalize = FALSE),
              RF.dist(core_unrooted, nine_unrooted, normalize = TRUE))
  )
  neighbours <- bind_rows(nearest_rows(core_unrooted, metadata, "Core genome"),
                          nearest_rows(nine_unrooted, metadata, "Nine locus"))
  core_dist <- cophenetic.phylo(core_unrooted)
  assigned <- metadata |> filter(!is.na(ST))
  st_summary <- assigned |>
    count(ST, name = "Genome_count") |>
    filter(Genome_count >= 2) |>
    rowwise() |>
    mutate(
      Within_ST_pair_count = choose(Genome_count, 2),
      Within_ST_patristic_min = {
        tips <- assigned$Genome[assigned$ST == ST]
        min(core_dist[tips, tips][lower.tri(core_dist[tips, tips])])
      },
      Within_ST_patristic_median = {
        tips <- assigned$Genome[assigned$ST == ST]
        median(core_dist[tips, tips][lower.tri(core_dist[tips, tips])])
      },
      Within_ST_patristic_max = {
        tips <- assigned$Genome[assigned$ST == ST]
        max(core_dist[tips, tips][lower.tri(core_dist[tips, tips])])
      },
      Genomes_with_same_ST_nearest_neighbour = {
        tips <- assigned$Genome[assigned$ST == ST]
        sum(vapply(tips, function(tip) {
          candidate <- core_dist[tip, setdiff(colnames(core_dist), tip)]
          tied <- names(candidate)[abs(candidate - min(candidate)) <= 1e-12]
          any(metadata$ST[match(tied, metadata$Genome)] == ST, na.rm = TRUE)
        }, logical(1)))
      }
    ) |>
    ungroup()
  table_dir <- file.path(root, "results/tables/core_genome")
  safe_write(comparison, file.path(table_dir, "core_nine_locus_tree_comparison_72.csv"), args$overwrite)
  safe_write(neighbours, file.path(table_dir, "reproductive_phylogeny_neighbour_comparison_72.csv"), args$overwrite)
  safe_write(st_summary, file.path(table_dir, "st_core_resolution_summary_72.csv"), args$overwrite)

  # Midpoint rooting is used only to make the display legible; no analysis uses this object.
  display_tree <- midpoint(core_unrooted)
  display_data <- metadata |>
    mutate(ST_label = if_else(is.na(ST), "ST unassigned", paste0("ST", ST)),
           HLGR_short = if_else(grepl("not detected", HLGR_proxy, fixed = TRUE), "Proxy not detected", "Proxy detected"))
  plot <- ggtree(display_tree) %<+% display_data +
    geom_tippoint(aes(color = Source, shape = HLGR_short), size = 2) +
    geom_tiplab(aes(label = paste0(Genome, " | ", ST_label)), size = 2.1, offset = 0.002) +
    scale_color_manual(values = c(Reproductive = "#B2182B", Bacteraemia = "#2166AC")) +
    labs(title = "Core-genome phylogeny of 72 E. faecalis genomes",
         subtitle = "Midpoint rooted for display only; all analyses use the unrooted topology",
         color = "Source", shape = "HLGR-associated genotype proxy") +
    theme_tree2() + theme(legend.position = "bottom")
  figure_dir <- file.path(root, "results/figures/core_genome_phylogeny")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(figure_dir, "core_genome_source_st_hlgr_display_72.png")
  pdf_path <- file.path(figure_dir, "core_genome_source_st_hlgr_display_72.pdf")
  if ((!args$overwrite) && (file.exists(png_path) || file.exists(pdf_path))) stop("Refusing to overwrite core tree figure", call. = FALSE)
  ggsave(png_path, plot, width = 14, height = 18, dpi = 300)
  ggsave(pdf_path, plot, width = 14, height = 18)
}

main()
