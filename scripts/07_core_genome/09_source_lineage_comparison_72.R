#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(phangorn)
})

get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) stop("Run this script with Rscript so the project root can be resolved.", call. = FALSE)
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
}

parse_support_label <- function(label) {
  if (is.na(label) || !nzchar(trimws(label))) return(NA_real_)
  token <- trimws(label)
  if (grepl("/", token, fixed = TRUE)) token <- strsplit(token, "/", fixed = TRUE)[[1]][1]
  out <- suppressWarnings(as.numeric(token))
  if (is.na(out)) return(NA_real_)
  as.numeric(out)
}

get_descendant_tips <- function(tree, node_id) {
  child_nodes <- tree$edge[tree$edge[, 1] == node_id, 2L]
  if (length(child_nodes) == 0L) {
    return(character())
  }

  out <- character()
  for (child in child_nodes) {
    if (child <= length(tree$tip.label)) {
      out <- c(out, tree$tip.label[child])
    } else {
      out <- c(out, get_descendant_tips(tree, child))
    }
  }
  unique(out)
}

get_internal_clade_candidates <- function(tree, minimum_support = 80, minimum_size = 3L) {
  if (length(tree$node.label) == 0L) return(tibble())

  internal_nodes <- which(!is.na(tree$node.label)) + length(tree$tip.label)
  if (length(internal_nodes) == 0L) return(tibble())

  out <- purrr::map_dfr(internal_nodes, function(node_id) {
    support <- parse_support_label(tree$node.label[node_id - length(tree$tip.label)])
    tip_labels <- get_descendant_tips(tree, node_id)
    tibble(
      node_id = node_id,
      support = support,
      n_tips = length(tip_labels),
      tip_labels = list(tip_labels)
    )
  }) |>
    filter(!is.na(support), support >= minimum_support, n_tips >= minimum_size) |>
    arrange(desc(support), desc(n_tips), node_id)

  if (nrow(out) == 0L) return(out)

  keep <- rep(TRUE, nrow(out))
  for (i in seq_len(nrow(out))) {
    if (!keep[i]) next
    tip_set_i <- unlist(out$tip_labels[[i]])
    for (j in seq_len(nrow(out))) {
      if (i == j || !keep[j]) next
      tip_set_j <- unlist(out$tip_labels[[j]])
      if (length(tip_set_j) > length(tip_set_i) && all(tip_set_i %in% tip_set_j)) {
        keep[i] <- FALSE
        break
      }
    }
  }

  out[keep, ] |>
    mutate(
      clade_id = paste0("clade_", row_number()),
      .before = 1
    )
}

build_clade_summary <- function(tree, metadata, tree_name, minimum_support = 80, minimum_size = 3L) {
  candidates <- get_internal_clade_candidates(tree, minimum_support, minimum_size)
  if (nrow(candidates) == 0L) {
    return(tibble(
      tree_type = tree_name,
      clade_id = character(),
      total_isolates = integer(),
      reproductive_isolates = integer(),
      bacteraemia_isolates = integer(),
      percent_reproductive = numeric(),
      percent_bacteraemia = numeric(),
      reproductive_strain_accessions = character(),
      support = numeric(),
      n_tips = integer(),
      tip_labels = list()
    ))
  }

  source_lookup <- metadata |> select(Genome, Source) |> tibble::deframe()
  clade_rows <- map_dfr(seq_len(nrow(candidates)), function(idx) {
    tips <- unlist(candidates$tip_labels[[idx]])
    source_vec <- metadata$Source[match(tips, metadata$Genome)]
    total <- length(tips)
    reproductive_n <- sum(source_vec == "Reproductive", na.rm = TRUE)
    bacteraemia_n <- sum(source_vec == "Bacteraemia", na.rm = TRUE)
    reproductive_names <- metadata |> filter(Genome %in% tips, Source == "Reproductive") |> pull(strain)
    reproductive_accessions <- metadata |> filter(Genome %in% tips, Source == "Reproductive") |> pull(Genome)
    tibble(
      tree_type = tree_name,
      clade_id = candidates$clade_id[idx],
      total_isolates = total,
      reproductive_isolates = reproductive_n,
      bacteraemia_isolates = bacteraemia_n,
      percent_reproductive = if (total == 0L) 0 else 100 * reproductive_n / total,
      percent_bacteraemia = if (total == 0L) 0 else 100 * bacteraemia_n / total,
      reproductive_strain_accessions = paste0(reproductive_accessions, collapse = "; "),
      support = candidates$support[idx],
      n_tips = total,
      tip_labels = list(tips)
    )
  })

  clade_rows
}

build_tip_clade_map <- function(tree, clade_table) {
  if (nrow(clade_table) == 0L) {
    return(tibble(Genome = tree$tip.label, clade_id = NA_character_))
  }

  mapped <- map_dfr(seq_len(nrow(clade_table)), function(i) {
    tibble(
      Genome = unlist(clade_table$tip_labels[[i]]),
      clade_id = clade_table$clade_id[i]
    )
  })

  all_tips <- tibble(Genome = tree$tip.label)
  all_tips |> left_join(mapped, by = "Genome")
}

nearest_neighbour_summary <- function(tree, metadata, tree_name) {
  distances <- cophenetic.phylo(unroot(tree))
  reproductive <- metadata |> filter(Source == "Reproductive") |> pull(Genome)
  if (length(reproductive) == 0L) return(tibble())

  map_dfr(reproductive, function(genome) {
    candidates <- distances[genome, setdiff(colnames(distances), genome)]
    nearest <- min(candidates)
    tied <- names(candidates)[abs(candidates - nearest) <= 1e-12]
    neighbour <- sort(tied)
    neighbour_source <- metadata$Source[match(neighbour, metadata$Genome)]
    tibble(
      Genome = genome,
      Tree = tree_name,
      Nearest_neighbour = if (length(neighbour) > 0L) paste(neighbour, collapse = "; ") else NA_character_,
      Nearest_neighbour_source = if (length(neighbour_source) > 0L) paste(neighbour_source, collapse = "; ") else NA_character_,
      Patristic_distance = if (length(neighbour) > 0L) paste(format(unname(candidates[neighbour]), digits = 6, scientific = FALSE), collapse = "; ") else NA_character_,
      Num_tied_neighbours = length(neighbour)
    )
  }) |>
    left_join(metadata |> select(Genome, strain, ST), by = "Genome")
}

make_source_tree_plot <- function(tree, metadata, tree_name, plot_label) {
  display_tree <- phangorn::midpoint(tree)
  plot_data <- metadata |>
    mutate(
      tip_label = if_else(Source == "Reproductive", paste0(strain, " | ", Genome), Genome),
      Source = factor(Source, levels = c("Reproductive", "Bacteraemia"))
    )

  ggtree(display_tree) %<+% plot_data +
    geom_tippoint(aes(colour = Source), size = 2.2) +
    geom_tiplab(aes(label = tip_label), size = 2.1, offset = 0.0003, align = TRUE, linesize = 0.15) +
    scale_colour_manual(values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"), name = "Source group") +
    theme_tree2() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = paste(plot_label, "source-annotated phylogeny"),
      subtitle = "Midpoint rooting is used only for display; all interpretation uses the existing unrooted topology.",
      colour = "Source group"
    )
}

build_cross_tree_table <- function(multilocus_clades, core_clades, metadata, nn_multilocus, nn_core) {
  rep_metadata <- metadata |> filter(Source == "Reproductive") |> select(Genome, strain, ST)

  ml_map <- multilocus_clades |> filter(!is.na(clade_id)) |> select(Genome, multilocus_clade = clade_id)
  core_map <- core_clades |> filter(!is.na(clade_id)) |> select(Genome, core_clade = clade_id)

  rep_metadata |>
    left_join(ml_map, by = "Genome") |>
    left_join(core_map, by = "Genome") |>
    left_join(nn_multilocus |> select(Genome, nearest_neighbour_ml = Nearest_neighbour, nearest_neighbour_ml_source = Nearest_neighbour_source, distance_ml = Patristic_distance), by = "Genome") |>
    left_join(nn_core |> select(Genome, nearest_neighbour_core = Nearest_neighbour, nearest_neighbour_core_source = Nearest_neighbour_source, distance_core = Patristic_distance), by = "Genome")
}

build_summary_table <- function(clade_summary, neighbour_table) {
  clade_summary |>
    group_by(tree_type) |>
    summarise(
      n_clades = n(),
      n_clades_with_reproductive = sum(reproductive_isolates > 0),
      n_clades_with_mixed_sources = sum(reproductive_isolates > 0 & bacteraemia_isolates > 0),
      n_clades_reproductive_only = sum(reproductive_isolates > 0 & bacteraemia_isolates == 0),
      reproductive_isolates = sum(reproductive_isolates),
      .groups = "drop"
    ) |>
    left_join(
      neighbour_table |>
        group_by(Tree) |>
        summarise(
          nearest_bacteraemia = sum(grepl("Bacteraemia", Nearest_neighbour_source, fixed = TRUE)),
          nearest_reproductive = sum(grepl("Reproductive", Nearest_neighbour_source, fixed = TRUE)),
          .groups = "drop",
          .by = NULL
        ) |>
        rename(tree_type = Tree),
      by = c("tree_type" = "tree_type")
    )
}

main <- function() {
  root <- get_project_root()
  metadata_path <- file.path(root, "results", "tables", "mlst_amr_phylogeny", "integrated_genome_mlst_amr_72.csv")
  ml_tree_path <- file.path(root, "analysis", "phylogenomics", "72_genomes_9_locus_observed_indels", "iqtree", "efaecalis_72_genomes_9_locus_observed_indels.treefile")
  core_tree_path <- file.path(root, "results", "core_genome", "iqtree_primary", "full_core_72_gtrg4_20260819T042202Z.treefile")

  metadata_raw <- read_csv(metadata_path, show_col_types = FALSE, na = c("", "NA"))

  metadata <- metadata_raw |>
    transmute(
      Genome = as.character(Genome),
      Source = as.character(Source),
      strain = if ("strain" %in% names(metadata_raw)) as.character(strain) else Genome,
      ST = if ("ST" %in% names(metadata_raw)) as.character(ST) else NA_character_
    )

  ml_tree <- read.tree(ml_tree_path)
  core_tree <- read.tree(core_tree_path)

  if (length(ml_tree$tip.label) != 72L || length(core_tree$tip.label) != 72L) stop("Expected 72 tips in both trees.", call. = FALSE)
  if (!identical(sort(ml_tree$tip.label), sort(metadata$Genome))) stop("MLST tree tips do not match the metadata accession set.", call. = FALSE)
  if (!identical(sort(core_tree$tip.label), sort(metadata$Genome))) stop("Core-genome tree tips do not match the metadata accession set.", call. = FALSE)

  fig_dir <- file.path(root, "results", "figures", "source_lineage_comparison_72")
  tab_dir <- file.path(root, "results", "tables", "source_lineage_comparison_72")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  ml_clades <- get_internal_clade_candidates(ml_tree)
  core_clades <- get_internal_clade_candidates(core_tree)

  ml_clade_summary <- build_clade_summary(ml_tree, metadata, "multilocus")
  core_clade_summary <- build_clade_summary(core_tree, metadata, "core_genome")

  ml_tip_clades <- build_tip_clade_map(ml_tree, ml_clade_summary)
  core_tip_clades <- build_tip_clade_map(core_tree, core_clade_summary)

  ml_clade_summary <- ml_clade_summary |>
    mutate(
      tip_labels = map(tip_labels, ~ unlist(.x))
    )
  core_clade_summary <- core_clade_summary |>
    mutate(
      tip_labels = map(tip_labels, ~ unlist(.x))
    )

  ml_summary_tibble <- ml_clade_summary |>
    mutate(tree_type = "multilocus") |>
    select(tree_type, clade_id, total_isolates, reproductive_isolates, bacteraemia_isolates, percent_reproductive, percent_bacteraemia, reproductive_strain_accessions)
  core_summary_tibble <- core_clade_summary |>
    mutate(tree_type = "core_genome") |>
    select(tree_type, clade_id, total_isolates, reproductive_isolates, bacteraemia_isolates, percent_reproductive, percent_bacteraemia, reproductive_strain_accessions)

  major_clade_table <- bind_rows(ml_summary_tibble, core_summary_tibble)
  write_tsv(major_clade_table, file.path(tab_dir, "major_clade_source_composition_72.tsv"))

  clade_stack <- major_clade_table |>
    mutate(
      clade_label = paste(tree_type, clade_id, sep = "_"),
      reproductive_isolates = as.integer(reproductive_isolates),
      bacteraemia_isolates = as.integer(bacteraemia_isolates)
    ) |>
    pivot_longer(
      cols = c(reproductive_isolates, bacteraemia_isolates),
      names_to = "source_group",
      values_to = "count"
    ) |>
    mutate(
      source_group = case_when(
        source_group == "reproductive_isolates" ~ "Reproductive",
        source_group == "bacteraemia_isolates" ~ "Bacteraemia",
        TRUE ~ source_group
      )
    ) |>
    arrange(tree_type, clade_id)

  clade_stack_plot <- ggplot(clade_stack, aes(x = clade_id, y = count, fill = source_group)) +
    geom_col(position = "stack", width = 0.8) +
    facet_wrap(~ tree_type, scales = "free_x") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"), name = "Source group") +
    labs(
      title = "Major-clade source composition",
      x = "Major clade",
      y = "Isolates",
      fill = "Source group"
    )
  ggsave(file.path(fig_dir, "major_clade_source_composition_72.png"), clade_stack_plot, width = 12, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, "major_clade_source_composition_72.pdf"), clade_stack_plot, width = 12, height = 8)

  ml_plot <- make_source_tree_plot(ml_tree, metadata, "multilocus", "Multilocus")
  core_plot <- make_source_tree_plot(core_tree, metadata, "core_genome", "Core-genome")
  ggsave(file.path(fig_dir, "multilocus_source_annotated_tree_72.png"), ml_plot, width = 14, height = 18, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "multilocus_source_annotated_tree_72.pdf"), ml_plot, width = 14, height = 18, bg = "white")
  ggsave(file.path(fig_dir, "core_genome_source_annotated_tree_72.png"), core_plot, width = 14, height = 18, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "core_genome_source_annotated_tree_72.pdf"), core_plot, width = 14, height = 18, bg = "white")

  ml_neighbours <- nearest_neighbour_summary(ml_tree, metadata, "multilocus")
  core_neighbours <- nearest_neighbour_summary(core_tree, metadata, "core_genome")
  write_tsv(ml_neighbours, file.path(tab_dir, "multilocus_reproductive_nearest_neighbour_72.tsv"))
  write_tsv(core_neighbours, file.path(tab_dir, "core_genome_reproductive_nearest_neighbour_72.tsv"))

  cross_tree <- build_cross_tree_table(ml_tip_clades, core_tip_clades, metadata, ml_neighbours, core_neighbours)
  write_tsv(cross_tree, file.path(tab_dir, "reproductive_strain_cross_tree_clade_comparison_72.tsv"))

  summary_rows <- bind_rows(
    ml_clade_summary |> summarise(
      tree_type = "multilocus",
      n_reproductive_isolates = sum(reproductive_isolates),
      n_mixed_source_clades = sum(reproductive_isolates > 0 & bacteraemia_isolates > 0),
      n_reproductive_only_clades = sum(reproductive_isolates > 0 & bacteraemia_isolates == 0),
      n_nearest_neighbour_bacteraemia = sum(grepl("Bacteraemia", ml_neighbours$Nearest_neighbour_source, fixed = TRUE)),
      n_nearest_neighbour_reproductive = sum(grepl("Reproductive", ml_neighbours$Nearest_neighbour_source, fixed = TRUE)),
      n_distinct_major_clades_with_reproductive = length(unique(ml_tip_clades$clade_id[ml_tip_clades$Genome %in% metadata$Genome[metadata$Source == "Reproductive"] & !is.na(ml_tip_clades$clade_id)]))
    ),
    core_clade_summary |> summarise(
      tree_type = "core_genome",
      n_reproductive_isolates = sum(reproductive_isolates),
      n_mixed_source_clades = sum(reproductive_isolates > 0 & bacteraemia_isolates > 0),
      n_reproductive_only_clades = sum(reproductive_isolates > 0 & bacteraemia_isolates == 0),
      n_nearest_neighbour_bacteraemia = sum(grepl("Bacteraemia", core_neighbours$Nearest_neighbour_source, fixed = TRUE)),
      n_nearest_neighbour_reproductive = sum(grepl("Reproductive", core_neighbours$Nearest_neighbour_source, fixed = TRUE)),
      n_distinct_major_clades_with_reproductive = length(unique(core_tip_clades$clade_id[core_tip_clades$Genome %in% metadata$Genome[metadata$Source == "Reproductive"] & !is.na(core_tip_clades$clade_id)]))
    )
  )
  write_tsv(summary_rows, file.path(tab_dir, "lineage_summary_72.tsv"))

  message("Generated lineage comparison outputs for 72 isolates in ", fig_dir)
}

if (!interactive()) main()
