#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c("ape", "ggtree", "treeio", "ggplot2", "dplyr", "readr", "stringr", "tibble")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them in the analysis environment before running this script.",
      call. = FALSE
    )
  }

  library(ape)
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("Unable to resolve project root from --file argument", call. = FALSE)
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (is.na(i)) return(default)
    if (i == length(args)) stop("Missing value for argument: ", flag, call. = FALSE)
    args[[i + 1]]
  }

  known <- c("--prefix", "--metadata", "--overwrite")
  flag_positions <- which(str_detect(args, "^--"))
  unknown <- setdiff(args[flag_positions], known)
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)

  list(
    prefix = value("--prefix", "analysis/core_genome/iqtree_primary_runs/full_core_72_gtrg4_20260819T042202Z"),
    metadata = value("--metadata", "data/metadata/curated_metadata_72_genomes.csv"),
    overwrite = "--overwrite" %in% args
  )
}

stop_if_exists <- function(paths, overwrite) {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0 && !overwrite) {
    stop(
      "Refusing to overwrite existing output(s): ",
      paste(existing, collapse = ", "),
      ". Re-run with --overwrite to replace.",
      call. = FALSE
    )
  }
}

save_plot_triplet <- function(plot, stem, output_dir, overwrite, width, height, dpi = 300) {
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  stop_if_exists(c(png_path, pdf_path, svg_path), overwrite)

  ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot = plot, width = width, height = height, bg = "white")

  svg_ok <- TRUE
  tryCatch(
    {
      ggsave(svg_path, plot = plot, width = width, height = height, bg = "white")
    },
    error = function(e) {
      svg_ok <<- FALSE
      message("SVG export skipped for ", stem, ": ", conditionMessage(e))
    }
  )

  tibble(
    stem = stem,
    png_path = png_path,
    pdf_path = pdf_path,
    svg_path = if (svg_ok) svg_path else NA_character_,
    svg_created = svg_ok
  )
}

support_label_category <- function(x) {
  if (is.na(x)) return("unlabelled")
  if (x >= 95) return(">=95")
  if (x >= 80) return("80-94")
  if (x >= 50) return("50-79")
  "<50"
}

build_unrooted_plot <- function(tree, metadata_plot, show_tip_labels, branch_length_mode = c("retain", "suppress"), title, subtitle) {
  branch_length_mode <- match.arg(branch_length_mode)
  branch_arg <- if (branch_length_mode == "retain") "branch.length" else "none"
  p <- ggtree(tree, layout = "equal_angle", branch.length = branch_arg) %<+% metadata_plot +
    geom_tippoint(aes(colour = reproductive_bacteraemia_category), size = 1.9, alpha = 0.95) +
    scale_colour_manual(
      values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"),
      drop = FALSE,
      name = "Source group"
    ) +
    labs(title = title, subtitle = subtitle, colour = "Source group") +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(8, 14, 8, 8)
    )

  if (show_tip_labels) {
    p <- p + geom_tiplab(aes(label = label), size = 1.95, linesize = 0.2, offset = 0.00015)
  }
  if (branch_length_mode == "retain") {
    p <- p + geom_treescale(width = 0.005, x = 0, y = 0, fontsize = 2.7, linesize = 0.3)
  }
  p
}

main <- function() {
  args <- parse_args()
  root <- get_project_root()

  prefix <- file.path(root, args$prefix)
  treefile_path <- paste0(prefix, ".treefile")
  contree_path <- paste0(prefix, ".contree")
  metadata_path <- file.path(root, args$metadata)

  required_inputs <- c(treefile_path, contree_path, metadata_path)
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs) > 0) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

  fig_dir <- file.path(root, "results", "figures", "core_genome_phylogeny")
  tab_dir <- file.path(root, "results", "tables", "core_genome_phylogeny")
  tree_out_dir <- file.path(root, "results", "core_genome", "trees")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tree_out_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- read_csv(metadata_path, show_col_types = FALSE)
  required_cols <- c("assembly_accession", "reproductive_bacteraemia_category")
  if (!all(required_cols %in% names(metadata))) {
    stop("Metadata missing required columns: ", paste(setdiff(required_cols, names(metadata)), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(metadata$assembly_accession)) stop("Duplicate assembly_accession values in metadata", call. = FALSE)

  metadata <- metadata %>%
    transmute(
      assembly_accession = as.character(assembly_accession),
      reproductive_bacteraemia_category = as.character(reproductive_bacteraemia_category)
    )

  if (!setequal(unique(metadata$reproductive_bacteraemia_category), c("Reproductive", "Bacteraemia"))) {
    stop("Metadata source groups must be exactly Reproductive and Bacteraemia", call. = FALSE)
  }

  combined_72 <- read.tree(treefile_path)
  support_72 <- read.tree(contree_path)

  tip_set <- unique(combined_72$tip.label)
  if (length(combined_72$tip.label) != 72 || length(tip_set) != 72) {
    stop("Combined ML tree does not contain exactly 72 unique tips", call. = FALSE)
  }

  if (!setequal(combined_72$tip.label, metadata$assembly_accession)) {
    missing_in_meta <- setdiff(combined_72$tip.label, metadata$assembly_accession)
    missing_in_tree <- setdiff(metadata$assembly_accession, combined_72$tip.label)
    stop(
      "Tree/metadata mismatch. Missing in metadata: ",
      paste(missing_in_meta, collapse = ", "),
      " | Missing in tree: ",
      paste(missing_in_tree, collapse = ", "),
      call. = FALSE
    )
  }

  support_internal <- suppressWarnings(as.numeric(support_72$node.label))
  ml_internal <- suppressWarnings(as.numeric(combined_72$node.label))
  support_source_qc <- tibble(
    tree_file = c(treefile_path, contree_path),
    internal_node_count = c(Nnode(combined_72), Nnode(support_72)),
    labelled_support_nodes = c(sum(!is.na(ml_internal)), sum(!is.na(support_internal)))
  )
  write_tsv(support_source_qc, file.path(tab_dir, "core_tree_support_source_qc.tsv"))

  reproductive_ids <- metadata %>%
    filter(reproductive_bacteraemia_category == "Reproductive") %>%
    pull(assembly_accession)
  bacteraemia_ids <- metadata %>%
    filter(reproductive_bacteraemia_category == "Bacteraemia") %>%
    pull(assembly_accession)

  if (length(reproductive_ids) != 14) stop("Expected 14 reproductive genomes, found: ", length(reproductive_ids), call. = FALSE)
  if (length(bacteraemia_ids) != 58) stop("Expected 58 bacteraemia genomes, found: ", length(bacteraemia_ids), call. = FALSE)

  reproductive_pruned_14 <- drop.tip(combined_72, bacteraemia_ids)
  bacteraemia_pruned_58 <- drop.tip(combined_72, reproductive_ids)
  reproductive_support_14 <- drop.tip(support_72, bacteraemia_ids)
  bacteraemia_support_58 <- drop.tip(support_72, reproductive_ids)

  tree_sets <- list(
    combined_72 = combined_72$tip.label,
    reproductive_pruned_14 = reproductive_pruned_14$tip.label,
    bacteraemia_pruned_58 = bacteraemia_pruned_58$tip.label
  )
  tree_counts <- vapply(tree_sets, function(x) length(unique(x)), integer(1))
  if (!identical(unname(tree_counts), c(72L, 14L, 58L))) {
    stop("Unexpected pruned tip counts: ", paste(names(tree_counts), tree_counts, collapse = ", "), call. = FALSE)
  }

  union_set <- union(tree_sets$reproductive_pruned_14, tree_sets$bacteraemia_pruned_58)
  intersection_set <- intersect(tree_sets$reproductive_pruned_14, tree_sets$bacteraemia_pruned_58)
  if (!setequal(union_set, tree_sets$combined_72)) stop("Pruned-tip union does not match combined-tip set", call. = FALSE)
  if (length(intersection_set) != 0) stop("Pruned-tip intersection must be empty", call. = FALSE)

  tree_membership <- bind_rows(
    tibble(tree_object = "combined_72", accession = sort(unique(tree_sets$combined_72))),
    tibble(tree_object = "reproductive_pruned_14", accession = sort(unique(tree_sets$reproductive_pruned_14))),
    tibble(tree_object = "bacteraemia_pruned_58", accession = sort(unique(tree_sets$bacteraemia_pruned_58)))
  )
  write_tsv(tree_membership, file.path(tab_dir, "core_tree_tip_membership.tsv"))

  tree_audit <- tree_membership %>%
    group_by(tree_object) %>%
    summarise(
      tip_count = n(),
      unique_tip_count = n_distinct(accession),
      accession_set = paste(sort(accession), collapse = ";"),
      .groups = "drop"
    )
  write_tsv(tree_audit, file.path(tab_dir, "core_tree_tip_audit.tsv"))

  combined_tree_out <- file.path(tree_out_dir, "core_72_combined_unrooted.treefile")
  reproductive_tree_out <- file.path(tree_out_dir, "core_14_reproductive_pruned_unrooted.treefile")
  bacteraemia_tree_out <- file.path(tree_out_dir, "core_58_bacteraemia_pruned_unrooted.treefile")
  stop_if_exists(c(combined_tree_out, reproductive_tree_out, bacteraemia_tree_out), args$overwrite)
  write.tree(combined_72, file = combined_tree_out)
  write.tree(reproductive_pruned_14, file = reproductive_tree_out)
  write.tree(bacteraemia_pruned_58, file = bacteraemia_tree_out)

  metadata_plot <- metadata %>% rename(label = assembly_accession)

  figure_records <- list()
  figure_source <- tibble(
    figure_stem = character(),
    tree_source_path = character(),
    tip_count = integer(),
    tree_mode = character(),
    branch_length_handling = character()
  )

  add_figure <- function(plot, stem, width, height, tree_source_path, tip_count, tree_mode, branch_length_handling) {
    rec <- save_plot_triplet(plot, stem, fig_dir, args$overwrite, width, height)
    figure_records[[length(figure_records) + 1]] <<- rec
    figure_source <<- bind_rows(
      figure_source,
      tibble(
        figure_stem = stem,
        tree_source_path = tree_source_path,
        tip_count = tip_count,
        tree_mode = tree_mode,
        branch_length_handling = branch_length_handling
      )
    )
  }

  add_figure(
    build_unrooted_plot(
      combined_72, metadata_plot, TRUE, "retain",
      "combined_72_unrooted_phylogram",
      "Unrooted ML phylogram (branch lengths retained)"
    ),
    "combined_72_unrooted_phylogram",
    15, 15, treefile_path, 72L, "phylogram", "retained"
  )
  add_figure(
    build_unrooted_plot(
      combined_72, metadata_plot, FALSE, "retain",
      "combined_72_unrooted_phylogram_presentation",
      "Unrooted ML phylogram (branch lengths retained; labels suppressed for readability)"
    ),
    "combined_72_unrooted_phylogram_presentation",
    13, 13, treefile_path, 72L, "phylogram", "retained"
  )
  add_figure(
    build_unrooted_plot(
      combined_72, metadata_plot, TRUE, "suppress",
      "combined_72_unrooted_cladogram",
      "Topology-only unrooted cladogram (branch-length scaling suppressed)"
    ),
    "combined_72_unrooted_cladogram",
    15, 15, treefile_path, 72L, "cladogram", "suppressed"
  )

  add_figure(
    build_unrooted_plot(
      reproductive_pruned_14,
      metadata_plot %>% filter(reproductive_bacteraemia_category == "Reproductive"),
      TRUE, "retain",
      "reproductive_14_pruned_unrooted_phylogram",
      "Pruned from 72-tip ML tree (branch lengths retained)"
    ),
    "reproductive_14_pruned_unrooted_phylogram",
    10, 9, treefile_path, 14L, "phylogram", "retained"
  )
  add_figure(
    build_unrooted_plot(
      reproductive_pruned_14,
      metadata_plot %>% filter(reproductive_bacteraemia_category == "Reproductive"),
      TRUE, "suppress",
      "reproductive_14_pruned_unrooted_cladogram",
      "Pruned from 72-tip ML tree, topology-only (branch-length scaling suppressed)"
    ),
    "reproductive_14_pruned_unrooted_cladogram",
    10, 9, treefile_path, 14L, "cladogram", "suppressed"
  )
  add_figure(
    build_unrooted_plot(
      bacteraemia_pruned_58,
      metadata_plot %>% filter(reproductive_bacteraemia_category == "Bacteraemia"),
      TRUE, "retain",
      "bacteraemia_58_pruned_unrooted_phylogram",
      "Pruned from 72-tip ML tree (branch lengths retained)"
    ),
    "bacteraemia_58_pruned_unrooted_phylogram",
    13, 11, treefile_path, 58L, "phylogram", "retained"
  )
  add_figure(
    build_unrooted_plot(
      bacteraemia_pruned_58,
      metadata_plot %>% filter(reproductive_bacteraemia_category == "Bacteraemia"),
      TRUE, "suppress",
      "bacteraemia_58_pruned_unrooted_cladogram",
      "Pruned from 72-tip ML tree, topology-only (branch-length scaling suppressed)"
    ),
    "bacteraemia_58_pruned_unrooted_cladogram",
    13, 11, treefile_path, 58L, "cladogram", "suppressed"
  )

  # Support-aware supplemental figure from the support-labelled consensus tree.
  support_plot <- ggtree(support_72, layout = "equal_angle", branch.length = "none") %<+% metadata_plot +
    geom_tippoint(aes(colour = reproductive_bacteraemia_category), size = 1.8) +
    geom_tiplab(aes(label = label), size = 1.9, linesize = 0.2, offset = 0.00015) +
    geom_text2(aes(
      label = ifelse(
        !isTip & !is.na(suppressWarnings(as.numeric(label))) & suppressWarnings(as.numeric(label)) >= 80,
        label,
        ""
      )
    ),
      hjust = -0.05, size = 1.9, colour = "black"
    ) +
    scale_colour_manual(values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"), drop = FALSE, name = "Source group") +
    labs(
      title = "combined_72_unrooted_cladogram_support_overlay",
      subtitle = "Support labels (UFBoot) displayed from .contree only; branch-length scaling suppressed"
    ) +
    theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold"), plot.subtitle = element_text(size = 9))
  add_figure(
    support_plot,
    "combined_72_unrooted_cladogram_support_overlay",
    15, 15, contree_path, 72L, "cladogram_with_support", "suppressed"
  )

  figure_output_manifest <- bind_rows(figure_records)
  write_tsv(figure_output_manifest, file.path(tab_dir, "core_tree_figure_output_paths.tsv"))
  write_tsv(figure_source, file.path(tab_dir, "core_tree_figure_source_table.tsv"))

  support_qc <- tibble(
    node_index = seq_along(support_72$node.label),
    node_label_raw = support_72$node.label,
    support_value = suppressWarnings(as.numeric(support_72$node.label))
  ) %>%
    mutate(support_category = vapply(support_value, support_label_category, character(1)))
  write_tsv(support_qc, file.path(tab_dir, "core_tree_support_qc.tsv"))

  support_summary <- support_qc %>%
    count(support_category, name = "node_count") %>%
    arrange(match(support_category, c(">=95", "80-94", "50-79", "<50", "unlabelled")))
  write_tsv(support_summary, file.path(tab_dir, "core_tree_support_qc_summary.tsv"))

  visualisation_qc_note <- c(
    "phylogram = quantitative branch lengths retained from the validated IQ-TREE ML tree",
    "cladogram = topology retained while branch-length scaling is suppressed for readability",
    "GCA_050472755.1 remains present with unmodified branch lengths in phylograms",
    "No independent subgroup tree inference was performed; subgroup trees were obtained by deterministic tip pruning."
  )
  write_lines(visualisation_qc_note, file.path(tab_dir, "core_tree_visualisation_qc_note.txt"))

  source_tree_sha_line <- system2("sha256sum", treefile_path, stdout = TRUE)
  source_tree_sha256 <- str_split_fixed(source_tree_sha_line, "\\s+", 2)[, 1]
  package_versions <- vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1))

  manifest_tbl <- tibble(
    field = c(
      "timestamp_utc",
      "git_commit",
      "source_iqtree_tree",
      "source_iqtree_tree_sha256",
      "metadata_file",
      "combined_tip_count",
      "reproductive_tip_count",
      "bacteraemia_tip_count",
      "pruning_method",
      "support_source_qc_table",
      "figure_source_table",
      "support_qc_table",
      "support_qc_summary_table",
      "tip_audit_table",
      "tip_membership_table",
      "tree_output_combined",
      "tree_output_reproductive",
      "tree_output_bacteraemia",
      "branch_length_mode_note",
      "package_versions",
      "figure_output_paths"
    ),
    value = c(
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      Sys.getenv("GIT_COMMIT", unset = ""),
      treefile_path,
      source_tree_sha256,
      metadata_path,
      "72",
      "14",
      "58",
      "ape::drop.tip",
      file.path(tab_dir, "core_tree_support_source_qc.tsv"),
      file.path(tab_dir, "core_tree_figure_source_table.tsv"),
      file.path(tab_dir, "core_tree_support_qc.tsv"),
      file.path(tab_dir, "core_tree_support_qc_summary.tsv"),
      file.path(tab_dir, "core_tree_tip_audit.tsv"),
      file.path(tab_dir, "core_tree_tip_membership.tsv"),
      combined_tree_out,
      reproductive_tree_out,
      bacteraemia_tree_out,
      "phylogram retained, cladogram suppressed",
      paste(paste0(names(package_versions), "=", package_versions), collapse = ";"),
      paste(figure_output_manifest$stem, collapse = ";")
    )
  )
  manifest_path <- file.path(tab_dir, "core_tree_visualisation_manifest.tsv")
  write_tsv(manifest_tbl, manifest_path)

  message("Core-genome tree visualisation complete.")
}

if (!interactive() && Sys.getenv("EFAECALIS_SKIP_MAIN", "0") != "1") {
  main()
}
