#!/usr/bin/env Rscript

# Build accepted-detection matrices while making every unresolved call prominent.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

find_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
}

parse_args <- function() list(overwrite = "--overwrite" %in% commandArgs(trailingOnly = TRUE))

safe_write <- function(x, path, overwrite) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  temporary <- paste0(path, ".tmp")
  write_csv(x, temporary, na = "NA")
  if (!file.rename(temporary, path)) stop("Could not atomically write ", path, call. = FALSE)
}

composite_ebp_status <- function(a, b, c) {
  values <- c(a, b, c)
  if (all(values == "accepted_present")) return("accepted_present")
  if (any(values == "ambiguous_multiple_hit")) return("ambiguous_multiple_hit")
  if (any(values == "review_required")) return("review_required")
  if (any(values == "flagged_partial")) return("flagged_partial")
  "not_detected"
}

main <- function() {
  args <- parse_args()
  root <- find_root()
  status_path <- file.path(root, "data/processed/virulence_adherence/virulence_resolved_target_status_by_genome_72.csv")
  integrated_path <- file.path(root, "results/tables/mlst_amr_phylogeny/integrated_genome_mlst_amr_72.csv")
  statuses <- read_csv(status_path, show_col_types = FALSE, na = c("", "NA"))
  integrated <- read_csv(integrated_path, show_col_types = FALSE, na = c("", "NA"))
  allowed <- c("accepted_present", "flagged_partial", "review_required", "ambiguous_multiple_hit", "not_detected")
  if (n_distinct(statuses$Genome) != 72L || nrow(integrated) != 72L || anyDuplicated(integrated$Genome)) {
    stop("Status and frozen integration inputs must contain 72 unique genomes", call. = FALSE)
  }
  if (!all(statuses$Detection_status %in% allowed)) stop("Unknown detection status", call. = FALSE)
  source_counts <- integrated |> count(Source)
  observed_counts <- setNames(source_counts$n, source_counts$Source)
  required_sources <- c("Reproductive", "Bacteraemia")
  if (length(intersect(required_sources, names(observed_counts))) != 2L || !all(required_sources %in% names(observed_counts))) {
    stop("The integrated metadata must include both Reproductive and Bacteraemia source groups.", call. = FALSE)
  }
  expected_targets <- c("ace", "efaA", "ebpA", "ebpB", "ebpC", "aggregation_substance_family_detected", "asa1_specific", "gelE", "sprE", "esp", "cylA")
  if (nrow(statuses) != 72L * length(expected_targets) ||
      !setequal(unique(statuses$Target_gene), expected_targets)) {
    stop("The resolved screening input must contain exactly the eleven approved targets for all 72 genomes", call. = FALSE)
  }
  source_check <- statuses |>
    distinct(Genome, Source) |>
    left_join(integrated |> select(Genome, Expected_source = Source), by = "Genome", relationship = "many-to-one")
  if (nrow(source_check) != 72L || any(is.na(source_check$Expected_source)) ||
      any(source_check$Source != source_check$Expected_source)) {
    stop("Virulence screening sources must match the frozen 72-genome integration table exactly", call. = FALSE)
  }
  grid <- expand_grid(Genome = integrated$Genome, Target_gene = expected_targets)
  if (nrow(anti_join(grid, statuses, by = c("Genome", "Target_gene"))) ||
      any(count(statuses, Genome, Target_gene)$n != 1L)) {
    stop("Every genome-target pair must have exactly one final status", call. = FALSE)
  }

  # Add the Ebp composite with conservative propagation of unresolved component calls.
  ebp <- statuses |>
    filter(Target_gene %in% c("ebpA", "ebpB", "ebpC")) |>
    select(Genome, Target_gene, Detection_status) |>
    pivot_wider(names_from = Target_gene, values_from = Detection_status) |>
    rowwise() |>
    transmute(Genome, Source = integrated$Source[match(Genome, integrated$Genome)],
              Target_gene = "Ebp_operon_complete",
              Detection_status = composite_ebp_status(ebpA, ebpB, ebpC),
              Status_reason = "Derived from ebpA, ebpB and ebpC final statuses") |>
    ungroup()
  analysis_status <- bind_rows(statuses, ebp)
  feature_order <- c(expected_targets, "Ebp_operon_complete")
  denominators <- tibble(Source = required_sources, Total = as.integer(observed_counts[required_sources]))
  counts <- analysis_status |>
    count(Target_gene, Source, Detection_status) |>
    complete(Target_gene = feature_order, Source = denominators$Source,
             Detection_status = allowed, fill = list(n = 0L)) |>
    pivot_wider(names_from = Detection_status, values_from = n, values_fill = 0L) |>
    left_join(denominators, by = "Source") |>
    mutate(Status_sum = accepted_present + flagged_partial + review_required + ambiguous_multiple_hit + not_detected)
  if (any(counts$Status_sum != counts$Total)) stop("Final status counts do not reconcile to group denominators", call. = FALSE)
  prevalence <- counts |>
    mutate(Accepted_percent = round(100 * accepted_present / Total, 2),
           Maximum_possible_count = accepted_present + flagged_partial + review_required + ambiguous_multiple_hit,
           Maximum_possible_percent_if_all_unresolved_confirmed = round(100 * Maximum_possible_count / Total, 2),
           Unresolved_excluded_from_binary_presence = flagged_partial + review_required + ambiguous_multiple_hit,
           Comparison_role = case_when(Target_gene == "aggregation_substance_family_detected" ~ "Primary aggregation-substance comparison",
                                       Target_gene == "asa1_specific" ~ "Secondary lower-confidence narrow breakdown",
                                       TRUE ~ "Other targeted determinant")) |>
    select(Gene = Target_gene, Source, Accepted_present_count = accepted_present, Total,
           Accepted_percent, Flagged_partial_count = flagged_partial,
           Review_required_count = review_required, Ambiguous_multiple_hit_count = ambiguous_multiple_hit,
           Not_detected_count = not_detected, Maximum_possible_count,
           Maximum_possible_percent_if_all_unresolved_confirmed,
           Unresolved_excluded_from_binary_presence, Comparison_role)

  prevalence <- prevalence |>
    left_join(
      prevalence |>
        select(Gene, Source, Accepted_percent, Maximum_possible_percent_if_all_unresolved_confirmed) |>
        pivot_wider(names_from = Source, values_from = c(Accepted_percent, Maximum_possible_percent_if_all_unresolved_confirmed), names_sep = "_"),
      by = "Gene"
    )

  gene_summary <- prevalence |>
    select(Gene, Accepted_percent_Reproductive, Accepted_percent_Bacteraemia,
          Maximum_possible_percent_if_all_unresolved_confirmed_Reproductive,
          Maximum_possible_percent_if_all_unresolved_confirmed_Bacteraemia) |>
    distinct() |>
    mutate(
      Near_universal_status = case_when(
        Accepted_percent_Reproductive >= 95 & Accepted_percent_Bacteraemia >= 95 ~ "Near-universal",
        (Accepted_percent_Reproductive < 95 & Maximum_possible_percent_if_all_unresolved_confirmed_Reproductive >= 95) |
         (Accepted_percent_Bacteraemia < 95 & Maximum_possible_percent_if_all_unresolved_confirmed_Bacteraemia >= 95) ~
         "Potentially near-universal pending unresolved calls",
        TRUE ~ "Not near-universal"
      )
    ) |>
    select(Gene, Near_universal_status)

  prevalence <- prevalence |> left_join(gene_summary, by = "Gene")

  # Construct the 72-row binary matrix from accepted calls only.
  matrix <- analysis_status |>
    transmute(Genome, Target_gene, Present = as.integer(Detection_status == "accepted_present")) |>
    pivot_wider(names_from = Target_gene, values_from = Present, values_fill = 0L) |>
    right_join(integrated |> select(Genome, Source), by = "Genome") |>
    select(Genome, Source, all_of(feature_order)) |>
    arrange(factor(Source, levels = c("Reproductive", "Bacteraemia")), Genome)
  if (nrow(matrix) != 72L || anyDuplicated(matrix$Genome)) stop("Binary matrix must contain 72 unique genomes", call. = FALSE)

  # Run gene-level Fisher tests after the accepted binary matrix has been built.
  n_rep <- denominators$Total[denominators$Source == "Reproductive"]
  n_bac <- denominators$Total[denominators$Source == "Bacteraemia"]
  fisher_rows <- lapply(feature_order, function(gene) {
    rep_present <- sum(matrix[[gene]][matrix$Source == "Reproductive"])
    bac_present <- sum(matrix[[gene]][matrix$Source == "Bacteraemia"])
    total_present <- rep_present + bac_present
    if (total_present < 2L) return(NULL)
    test <- fisher.test(base::matrix(c(rep_present, n_rep - rep_present, bac_present, n_bac - bac_present), nrow = 2, byrow = TRUE), alternative = "two.sided")
    tibble(Gene = gene, Comparison_role = case_when(gene == "aggregation_substance_family_detected" ~ "Primary aggregation-substance comparison",
                                                  gene == "asa1_specific" ~ "Secondary lower-confidence narrow breakdown",
                                                  TRUE ~ "Other targeted determinant"),
          Reproductive_count = rep_present, Reproductive_total = n_rep,
          Reproductive_percent = round(100 * rep_present / n_rep, 2),
          Bacteraemia_count = bac_present, Bacteraemia_total = n_bac,
          Bacteraemia_percent = round(100 * bac_present / n_bac, 2),
          Prevalence_difference_percentage_points = round(100 * rep_present / n_rep - 100 * bac_present / n_bac, 2),
          Odds_ratio = unname(test$estimate), Raw_p_value = test$p.value)
  })
  fisher <- bind_rows(fisher_rows) |>
    mutate(BH_adjusted_q_value = p.adjust(Raw_p_value, method = "BH"),
          BH_adjusted_p_value = BH_adjusted_q_value,
          Significant_q_lt_0.05 = BH_adjusted_q_value < 0.05) |>
    arrange(Raw_p_value)
  insufficient <- tibble(Gene = feature_order, Total_accepted_occurrences = vapply(feature_order, function(x) sum(matrix[[x]]), integer(1))) |>
    filter(Total_accepted_occurrences < 2L) |>
    mutate(Note = "insufficient occurrences for testing")
  unresolved <- prevalence |>
    filter(Review_required_count > 0 | Ambiguous_multiple_hit_count > 0 | Flagged_partial_count > 0) |>
    arrange(Gene, Source)
  integrated_output <- integrated |>
    left_join(matrix |> select(-Source), by = "Genome", relationship = "one-to-one")
  st_summary <- integrated_output |>
    filter(!is.na(ST)) |>
    group_by(ST) |>
    filter(n() >= 2L) |>
    summarise(Genome_count = n(), Reproductive_count = sum(Source == "Reproductive"),
             Bacteraemia_count = sum(Source == "Bacteraemia"),
             across(all_of(feature_order), list(count = sum, percent = ~round(100 * mean(.x), 2))), .groups = "drop")

  table_dir <- file.path(root, "results/tables/virulence_adherence")
  safe_write(matrix, file.path(table_dir, "virulence_adherence_presence_absence_72.csv"), args$overwrite)
  safe_write(prevalence, file.path(table_dir, "virulence_prevalence_with_unresolved_72.csv"), args$overwrite)
  safe_write(unresolved, file.path(table_dir, "virulence_unresolved_calls_72.csv"), args$overwrite)
  safe_write(fisher, file.path(table_dir, "virulence_fisher_tests_72.csv"), args$overwrite)
  safe_write(insufficient, file.path(table_dir, "virulence_insufficient_occurrences_72.csv"), args$overwrite)
  safe_write(integrated_output, file.path(table_dir, "integrated_genome_mlst_amr_virulence_72.csv"), args$overwrite)
  safe_write(st_summary, file.path(table_dir, "virulence_prevalence_by_st_72.csv"), args$overwrite)

  # Show accepted prevalence as the main estimate and the unresolved upper bound as a capped uncertainty interval.
  near_universal_summary <- prevalence |>
    filter(Near_universal_status == "Near-universal") |>
    select(
      Gene,
      Accepted_percent_Reproductive,
      Accepted_percent_Bacteraemia,
      Maximum_possible_percent_if_all_unresolved_confirmed_Reproductive,
      Maximum_possible_percent_if_all_unresolved_confirmed_Bacteraemia
    ) |>
    distinct() |>
    arrange(Gene)
  safe_write(near_universal_summary, file.path(table_dir, "virulence_near_universal_genes_72.csv"), args$overwrite)

  variable_genes <- prevalence |>
    filter(Near_universal_status != "Near-universal") |>
    pull(Gene) |>
    unique()

  variable_unresolved <- prevalence |>
    filter(Near_universal_status != "Near-universal") |>
    select(Gene, Source, Accepted_present_count, Flagged_partial_count, Review_required_count, Ambiguous_multiple_hit_count, Maximum_possible_count) |>
    arrange(Gene, factor(Source, levels = c("Reproductive", "Bacteraemia")))
  safe_write(variable_unresolved, file.path(table_dir, "virulence_variable_genes_unresolved_calls_72.csv"), args$overwrite)

  plot_data <- prevalence |>
    filter(Near_universal_status != "Near-universal") |>
    left_join(
      fisher |>
        filter(Significant_q_lt_0.05 == TRUE) |>
        select(Gene, BH_adjusted_q_value, Significant_q_lt_0.05),
      by = "Gene"
    ) |>
    mutate(
      Gene = factor(Gene, levels = unique(Gene[order(-abs(Accepted_percent_Reproductive - Accepted_percent_Bacteraemia))])),
      Source = factor(Source, levels = c("Reproductive", "Bacteraemia")),
      y_base = as.numeric(Gene),
      y_num = y_base + if_else(Source == "Reproductive", 0.12, -0.12),
      Has_interval = Maximum_possible_percent_if_all_unresolved_confirmed > Accepted_percent,
      significance_marker = if_else(Significant_q_lt_0.05 == TRUE & Source == "Reproductive", "*", NA_character_)
    )

  plot <- ggplot(plot_data, aes(x = Accepted_percent, y = y_num, colour = Source, shape = Source)) +
    geom_segment(aes(x = Accepted_percent, xend = Maximum_possible_percent_if_all_unresolved_confirmed,
                    y = y_num, yend = y_num), linewidth = 1.1, alpha = 0.75) +
    geom_point(size = 3.5) +
    geom_text(data = plot_data |> filter(!is.na(significance_marker)),
              aes(label = significance_marker, x = pmax(Accepted_percent, 0) + 2),
              vjust = 0.5, size = 4, colour = "black") +
    scale_colour_manual(values = c(Reproductive = "#D55E00", Bacteraemia = "#0072B2"), name = "Source") +
    scale_shape_manual(values = c(Reproductive = 15, Bacteraemia = 17), name = "Source") +
    scale_y_continuous(
      breaks = seq_along(levels(plot_data$Gene)),
      labels = levels(plot_data$Gene),
      limits = c(0.5, length(levels(plot_data$Gene)) + 0.5)
    ) +
    scale_x_continuous(limits = c(0, 120), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = "Accepted virulence prevalence among variable genes",
      subtitle = paste0("Reproductive (n = ", denominators$Total[denominators$Source == "Reproductive"], "); Bacteraemia (n = ", denominators$Total[denominators$Source == "Bacteraemia"], ")"),
      caption = "Point = accepted prevalence; horizontal capped interval = accepted prevalence to maximum possible prevalence if all unresolved calls are ultimately accepted.",
      x = "Prevalence (%)",
      y = NULL,
      colour = "Source",
      shape = "Source"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0), axis.text.y = element_text(size = 9))
  figure_dir <- file.path(root, "results/figures/virulence_adherence")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(figure_dir, "virulence_prevalence_variable_genes_72.png")
  pdf <- file.path(figure_dir, "virulence_prevalence_variable_genes_72.pdf")
  if ((!args$overwrite) && (file.exists(png) || file.exists(pdf))) stop("Refusing to overwrite variable-gene prevalence figure", call. = FALSE)
  ggsave(png, plot, width = 12, height = 7, dpi = 300)
  ggsave(pdf, plot, width = 12, height = 7)

  # Write a Results-ready summary before statistical interpretation.
  summary_path <- file.path(table_dir, "virulence_results_summary_72.md")
  if (file.exists(summary_path) && !args$overwrite) stop("Refusing to overwrite Results summary", call. = FALSE)
  lines <- c("# Targeted virulence/adherence results summary", "",
             "Prevalence below means accepted genomic detection only. Review-required and ambiguous calls were excluded from binary presence.",
             "The primary aggregation-substance comparison is `aggregation_substance_family_detected`; `asa1_specific` is a secondary, lower-confidence narrow breakdown.", "")
  for (gene in feature_order) {
    gene_rows <- prevalence |> filter(Gene == gene)
    lines <- c(lines, paste0("- ", gene, ": ", paste0(gene_rows$Source, " accepted ", gene_rows$Accepted_present_count,
      "/", gene_rows$Total, " (", gene_rows$Accepted_percent, "%); review-required ", gene_rows$Review_required_count,
      "; ambiguous ", gene_rows$Ambiguous_multiple_hit_count, "; partial ", gene_rows$Flagged_partial_count, collapse = "; "), "."))
  }
  if (nrow(unresolved)) lines <- c(lines, "", "**Warning:** unresolved or partial calls affect the targets listed above and must be considered before interpreting prevalence or significance.")
  writeLines(lines, summary_path)
  if (nrow(unresolved)) {
    warning("Review-required, ambiguous or partial calls are present. Inspect virulence_unresolved_calls_72.csv before statistical results.", call. = FALSE)
  }
}

# Permit the companion unit test to load pure helpers without running the analysis.
if (!identical(Sys.getenv("EFAECALIS_SKIP_MAIN", unset = "0"), "1")) main()
