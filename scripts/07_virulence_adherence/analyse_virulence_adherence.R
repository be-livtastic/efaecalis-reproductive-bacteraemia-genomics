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
  if (!identical(unname(observed_counts[c("Reproductive", "Bacteraemia")]), c(14L, 58L))) {
    stop("Frozen source groups must contain exactly 14 reproductive and 58 bacteraemia genomes", call. = FALSE)
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
  denominators <- tibble(Source = c("Reproductive", "Bacteraemia"), Total = c(14L, 58L))
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
           Maximum_possible_count = accepted_present + review_required + ambiguous_multiple_hit,
           Maximum_possible_percent_if_all_unresolved_confirmed = round(100 * Maximum_possible_count / Total, 2),
           Unresolved_excluded_from_binary_presence = review_required + ambiguous_multiple_hit,
           Comparison_role = case_when(Target_gene == "aggregation_substance_family_detected" ~ "Primary aggregation-substance comparison",
                                       Target_gene == "asa1_specific" ~ "Secondary lower-confidence narrow breakdown",
                                       TRUE ~ "Other targeted determinant")) |>
    select(Gene = Target_gene, Source, Accepted_present_count = accepted_present, Total,
           Accepted_percent, Flagged_partial_count = flagged_partial,
           Review_required_count = review_required, Ambiguous_multiple_hit_count = ambiguous_multiple_hit,
           Not_detected_count = not_detected, Maximum_possible_count,
           Maximum_possible_percent_if_all_unresolved_confirmed,
           Unresolved_excluded_from_binary_presence, Comparison_role)

  # Construct the 72-row binary matrix from accepted calls only.
  matrix <- analysis_status |>
    transmute(Genome, Target_gene, Present = as.integer(Detection_status == "accepted_present")) |>
    pivot_wider(names_from = Target_gene, values_from = Present, values_fill = 0L) |>
    right_join(integrated |> select(Genome, Source), by = "Genome") |>
    select(Genome, Source, all_of(feature_order)) |>
    arrange(factor(Source, levels = c("Reproductive", "Bacteraemia")), Genome)
  if (nrow(matrix) != 72L || anyDuplicated(matrix$Genome)) stop("Binary matrix must contain 72 unique genomes", call. = FALSE)

  # Run gene-level Fisher tests after the accepted binary matrix has been built.
  fisher_rows <- lapply(feature_order, function(gene) {
    rep_present <- sum(matrix[[gene]][matrix$Source == "Reproductive"])
    bac_present <- sum(matrix[[gene]][matrix$Source == "Bacteraemia"])
    total_present <- rep_present + bac_present
    if (total_present < 2L) return(NULL)
    test <- fisher.test(base::matrix(c(rep_present, 14L - rep_present, bac_present, 58L - bac_present), nrow = 2, byrow = TRUE), alternative = "two.sided")
    tibble(Gene = gene, Comparison_role = case_when(gene == "aggregation_substance_family_detected" ~ "Primary aggregation-substance comparison",
                                                   gene == "asa1_specific" ~ "Secondary lower-confidence narrow breakdown",
                                                   TRUE ~ "Other targeted determinant"),
           Reproductive_count = rep_present, Reproductive_total = 14L,
           Reproductive_percent = round(100 * rep_present / 14, 2),
           Bacteraemia_count = bac_present, Bacteraemia_total = 58L,
           Bacteraemia_percent = round(100 * bac_present / 58, 2),
           Prevalence_difference_percentage_points = round(100 * rep_present / 14 - 100 * bac_present / 58, 2),
           Odds_ratio = unname(test$estimate), Raw_p_value = test$p.value)
  })
  fisher <- bind_rows(fisher_rows) |>
    mutate(BH_adjusted_p_value = p.adjust(Raw_p_value, method = "BH"), Significant_q_lt_0.05 = BH_adjusted_p_value < 0.05) |>
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

  # Show accepted prevalence and the unresolved upper bound with explicit counts on every affected bar.
  plot_data <- prevalence |>
    mutate(Gene = factor(Gene, levels = rev(feature_order)),
           Unresolved_label = if_else(Review_required_count + Ambiguous_multiple_hit_count + Flagged_partial_count > 0,
             paste0("review=", Review_required_count, "; ambiguous=", Ambiguous_multiple_hit_count,
                    "; partial=", Flagged_partial_count), ""))
  plot <- ggplot(plot_data, aes(Accepted_percent, Gene, fill = Source)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_segment(aes(x = Accepted_percent, xend = Maximum_possible_percent_if_all_unresolved_confirmed,
                     yend = Gene, color = Source), linewidth = 1.1, position = position_dodge(width = 0.75)) +
    geom_text(aes(label = Unresolved_label), position = position_dodge(width = 0.75),
              hjust = -0.05, size = 2.6, color = "black") +
    scale_fill_manual(values = c(Reproductive = "#B2182B", Bacteraemia = "#2166AC")) +
    scale_color_manual(values = c(Reproductive = "#B2182B", Bacteraemia = "#2166AC"), guide = "none") +
    scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
    labs(title = "Accepted-detection prevalence of targeted virulence/adherence determinants",
         subtitle = "Thin extensions show the maximum possible percentage if all unresolved calls were confirmed",
         caption = "Review-required and ambiguous calls are excluded from binary presence; accepted prevalence may therefore underestimate genomic prevalence pending review.",
         x = "Accepted-detection prevalence (%)", y = NULL, fill = "Source") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
  figure_dir <- file.path(root, "results/figures/virulence_adherence")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(figure_dir, "virulence_prevalence_with_unresolved_72.png")
  pdf <- file.path(figure_dir, "virulence_prevalence_with_unresolved_72.pdf")
  if ((!args$overwrite) && (file.exists(png) || file.exists(pdf))) stop("Refusing to overwrite prevalence figure", call. = FALSE)
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
