#!/usr/bin/env Rscript

# Integrate formal MLST, frozen AMR outputs, metadata, and the unrooted nine-locus tree.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(ape)
})

PRIMARY_PROXY <- "aac(6')-Ie/aph(2'')-Ia"
EXPECTED_GROUPS <- c(Reproductive = 14L, Bacteraemia = 58L)
PERMUTATIONS <- 100000L
MANTEL_PERMUTATIONS <- 9999L
RANDOM_SEED <- 20260816L

# Resolve project-relative frozen inputs and new output paths.
get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  derived <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  derived
}

# Parse the overwrite safeguard without accepting hidden analysis options.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  unknown <- setdiff(args, "--overwrite")
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  list(overwrite = "--overwrite" %in% args)
}

# Refuse partial reruns that could mix upstream states.
preflight_outputs <- function(paths, overwrite) {
  if (anyDuplicated(normalizePath(paths, mustWork = FALSE))) stop("Duplicate output path in registry", call. = FALSE)
  existing <- paths[file.exists(paths)]
  if (length(existing) && !overwrite) stop("Refusing to overwrite: ", paste(existing, collapse = ", "), call. = FALSE)
  walk(unique(dirname(paths)), dir.create, recursive = TRUE, showWarnings = FALSE)
}

# Validate one-to-one accession joins against the canonical 72-genome manifest.
validate_genomes <- function(data, genome_col, label) {
  genomes <- data[[genome_col]]
  if (nrow(data) != 72L || n_distinct(genomes) != 72L || anyNA(genomes)) stop(label, " must contain 72 unique genomes", call. = FALSE)
}

# Calculate the prespecified within-group positive-pair statistic.
pair_concentration_statistic <- function(groups, positive) {
  counts <- table(groups[positive])
  sum(counts * (counts - 1) / 2)
}

# Run the one-sided exploratory concentration permutation with fixed labels and group sizes.
run_group_permutation <- function(groups, positive, analysis_name) {
  keep <- !is.na(groups) & groups != "" & groups != "Unassigned"
  groups <- as.character(groups[keep]); positive <- as.logical(positive[keep])
  n_positive <- sum(positive); n_genomes <- length(positive)
  if (n_positive < 2L) stop("Too few assigned proxy-positive genomes for ", analysis_name, call. = FALSE)
  observed <- pair_concentration_statistic(groups, positive)
  set.seed(RANDOM_SEED)
  null <- replicate(PERMUTATIONS, {
    permuted <- rep(FALSE, n_genomes)
    permuted[sample.int(n_genomes, n_positive)] <- TRUE
    pair_concentration_statistic(groups, permuted)
  })
  null_sd <- sd(null)
  tibble(
    Analysis = analysis_name,
    Interpretation = "Exploratory population-structure test; not a clinical association test",
    Statistic = "sum_s choose(k_s, 2)", Observed_T = observed,
    Null_mean = mean(null), Null_SD = null_sd,
    Standardized_effect = if_else(null_sd > 0, (observed - mean(null)) / null_sd, NA_real_),
    Empirical_one_sided_p = (1 + sum(null >= observed)) / (PERMUTATIONS + 1),
    Assigned_genomes = n_genomes, Assigned_proxy_positive = n_positive,
    Permutations = PERMUTATIONS, Seed = RANDOM_SEED, Extreme_tail = "T_perm >= T_obs"
  )
}

# Calculate mean carrier-to-carrier distance on the original unrooted tree.
mean_selected_distance <- function(distance_matrix, selected) {
  n <- length(selected)
  sum(distance_matrix[selected, selected, drop = FALSE]) / (n * (n - 1))
}

# Run unrestricted and source-stratified one-sided phylogenetic clustering permutations.
run_phylogenetic_permutations <- function(distance_matrix, metadata) {
  positive <- metadata$HLGR_proxy == "HLGR-associated genotype proxy detected"
  observed_indices <- which(positive)
  observed <- mean_selected_distance(distance_matrix, observed_indices)
  set.seed(RANDOM_SEED)
  unrestricted <- replicate(PERMUTATIONS, mean_selected_distance(distance_matrix, sample.int(nrow(metadata), sum(positive))))
  source_indices <- split(seq_len(nrow(metadata)), metadata$Source)
  source_positive <- table(metadata$Source[positive])
  set.seed(RANDOM_SEED)
  stratified <- replicate(PERMUTATIONS, {
    selected <- unlist(imap(source_indices, ~ sample(.x, source_positive[[.y]])), use.names = FALSE)
    mean_selected_distance(distance_matrix, selected)
  })
  bind_rows(
    tibble(Null_model = "Unrestricted", Null_values = list(unrestricted)),
    tibble(Null_model = "Source-stratified (preserves 5/14 and 24/58)", Null_values = list(stratified))
  ) |>
    mutate(
      Hypothesis = "Proxy-positive genomes have smaller mean pairwise patristic distance than random",
      Observed_mean_pairwise_distance = observed,
      Null_mean = map_dbl(Null_values, mean), Null_SD = map_dbl(Null_values, sd),
      Standardized_effect = (observed - Null_mean) / Null_SD,
      Empirical_one_sided_p = map_dbl(Null_values, ~ (1 + sum(.x <= observed)) / (PERMUTATIONS + 1)),
      Proxy_positive_genomes = sum(positive), Permutations = PERMUTATIONS,
      Seed = RANDOM_SEED, Extreme_tail = "D_perm <= D_obs"
    ) |>
    select(-Null_values)
}

# Calculate genome-pair Jaccard distances with all-zero pairs treated as identical.
jaccard_matrix <- function(binary_matrix) {
  n <- nrow(binary_matrix)
  output <- matrix(0, n, n, dimnames = list(rownames(binary_matrix), rownames(binary_matrix)))
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    union <- sum(binary_matrix[i, ] | binary_matrix[j, ])
    distance <- if (union == 0L) 0 else 1 - sum(binary_matrix[i, ] & binary_matrix[j, ]) / union
    output[i, j] <- output[j, i] <- distance
  }
  output
}

# Rank a symmetric distance matrix once for an efficient Spearman Mantel permutation.
rank_distance_matrix <- function(distance_matrix) {
  upper <- upper.tri(distance_matrix)
  ranked <- matrix(0, nrow(distance_matrix), ncol(distance_matrix), dimnames = dimnames(distance_matrix))
  ranked[upper] <- rank(distance_matrix[upper], ties.method = "average")
  ranked <- ranked + t(ranked)
  ranked
}

# Run a two-sided exploratory Mantel-style label permutation.
run_mantel <- function(amr_distance, phylogenetic_distance) {
  upper <- upper.tri(amr_distance)
  amr_rank <- rank_distance_matrix(amr_distance)
  phylo_vector <- rank_distance_matrix(phylogenetic_distance)[upper]
  observed <- cor(amr_rank[upper], phylo_vector)
  set.seed(RANDOM_SEED)
  null <- replicate(MANTEL_PERMUTATIONS, {
    permutation <- sample.int(nrow(amr_rank))
    cor(amr_rank[permutation, permutation][upper], phylo_vector)
  })
  tibble(
    Analysis = "Exploratory AMR Jaccard versus unrooted nine-locus patristic distance",
    Statistic = "Spearman rank correlation implemented as Pearson correlation of ranked distances",
    Observed_correlation = observed, Null_mean = mean(null), Null_SD = sd(null),
    Empirical_two_sided_p = (1 + sum(abs(null) >= abs(observed))) / (MANTEL_PERMUTATIONS + 1),
    Permutations = MANTEL_PERMUTATIONS, Seed = RANDOM_SEED
  )
}

# Preserve all pairwise source, ST, and study overlap for transparent descriptive summaries.
build_pairwise_audit <- function(metadata, amr_distance, phylogenetic_distance) {
  pairs <- t(combn(seq_len(nrow(metadata)), 2))
  map_dfr(seq_len(nrow(pairs)), function(index) {
    i <- pairs[index, 1]; j <- pairs[index, 2]
    assigned_pair <- metadata$ST[[i]] != "Unassigned" && metadata$ST[[j]] != "Unassigned"
    tibble(
      Genome_1 = metadata$Genome[[i]], Genome_2 = metadata$Genome[[j]],
      AMR_Jaccard_distance = amr_distance[i, j], AMR_Jaccard_similarity = 1 - amr_distance[i, j],
      Nine_locus_patristic_distance = phylogenetic_distance[i, j],
      ST_pair_evaluable = assigned_pair,
      ST_relationship = if (!assigned_pair) "ST comparison unavailable" else if (metadata$ST[[i]] == metadata$ST[[j]]) "Same ST" else "Different ST",
      Source_relationship = if (metadata$Source[[i]] == metadata$Source[[j]]) "Same source" else "Different source",
      Study_relationship = if (metadata$Study_ID[[i]] == metadata$Study_ID[[j]]) "Same study" else "Different study"
    )
  })
}

# Summarise non-independent genome pairs descriptively without naive inferential tests.
summarise_pair_group <- function(data, label, condition) {
  values <- data$AMR_Jaccard_similarity[condition]
  tibble(
    Stratum = label, Pair_count = length(values), Median_similarity = if (length(values)) median(values) else NA_real_,
    IQR_similarity = if (length(values)) IQR(values) else NA_real_,
    Minimum_similarity = if (length(values)) min(values) else NA_real_,
    Maximum_similarity = if (length(values)) max(values) else NA_real_
  )
}

# Write the complete integrated analysis only after every cross-pipeline gate passes.
main <- function() {
  args <- parse_args(); root <- get_project_root()
  paths <- list(
    assignments = file.path(root, "data", "processed", "mlst", "formal_mlst_assignments_72.csv"),
    amr_matrix = file.path(root, "results", "tables", "amr", "amr_gene_presence_absence.csv"),
    functional = file.path(root, "data", "processed", "amr", "amr_functional_hits_with_contigs_72.csv"),
    metadata = file.path(root, "data", "metadata", "sample_metadata_source_72.tsv"),
    proxy_config = file.path(root, "config", "hlgr_proxy_determinants.tsv"),
    tree = file.path(root, "analysis", "phylogenomics", "72_genomes_9_locus_observed_indels", "iqtree", "efaecalis_72_genomes_9_locus_observed_indels.treefile")
  )
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing)) stop("Missing required frozen input(s): ", paste(missing, collapse = ", "), call. = FALSE)
  output_dir <- file.path(root, "results", "tables", "mlst_amr_phylogeny")
  outputs <- set_names(file.path(output_dir, c(
    "integrated_genome_mlst_amr_72.csv", "hlgr_proxy_supporting_hits_72.csv", "hlgr_proxy_review_audit_72.csv",
    "st_hlgr_proxy_prevalence_72.csv", "st_concentration_permutation_exploratory_72.csv",
    "phylogenetic_clustering_permutations_72.csv", "amr_phylogeny_mantel_exploratory_72.csv",
    "amr_pairwise_relationship_audit_72.csv", "amr_pairwise_similarity_summary_72.csv",
    "mlst_amr_phylogeny_qc_summary_72.csv"
  )), c("integrated", "proxy_hits", "proxy_review", "st_summary", "st_permutation", "phylo_permutation", "mantel", "pair_audit", "pair_summary", "qc"))
  preflight_outputs(outputs, args$overwrite)

  assignments <- read_csv(paths$assignments, col_types = cols(.default = col_character()))
  amr <- read_csv(paths$amr_matrix, col_types = cols(Genome = col_character(), Source = col_character(), .default = col_integer()))
  functional <- read_csv(paths$functional, col_types = cols(.default = col_character()))
  metadata <- read_tsv(paths$metadata, col_types = cols(.default = col_character()), na = c("", "NA"))
  proxy_config <- read_tsv(paths$proxy_config, col_types = cols(.default = col_character()))
  validate_genomes(assignments, "Genome", "Formal MLST assignments"); validate_genomes(amr, "Genome", "Frozen AMR matrix")
  if (n_distinct(metadata$assembly_accession) != 72L || nrow(metadata) != 72L) stop("Source metadata must contain 72 unique accessions", call. = FALSE)

  watchlist <- proxy_config |> filter(Decision == "review_required") |> pull(Gene_symbol)
  review_hits <- functional |>
    filter(Gene != PRIMARY_PROXY, toupper(Class) == "AMINOGLYCOSIDE", Gene %in% watchlist | str_detect(toupper(Subclass), "GENTAMICIN")) |>
    mutate(Review_status = "HLGR_proxy_review_required")
  write_csv(review_hits, outputs[["proxy_review"]])
  if (nrow(review_hits)) stop("Accepted non-primary determinant(s) require HLGR-proxy review; see ", outputs[["proxy_review"]], call. = FALSE)

  proxy_hits <- functional |> filter(Gene == PRIMARY_PROXY)
  matrix_carriers <- amr |> filter(.data[[PRIMARY_PROXY]] == 1L) |> pull(Genome)
  hit_carriers <- sort(unique(proxy_hits$Genome))
  if (!identical(sort(matrix_carriers), hit_carriers)) stop("Frozen AMR matrix and accepted-hit carrier accessions disagree", call. = FALSE)
  source_counts <- amr |> filter(Genome %in% matrix_carriers) |> count(Source) |> deframe()
  if (length(matrix_carriers) != 29L || source_counts[["Reproductive"]] != 5L || source_counts[["Bacteraemia"]] != 24L) {
    stop("Expected 29 accepted aac(6')-Ie/aph(2'')-Ia carriers from the validated AMR analysis, but detected ", length(matrix_carriers), ". Input state may have changed. Refusing integration until reviewed.", call. = FALSE)
  }
  write_csv(proxy_hits, outputs[["proxy_hits"]])

  metadata_clean <- metadata |>
    transmute(
      Genome = assembly_accession, Source_metadata = dataset_category, Study_ID = bioproject_accession,
      Study_title = study_title, Associated_publication = associated_paper, Country = country, Year = collection_year
    )
  if (anyNA(metadata_clean$Study_ID) || any(metadata_clean$Study_ID == "")) stop("Study_ID/BioProject must be complete", call. = FALSE)
  burden <- amr |> mutate(AMR_burden = rowSums(across(-c(Genome, Source)))) |> select(Genome, Source_amr = Source, AMR_burden)
  other_amino <- functional |> filter(toupper(Class) == "AMINOGLYCOSIDE", Gene != PRIMARY_PROXY) |>
    distinct(Genome, Gene) |> arrange(Genome, Gene) |> summarise(Other_aminoglycoside_genes = paste(Gene, collapse = "; "), .by = Genome)
  integrated <- assignments |>
    left_join(metadata_clean, by = "Genome") |> left_join(burden, by = "Genome") |> left_join(other_amino, by = "Genome") |>
    mutate(
      HLGR_proxy = if_else(Genome %in% matrix_carriers, "HLGR-associated genotype proxy detected", "HLGR-associated genotype proxy not detected"),
      Other_aminoglycoside_genes = replace_na(Other_aminoglycoside_genes, "None detected")
    )
  if (anyNA(integrated$Study_ID) || any(integrated$Source != integrated$Source_metadata) || any(integrated$Source != integrated$Source_amr)) stop("Metadata/MLST/AMR source join disagreement", call. = FALSE)
  integrated <- integrated |> select(Genome, Source, Study_ID, Study_title, Associated_publication, Country, Year,
                                     all_of(c("gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL")), Allele_profile,
                                     Dataset_profile, Dataset_profile_complete, MLST_status, MLST_status_reason, ST,
                                     PubMLST_CC_if_officially_available, HLGR_proxy, Other_aminoglycoside_genes, AMR_burden) |> arrange(Genome)
  validate_genomes(integrated, "Genome", "Integrated table")

  assigned_carriers <- sum(integrated$ST != "Unassigned" & integrated$HLGR_proxy == "HLGR-associated genotype proxy detected")
  st_summary <- integrated |> filter(ST != "Unassigned") |>
    summarise(
      Genomes = n(), Source_composition = paste(names(table(Source)), as.integer(table(Source)), sep = "=", collapse = "; "),
      Study_composition = paste(names(table(Study_ID)), as.integer(table(Study_ID)), sep = "=", collapse = "; "),
      Proxy_positive = sum(HLGR_proxy == "HLGR-associated genotype proxy detected"),
      Proxy_negative = sum(HLGR_proxy != "HLGR-associated genotype proxy detected"),
      Within_ST_proxy_prevalence_percent = round(100 * Proxy_positive / Genomes, 2),
      Percent_of_formally_assigned_carriers = round(100 * Proxy_positive / assigned_carriers, 2), .by = ST
    ) |> arrange(desc(Proxy_positive), desc(Genomes), as.integer(ST))
  st_permutations <- bind_rows(
    run_group_permutation(integrated$ST, integrated$HLGR_proxy == "HLGR-associated genotype proxy detected", "Formal PubMLST ST concentration"),
    run_group_permutation(if_else(integrated$Dataset_profile_complete == "TRUE", integrated$Dataset_profile, NA_character_), integrated$HLGR_proxy == "HLGR-associated genotype proxy detected", "Complete seven-locus dataset-profile sensitivity")
  )

  tree <- read.tree(paths$tree)
  if (length(tree$tip.label) != 72L || n_distinct(tree$tip.label) != 72L || !setequal(tree$tip.label, integrated$Genome)) stop("Frozen unrooted tree tips do not match 72 integrated genomes", call. = FALSE)
  integrated_tree_order <- integrated |> slice(match(tree$tip.label, Genome))
  phylo_distance <- cophenetic.phylo(tree)
  phylo_distance <- phylo_distance[integrated_tree_order$Genome, integrated_tree_order$Genome]
  phylo_permutations <- run_phylogenetic_permutations(phylo_distance, integrated_tree_order)

  amr_tree_order <- amr |> slice(match(tree$tip.label, Genome))
  binary <- as.matrix(amr_tree_order |> select(-Genome, -Source)) > 0
  rownames(binary) <- amr_tree_order$Genome
  amr_distance <- jaccard_matrix(binary)
  mantel <- run_mantel(amr_distance, phylo_distance)
  pair_audit <- build_pairwise_audit(integrated_tree_order, amr_distance, phylo_distance)
  pair_summary <- bind_rows(
    summarise_pair_group(pair_audit, "Same ST", pair_audit$ST_relationship == "Same ST"),
    summarise_pair_group(pair_audit, "Different ST", pair_audit$ST_relationship == "Different ST"),
    summarise_pair_group(pair_audit, "Same source", pair_audit$Source_relationship == "Same source"),
    summarise_pair_group(pair_audit, "Different source", pair_audit$Source_relationship == "Different source"),
    summarise_pair_group(pair_audit, "Same ST + same source", pair_audit$ST_relationship == "Same ST" & pair_audit$Source_relationship == "Same source"),
    summarise_pair_group(pair_audit, "Same ST + different source", pair_audit$ST_relationship == "Same ST" & pair_audit$Source_relationship == "Different source"),
    summarise_pair_group(pair_audit, "Different ST + same source", pair_audit$ST_relationship == "Different ST" & pair_audit$Source_relationship == "Same source"),
    summarise_pair_group(pair_audit, "Different ST + different source", pair_audit$ST_relationship == "Different ST" & pair_audit$Source_relationship == "Different source"),
    summarise_pair_group(pair_audit, "Same study", pair_audit$Study_relationship == "Same study"),
    summarise_pair_group(pair_audit, "Different study", pair_audit$Study_relationship == "Different study")
  ) |> mutate(Note = "Descriptive only: genome pairs are non-independent and ST, source, and study overlap")
  qc <- tibble(
    Metric = c("genomes", "reproductive_genomes", "bacteraemia_genomes", "formal_st_assigned", "formal_st_unassigned",
               "complete_dataset_profiles", "proxy_positive_total", "proxy_positive_reproductive", "proxy_positive_bacteraemia",
               "tree_tips_matched", "amr_genes", "raw_amrfinder_rows_reinterpreted"),
    Value = as.character(c(72, 14, 58, sum(integrated$ST != "Unassigned"), sum(integrated$ST == "Unassigned"),
                           sum(integrated$Dataset_profile_complete == "TRUE"), 29, 5, 24, 72, ncol(binary), 0))
  )
  write_csv(integrated, outputs[["integrated"]]); write_csv(st_summary, outputs[["st_summary"]]); write_csv(st_permutations, outputs[["st_permutation"]])
  write_csv(phylo_permutations, outputs[["phylo_permutation"]]); write_csv(mantel, outputs[["mantel"]]); write_csv(pair_audit, outputs[["pair_audit"]])
  write_csv(pair_summary, outputs[["pair_summary"]]); write_csv(qc, outputs[["qc"]])
  print(qc, n = Inf)
}

# Run only during a deliberate Rscript invocation.
if (!interactive() && Sys.getenv("EFAECALIS_SKIP_MAIN", "0") != "1") main()
