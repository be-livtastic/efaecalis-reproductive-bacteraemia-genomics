# ================================================================
# AMRFinderPlus analysis for E. faecalis dissertation datasets
#
# Outputs:
#   - Copies of raw TSV files
#   - CSV conversions of raw TSV files
#   - Cleaned AMRFinderPlus results
#   - Confidence-separated hit tables
#   - Low-coverage hit tables (<80% reference coverage)
#   - Per-genome AMR summaries
#   - Drug-class, dataset-comparison, and selected-category plots
# ================================================================

# ------------------------------
# 1. Install/load packages
# ------------------------------
required_packages <- c(
  "tidyverse",
  "janitor",
  "fs",
  "scales"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ------------------------------
# 2. Project settings
# ------------------------------

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath(
    file.path(getwd(), "scripts", "03_amr_analysis",
              "amrfinder_analysis_72_genomes.R"),
    mustWork = TRUE
  )
}

default_project_dir <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)
PROJECT_DIR <- Sys.getenv(
  "EFAECALIS_PROJECT_ROOT",
  unset = default_project_dir
)
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

ALLOW_OVERWRITE <- tolower(
  Sys.getenv("EFAECALIS_ALLOW_OVERWRITE", unset = "false")
) %in% c("1", "true", "yes")

REPRODUCTIVE_DIR <- file.path(
  PROJECT_DIR, "data", "interim", "amrfinder", "reproductive"
)
BACTERAEMIA_DIR <- file.path(
  PROJECT_DIR, "data", "interim", "amrfinder", "bacteraemia"
)
OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "amr")

# Coverage threshold requested for flagging low-coverage AMR hits.
LOW_COVERAGE_THRESHOLD <- 80

# Optional: specify the cleaned name of a column containing genome/isolate IDs.
# Leave as NULL to let the script detect a likely column automatically.
# Examples after janitor::clean_names():
#   "genome", "genome_id", "isolate", "assembly_accession", "sample"
GENOME_ID_COLUMN <- NULL

# Refuse to overwrite a previous run unless explicitly enabled. This check is
# deliberately performed before any output directory is created.
if (dir.exists(OUTPUT_DIR) &&
    length(list.files(OUTPUT_DIR, all.files = TRUE, no.. = TRUE)) > 0 &&
    !ALLOW_OVERWRITE) {
  stop(
    "Output directory is not empty: ", OUTPUT_DIR, "\n",
    "Choose a new output directory or set ",
    "EFAECALIS_ALLOW_OVERWRITE=true after reviewing its contents."
  )
}

# Create organised output folders safely.
output_folders <- c(
  OUTPUT_DIR,
  file.path(OUTPUT_DIR, "raw_tsv"),
  file.path(OUTPUT_DIR, "raw_csv"),
  file.path(OUTPUT_DIR, "cleaned_tables"),
  file.path(OUTPUT_DIR, "summary_tables"),
  file.path(OUTPUT_DIR, "graphs")
)

purrr::walk(output_folders, fs::dir_create)

# Stop early if input folders are not found.
if (!dir.exists(REPRODUCTIVE_DIR)) {
  stop("Reproductive dataset folder not found: ", REPRODUCTIVE_DIR)
}

if (!dir.exists(BACTERAEMIA_DIR)) {
  stop("Bacteraemia dataset folder not found: ", BACTERAEMIA_DIR)
}

# ------------------------------
# 3. Helper functions
# ------------------------------

# Return the first available column from a list of possible names.
first_existing_column <- function(data, candidates) {
  match <- candidates[candidates %in% names(data)]
  if (length(match) == 0) return(NULL)
  match[[1]]
}

# Convert values such as "98.6", "98.6%", or numeric values to numeric.
to_numeric_percent <- function(x) {
  readr::parse_number(as.character(x))
}

# Create a readable genome ID from a filename.
genome_id_from_filename <- function(path) {
  basename(path) |>
    tools::file_path_sans_ext() |>
    stringr::str_remove("(?i)([_-]?(amrfinderplus|amrfinder|results?|output|hits?|all))+$")
}

# Identify a likely genome/isolate column.
detect_genome_column <- function(data) {
  if (!is.null(GENOME_ID_COLUMN) && GENOME_ID_COLUMN %in% names(data)) {
    return(GENOME_ID_COLUMN)
  }
  
  first_existing_column(
    data,
    c(
      "genome_id", "genome", "isolate_id", "isolate", "assembly_accession",
      "assembly", "biosample_accession", "biosample", "sample_id", "sample",
      "strain", "file_name", "filename", "source_file"
    )
  )
}

# Read every TSV in one dataset folder, copy the raw files, and add metadata.
read_amrfinder_dataset <- function(input_dir, dataset_name) {
  tsv_files <- fs::dir_ls(
    input_dir,
    recurse = TRUE,
    type = "file",
    regexp = "(?i)\\.tsv$"
  )
  
  if (length(tsv_files) == 0) {
    stop("No TSV files were found in: ", input_dir)
  }
  
  message("Reading ", length(tsv_files), " TSV file(s) for ", dataset_name, ".")
  
  purrr::map_dfr(tsv_files, function(file_path) {
    # Preserve an untouched copy of each raw TSV.
    copied_name <- paste0(dataset_name, "__", basename(file_path))
    fs::file_copy(
      file_path,
      file.path(OUTPUT_DIR, "raw_tsv", copied_name),
      overwrite = ALLOW_OVERWRITE
    )
    
    raw_data <- readr::read_tsv(
      file_path,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    )
    
    # Convert the original content to CSV before cleaning column names.
    readr::write_csv(
      raw_data,
      file.path(
        OUTPUT_DIR,
        "raw_csv",
        paste0(dataset_name, "__", tools::file_path_sans_ext(basename(file_path)), ".csv")
      ),
      na = ""
    )
    
    cleaned <- raw_data |>
      janitor::clean_names()
    
    genome_column <- detect_genome_column(cleaned)
    fallback_id <- genome_id_from_filename(file_path)
    
    if (is.null(genome_column)) {
      cleaned <- cleaned |>
        mutate(genome_id = fallback_id)
    } else {
      cleaned <- cleaned |>
        mutate(
          genome_id = as.character(.data[[genome_column]]),
          genome_id = na_if(genome_id, ""),
          genome_id = replace_na(genome_id, fallback_id)
        )
    }
    
    cleaned |>
      mutate(
        dataset = dataset_name,
        source_file = basename(file_path),
        source_path = as.character(file_path),
        .before = 1
      )
  })
}

# Standardise the important AMRFinderPlus columns while retaining all raw fields.
standardise_amrfinder_columns <- function(data) {
  gene_col <- first_existing_column(
    data,
    c("gene_symbol", "gene", "name", "sequence_name")
  )
  
  sequence_name_col <- first_existing_column(
    data,
    c("sequence_name", "name", "name_of_closest_sequence", "gene_symbol")
  )
  
  class_col <- first_existing_column(
    data,
    c("class", "drug_class", "antimicrobial_class")
  )
  
  subclass_col <- first_existing_column(
    data,
    c("subclass", "drug_subclass", "antimicrobial_subclass")
  )
  
  method_col <- first_existing_column(
    data,
    c("method", "detection_method")
  )
  
  coverage_col <- first_existing_column(
    data,
    c(
      "coverage_of_reference_sequence",
      "percent_coverage_of_reference_sequence",
      "reference_coverage",
      "coverage",
      "pct_coverage_of_reference_sequence"
    )
  )
  
  identity_col <- first_existing_column(
    data,
    c(
      "identity_to_reference_sequence",
      "percent_identity_to_reference_sequence",
      "reference_identity",
      "identity",
      "pct_identity_to_reference_sequence"
    )
  )
  
  element_type_col <- first_existing_column(data, c("element_type", "type"))
  element_subtype_col <- first_existing_column(data, c("element_subtype", "subtype"))
  contig_col <- first_existing_column(data, c("contig_id", "contig", "sequence_id"))
  
  # Warn clearly when essential columns are missing.
  missing_essential <- c(
    gene = is.null(gene_col),
    method = is.null(method_col),
    coverage = is.null(coverage_col)
  )
  
  if (any(missing_essential)) {
    warning(
      "Could not automatically detect these columns: ",
      paste(names(missing_essential)[missing_essential], collapse = ", "),
      ". Check names(combined_raw) and update the candidate lists if needed."
    )
  }
  
  data |>
    mutate(
      gene_symbol_std = if (!is.null(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
      sequence_name_std = if (!is.null(sequence_name_col)) as.character(.data[[sequence_name_col]]) else NA_character_,
      drug_class_std = if (!is.null(class_col)) as.character(.data[[class_col]]) else NA_character_,
      drug_subclass_std = if (!is.null(subclass_col)) as.character(.data[[subclass_col]]) else NA_character_,
      method_std = if (!is.null(method_col)) stringr::str_to_upper(as.character(.data[[method_col]])) else NA_character_,
      reference_coverage_pct = if (!is.null(coverage_col)) to_numeric_percent(.data[[coverage_col]]) else NA_real_,
      reference_identity_pct = if (!is.null(identity_col)) to_numeric_percent(.data[[identity_col]]) else NA_real_,
      element_type_std = if (!is.null(element_type_col)) as.character(.data[[element_type_col]]) else NA_character_,
      element_subtype_std = if (!is.null(element_subtype_col)) as.character(.data[[element_subtype_col]]) else NA_character_,
      contig_id_std = if (!is.null(contig_col)) as.character(.data[[contig_col]]) else NA_character_,
      gene_label = dplyr::coalesce(gene_symbol_std, sequence_name_std, "Unlabelled hit")
    )
}

# Assign confidence and functional-status labels.
classify_hits <- function(data) {
  data |>
    mutate(
      method_std = replace_na(method_std, "UNKNOWN"),
      
      confidence_group = case_when(
        stringr::str_detect(method_std, "INTERNAL[_ ]?STOP") ~ "INTERNAL_STOP",
        method_std %in% c("EXACT", "EXACTX") |
          stringr::str_detect(method_std, "^EXACTX?$") ~ "EXACT / EXACTX",
        method_std %in% c("BLAST", "BLASTX") |
          stringr::str_detect(method_std, "^BLASTX?$") ~ "BLASTX",
        method_std %in% c("PARTIAL", "PARTIALX") |
          stringr::str_detect(method_std, "PARTIAL") ~ "PARTIALX",
        TRUE ~ "OTHER / UNKNOWN"
      ),
      
      functional_amr_hit = confidence_group != "INTERNAL_STOP",
      
      coverage_flag = case_when(
        is.na(reference_coverage_pct) ~ "Coverage unavailable",
        reference_coverage_pct < LOW_COVERAGE_THRESHOLD ~ paste0(
          "Low coverage (<", LOW_COVERAGE_THRESHOLD, "%)"
        ),
        TRUE ~ paste0("Coverage >=", LOW_COVERAGE_THRESHOLD, "%")
      ),
      
      low_coverage = !is.na(reference_coverage_pct) &
        reference_coverage_pct < LOW_COVERAGE_THRESHOLD
    )
}

# Add biologically relevant category flags.
add_gene_category_flags <- function(data) {
  data |>
    mutate(
      searchable_annotation = stringr::str_to_lower(
        paste(
          replace_na(gene_symbol_std, ""),
          replace_na(sequence_name_std, ""),
          replace_na(drug_class_std, ""),
          replace_na(drug_subclass_std, "")
        )
      ),
      
      streptogramin_lincosamide_flag = stringr::str_detect(
        searchable_annotation,
        paste0(
          "streptogramin|lincosamide|lincomycin|clindamycin|",
          "(^|[^a-z])lsa[a-z0-9()'_-]*|(^|[^a-z])lnu[a-z0-9()'_-]*|",
          "(^|[^a-z])vat[a-z0-9()'_-]*|(^|[^a-z])vga[a-z0-9()'_-]*"
        )
      ),
      
      tetracycline_flag = stringr::str_detect(
        searchable_annotation,
        "tetracycline|(^|[^a-z])tet[a-z0-9()'_-]*"
      ),
      
      aminoglycoside_flag = stringr::str_detect(
        searchable_annotation,
        paste0(
          "aminoglycoside|gentamicin|kanamycin|amikacin|tobramycin|streptomycin|",
          "(^|[^a-z])aac[a-z0-9()'()_-]*|(^|[^a-z])aph[a-z0-9()'()_-]*|",
          "(^|[^a-z])ant[a-z0-9()'()_-]*|(^|[^a-z])aad[a-z0-9()'()_-]*"
        )
      ),
      
      gentamicin_related_flag = stringr::str_detect(
        searchable_annotation,
        paste0(
          "gentamicin|aac\\(6['’′]?\\)[- ]?ie[- ]?aph\\(2['’′]{2}?\\)[- ]?ia|",
          "aac\\(6['’′]?\\)[- ]?ii|aph\\(2['’′]{2}?\\)|aac6|aph2"
        )
      )
    )
}

# Collapse distinct, non-missing values into a semicolon-separated string.
collapse_distinct <- function(x) {
  values <- x |>
    as.character() |>
    stringr::str_trim() |>
    na_if("") |>
    stats::na.omit() |>
    unique() |>
    sort()
  
  if (length(values) == 0) return(NA_character_)
  paste(values, collapse = "; ")
}

# Collapse genes only where a selected logical flag is TRUE.
collapse_flagged_genes <- function(gene, flag) {
  collapse_distinct(gene[replace_na(flag, FALSE)])
}

# Create one row per genome.
create_genome_summary <- function(data) {
  functional <- data |>
    filter(functional_amr_hit)
  
  all_genomes <- data |>
    distinct(dataset, genome_id)
  
  summaries <- functional |>
    group_by(dataset, genome_id) |>
    summarise(
      functional_amr_hit_count = n(),
      unique_amr_gene_count = n_distinct(gene_label[!is.na(gene_label)]),
      amr_genes_present = collapse_distinct(gene_label),
      amr_drug_classes = collapse_distinct(drug_class_std),
      amr_drug_subclasses = collapse_distinct(drug_subclass_std),
      
      exact_exactx_hit_count = sum(confidence_group == "EXACT / EXACTX", na.rm = TRUE),
      blastx_hit_count = sum(confidence_group == "BLASTX", na.rm = TRUE),
      partialx_hit_count = sum(confidence_group == "PARTIALX", na.rm = TRUE),
      other_hit_count = sum(confidence_group == "OTHER / UNKNOWN", na.rm = TRUE),
      low_coverage_functional_hit_count = sum(low_coverage, na.rm = TRUE),
      
      streptogramin_lincosamide_gene_count = n_distinct(
        gene_label[streptogramin_lincosamide_flag & !is.na(gene_label)]
      ),
      streptogramin_lincosamide_genes = collapse_flagged_genes(
        gene_label, streptogramin_lincosamide_flag
      ),
      
      tetracycline_gene_count = n_distinct(
        gene_label[tetracycline_flag & !is.na(gene_label)]
      ),
      tetracycline_genes = collapse_flagged_genes(gene_label, tetracycline_flag),
      
      aminoglycoside_gene_count = n_distinct(
        gene_label[aminoglycoside_flag & !is.na(gene_label)]
      ),
      aminoglycoside_genes = collapse_flagged_genes(gene_label, aminoglycoside_flag),
      
      gentamicin_related_gene_count = n_distinct(
        gene_label[gentamicin_related_flag & !is.na(gene_label)]
      ),
      gentamicin_related_genes = collapse_flagged_genes(
        gene_label, gentamicin_related_flag
      ),
      
      .groups = "drop"
    )
  
  internal_stop_counts <- data |>
    group_by(dataset, genome_id) |>
    summarise(
      internal_stop_hit_count = sum(confidence_group == "INTERNAL_STOP", na.rm = TRUE),
      internal_stop_genes = collapse_flagged_genes(
        gene_label,
        confidence_group == "INTERNAL_STOP"
      ),
      .groups = "drop"
    )
  
  all_genomes |>
    left_join(summaries, by = c("dataset", "genome_id")) |>
    left_join(internal_stop_counts, by = c("dataset", "genome_id")) |>
    mutate(
      across(
        ends_with("_count"),
        ~ replace_na(as.integer(.x), 0L)
      )
    ) |>
    arrange(dataset, desc(unique_amr_gene_count), genome_id)
}

# Save a dataframe as CSV.
save_csv <- function(data, filename, folder = "cleaned_tables") {
  readr::write_csv(
    data,
    file.path(OUTPUT_DIR, folder, filename),
    na = ""
  )
}

# ------------------------------
# 4. Import both datasets
# ------------------------------
reproductive_raw <- read_amrfinder_dataset(
  REPRODUCTIVE_DIR,
  "Reproductive"
)

bacteraemia_raw <- read_amrfinder_dataset(
  BACTERAEMIA_DIR,
  "Bacteraemia"
)

combined_raw <- bind_rows(reproductive_raw, bacteraemia_raw)

# ------------------------------
# 5. Clean and classify results
# ------------------------------
combined_clean <- combined_raw |>
  standardise_amrfinder_columns() |>
  classify_hits() |>
  add_gene_category_flags()

reproductive_clean <- combined_clean |>
  filter(dataset == "Reproductive")

bacteraemia_clean <- combined_clean |>
  filter(dataset == "Bacteraemia")

# Main confidence groups.
exact_exactx_hits <- combined_clean |>
  filter(confidence_group == "EXACT / EXACTX")

blastx_hits <- combined_clean |>
  filter(confidence_group == "BLASTX")

partialx_hits <- combined_clean |>
  filter(confidence_group == "PARTIALX")

internal_stop_hits <- combined_clean |>
  filter(confidence_group == "INTERNAL_STOP")

low_coverage_hits <- combined_clean |>
  filter(low_coverage)

# Functional AMR table excludes INTERNAL_STOP hits.
functional_amr_hits <- combined_clean |>
  filter(functional_amr_hit)

# ------------------------------
# 6. Save cleaned and separated hit tables
# ------------------------------
save_csv(combined_clean, "combined_72_amrfinder_cleaned.csv")
save_csv(reproductive_clean, "reproductive_14_amrfinder_cleaned.csv")
save_csv(bacteraemia_clean, "bacteraemia_58_amrfinder_cleaned.csv")

save_csv(functional_amr_hits, "combined_72_functional_amr_hits_excluding_internal_stop.csv")
save_csv(exact_exactx_hits, "combined_72_exact_exactx_hits.csv")
save_csv(blastx_hits, "combined_72_blastx_hits.csv")
save_csv(partialx_hits, "combined_72_partialx_hits.csv")
save_csv(internal_stop_hits, "combined_72_internal_stop_hits_not_counted_as_functional.csv")
save_csv(low_coverage_hits, paste0("combined_72_low_coverage_hits_below_", LOW_COVERAGE_THRESHOLD, "pct.csv"))

# Also save dataset-specific confidence/coverage tables.
for (dataset_name in c("Reproductive", "Bacteraemia")) {
  dataset_slug <- stringr::str_to_lower(dataset_name)
  dataset_data <- combined_clean |>
    filter(dataset == dataset_name)
  
  save_csv(
    dataset_data |> filter(functional_amr_hit),
    paste0(dataset_slug, "_functional_amr_hits.csv")
  )
  
  save_csv(
    dataset_data |> filter(confidence_group == "INTERNAL_STOP"),
    paste0(dataset_slug, "_internal_stop_hits.csv")
  )
  
  save_csv(
    dataset_data |> filter(low_coverage),
    paste0(dataset_slug, "_low_coverage_hits_below_", LOW_COVERAGE_THRESHOLD, "pct.csv")
  )
}

# ------------------------------
# 7. Create AMR summary tables
# ------------------------------
combined_summary <- create_genome_summary(combined_clean)

reproductive_summary <- combined_summary |>
  filter(dataset == "Reproductive")

bacteraemia_summary <- combined_summary |>
  filter(dataset == "Bacteraemia")

save_csv(
  reproductive_summary,
  "reproductive_14_amr_summary_by_genome.csv",
  folder = "summary_tables"
)

save_csv(
  bacteraemia_summary,
  "bacteraemia_58_amr_summary_by_genome.csv",
  folder = "summary_tables"
)

save_csv(
  combined_summary,
  "combined_72_amr_summary_by_genome.csv",
  folder = "summary_tables"
)

# Dataset-level overview table.
dataset_overview <- combined_summary |>
  group_by(dataset) |>
  summarise(
    genomes_in_summary = n(),
    genomes_with_functional_amr_hits = sum(unique_amr_gene_count > 0),
    median_unique_amr_genes = median(unique_amr_gene_count),
    mean_unique_amr_genes = mean(unique_amr_gene_count),
    minimum_unique_amr_genes = min(unique_amr_gene_count),
    maximum_unique_amr_genes = max(unique_amr_gene_count),
    genomes_with_low_coverage_hits = sum(low_coverage_functional_hit_count > 0),
    genomes_with_internal_stop_hits = sum(internal_stop_hit_count > 0),
    genomes_with_streptogramin_lincosamide_genes = sum(
      streptogramin_lincosamide_gene_count > 0
    ),
    genomes_with_tetracycline_genes = sum(tetracycline_gene_count > 0),
    genomes_with_aminoglycoside_genes = sum(aminoglycoside_gene_count > 0),
    genomes_with_gentamicin_related_genes = sum(gentamicin_related_gene_count > 0),
    .groups = "drop"
  )

save_csv(
  dataset_overview,
  "dataset_level_amr_overview.csv",
  folder = "summary_tables"
)

# Gene prevalence across genomes.
gene_prevalence <- functional_amr_hits |>
  distinct(dataset, genome_id, gene_label) |>
  count(dataset, gene_label, name = "genomes_with_gene") |>
  group_by(dataset) |>
  mutate(
    total_genomes = n_distinct(combined_summary$genome_id[combined_summary$dataset == first(dataset)]),
    prevalence_pct = 100 * genomes_with_gene / total_genomes
  ) |>
  ungroup() |>
  arrange(dataset, desc(genomes_with_gene), gene_label)

save_csv(
  gene_prevalence,
  "amr_gene_prevalence_by_dataset.csv",
  folder = "summary_tables"
)

# Drug-class prevalence across genomes.
drug_class_prevalence <- functional_amr_hits |>
  filter(!is.na(drug_class_std), drug_class_std != "") |>
  distinct(dataset, genome_id, drug_class_std) |>
  count(dataset, drug_class_std, name = "genomes_with_class") |>
  group_by(dataset) |>
  mutate(
    total_genomes = n_distinct(combined_summary$genome_id[combined_summary$dataset == first(dataset)]),
    prevalence_pct = 100 * genomes_with_class / total_genomes
  ) |>
  ungroup() |>
  arrange(dataset, desc(genomes_with_class), drug_class_std)

save_csv(
  drug_class_prevalence,
  "amr_drug_class_prevalence_by_dataset.csv",
  folder = "summary_tables"
)

# ------------------------------
# 8. Graphing functions
# ------------------------------
save_plot <- function(plot_object, filename, width = 11, height = 7) {
  ggplot2::ggsave(
    filename = file.path(OUTPUT_DIR, "graphs", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

plot_drug_classes <- function(data, dataset_label, filename) {
  plot_data <- data |>
    filter(
      functional_amr_hit,
      !is.na(drug_class_std),
      drug_class_std != ""
    ) |>
    distinct(genome_id, drug_class_std) |>
    count(drug_class_std, sort = TRUE, name = "genome_count") |>
    mutate(drug_class_std = forcats::fct_reorder(drug_class_std, genome_count))
  
  p <- ggplot(plot_data, aes(x = drug_class_std, y = genome_count)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("AMR drug-class distribution: ", dataset_label),
      subtitle = "Each drug class is counted once per genome",
      x = "AMR drug class",
      y = "Number of genomes"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  
  save_plot(p, filename)
}

# ------------------------------
# 9. Create retained graphs
# ------------------------------
plot_drug_classes(
  reproductive_clean,
  "14 reproductive genomes",
  "reproductive_14_amr_drug_classes.png"
)

plot_drug_classes(
  bacteraemia_clean,
  "58 bacteraemia genomes",
  "bacteraemia_58_amr_drug_classes.png"
)

plot_drug_classes(
  combined_clean,
  "combined 72 genomes",
  "combined_72_amr_drug_classes.png"
)

# Compare unique AMR gene burden between the two source groups.
p_comparison <- combined_summary |>
  ggplot(aes(x = dataset, y = unique_amr_gene_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.65) +
  labs(
    title = "Functional AMR gene burden by dataset",
    subtitle = "Each point represents one genome; INTERNAL_STOP hits are excluded",
    x = "Dataset",
    y = "Number of unique functional AMR genes per genome"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

save_plot(
  p_comparison,
  "combined_72_amr_gene_burden_dataset_comparison.png",
  width = 8,
  height = 6
)

# Presence of the requested AMR gene categories by dataset.
category_presence <- combined_summary |>
  transmute(
    dataset,
    genome_id,
    `Streptogramin / lincosamide` = streptogramin_lincosamide_gene_count > 0,
    Tetracycline = tetracycline_gene_count > 0,
    Aminoglycoside = aminoglycoside_gene_count > 0,
    `Gentamicin-related` = gentamicin_related_gene_count > 0
  ) |>
  pivot_longer(
    cols = -c(dataset, genome_id),
    names_to = "gene_category",
    values_to = "present"
  ) |>
  group_by(dataset, gene_category) |>
  summarise(
    genomes_with_category = sum(present),
    total_genomes = n(),
    prevalence_pct = 100 * genomes_with_category / total_genomes,
    .groups = "drop"
  )

save_csv(
  category_presence,
  "requested_gene_category_prevalence.csv",
  folder = "summary_tables"
)

p_categories <- ggplot(
  category_presence,
  aes(x = gene_category, y = prevalence_pct, fill = dataset)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Prevalence of selected AMR gene categories",
    subtitle = "Percentage of genomes containing at least one functional gene in each category",
    x = "AMR gene category",
    y = "Genomes containing category (%)",
    fill = "Dataset"
  ) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank()
  )

save_plot(
  p_categories,
  "combined_72_selected_amr_category_prevalence.png",
  width = 10,
  height = 7
)

# ------------------------------
# 10. Quality-control report
# ------------------------------
qc_report <- tibble(
  check = c(
    "Rows imported",
    "Reproductive rows imported",
    "Bacteraemia rows imported",
    "Unique genome IDs detected",
    "Reproductive genome IDs detected",
    "Bacteraemia genome IDs detected",
    "Functional AMR hits",
    "INTERNAL_STOP hits excluded from functional counts",
    paste0("Low-coverage hits below ", LOW_COVERAGE_THRESHOLD, "%"),
    "Rows with missing method",
    "Rows with missing reference coverage",
    "Rows with missing gene label"
  ),
  value = c(
    nrow(combined_clean),
    nrow(reproductive_clean),
    nrow(bacteraemia_clean),
    n_distinct(combined_clean$genome_id),
    n_distinct(reproductive_clean$genome_id),
    n_distinct(bacteraemia_clean$genome_id),
    sum(combined_clean$functional_amr_hit, na.rm = TRUE),
    sum(combined_clean$confidence_group == "INTERNAL_STOP", na.rm = TRUE),
    sum(combined_clean$low_coverage, na.rm = TRUE),
    sum(combined_clean$method_std == "UNKNOWN", na.rm = TRUE),
    sum(is.na(combined_clean$reference_coverage_pct)),
    sum(is.na(combined_clean$gene_label) | combined_clean$gene_label == "Unlabelled hit")
  )
)

save_csv(
  qc_report,
  "amrfinder_analysis_quality_control_report.csv",
  folder = "summary_tables"
)

# Warn if the expected 14 + 58 genome IDs were not detected.
reproductive_genome_n <- n_distinct(reproductive_clean$genome_id)
bacteraemia_genome_n <- n_distinct(bacteraemia_clean$genome_id)

if (reproductive_genome_n != 14) {
  warning(
    "Detected ", reproductive_genome_n,
    " reproductive genome ID(s), not 14. Check GENOME_ID_COLUMN or input filenames."
  )
}

if (bacteraemia_genome_n != 58) {
  warning(
    "Detected ", bacteraemia_genome_n,
    " bacteraemia genome ID(s), not 58. Check GENOME_ID_COLUMN or input filenames."
  )
}

# Save package and R session details for reproducibility.
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "amrfinder_analysis_session_info.txt")
)

message("\nAnalysis complete.")
message("Output folder: ", OUTPUT_DIR)
message("Detected reproductive genomes: ", reproductive_genome_n)
message("Detected bacteraemia genomes: ", bacteraemia_genome_n)
message("Detected combined genomes: ", n_distinct(combined_clean$genome_id))
