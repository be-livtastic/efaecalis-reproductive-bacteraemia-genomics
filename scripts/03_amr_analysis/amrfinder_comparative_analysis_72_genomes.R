#!/usr/bin/env Rscript

# ==============================================================================
# Comparative AMRFinderPlus analysis for 72 Enterococcus faecalis genomes
# ==============================================================================
#
# Statistics:
#   Two-sided Fisher's exact tests are applied to genes occurring in at least
#   two genomes. Benjamini-Hochberg correction is applied across tested genes.
# ==============================================================================

EXPECTED_AMRFINDER_COLUMNS <- c(
  "Name", "Protein id", "Contig id", "Start", "Stop", "Strand",
  "Element symbol", "Element name", "Scope", "Type", "Subtype", "Class",
  "Subclass", "Method", "Target length", "Reference sequence length",
  "% Coverage of reference", "% Identity to reference", "Alignment length",
  "Closest reference accession", "Closest reference name", "HMM accession",
  "HMM description"
)

EXPECTED_METHODS <- c("EXACTX", "BLASTX", "POINTX", "PARTIALX", "INTERNAL_STOP")
ACCEPTED_METHODS <- c("EXACTX", "BLASTX", "POINTX")
EXCLUDED_METHODS <- c("PARTIALX", "INTERNAL_STOP")
EXPECTED_GROUP_COUNTS <- c(Reproductive = 14L, Bacteraemia = 58L)
EXPECTED_RAW_ROWS <- 567L
EXPECTED_TYPE_COUNTS <- c(AMR = 523L, STRESS = 44L)
LOW_COVERAGE_THRESHOLD <- 80

# Heatmap display choice only; it does not filter analytical tables.
HEATMAP_TOP_N_GENES <- 30L

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2"
)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
  }

  candidate <- file.path(
    getwd(), "scripts", "03_amr_analysis",
    "amrfinder_comparative_analysis_72_genomes.R"
  )
  if (!file.exists(candidate)) {
    stop(
      "Cannot determine the script location. Run the script directly with Rscript ",
      "or set EFAECALIS_PROJECT_ROOT.",
      call. = FALSE
    )
  }
  normalizePath(candidate, mustWork = TRUE)
}

get_config <- function() {
  script_path <- get_script_path()
  derived_root <- normalizePath(
    file.path(dirname(script_path), "..", ".."),
    mustWork = TRUE
  )
  project_root <- Sys.getenv("EFAECALIS_PROJECT_ROOT", unset = derived_root)
  project_root <- normalizePath(project_root, mustWork = TRUE)

  list(
    project_root = project_root,
    reproductive_amr_dir = file.path(
      project_root, "local_archive", "large_outputs",
      "Efaecalis_14_AMRFinder", "05_amrfinder_results"
    ),
    bacteraemia_amr_dir = file.path(
      project_root, "local_archive", "large_outputs",
      "Efaecalis_58_AMRFinder", "05_amrfinder_results"
    ),
    metadata_path = file.path(
      project_root, "data", "accession_lists", "selected_72_accessions.tsv"
    ),
    strain_metadata_path = file.path(
      project_root, "data", "metadata", "curated_metadata_72_genomes.csv"
    ),
    processed_output_dir = file.path(project_root, "data", "processed", "amr"),
    table_output_dir = file.path(project_root, "results", "tables", "amr"),
    figure_output_dir = file.path(project_root, "results", "figures", "amr"),
    log_output_dir = file.path(project_root, "analysis", "logs", "amr")
  )
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  allowed <- c("--overwrite", "--heatmaps-only", "--help")
  unknown <- setdiff(args, allowed)
  if (length(unknown) > 0L) {
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if ("--help" %in% args) {
    cat(
      "Usage: Rscript scripts/03_amr_analysis/",
      "amrfinder_comparative_analysis_72_genomes.R ",
      "[--overwrite] [--heatmaps-only]\n",
      sep = ""
    )
    quit(save = "no", status = 0L)
  }
  list(
    overwrite = "--overwrite" %in% args,
    heatmaps_only = "--heatmaps-only" %in% args
  )
}

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Required R package(s) are not installed: ", paste(missing, collapse = ", "),
      ". Install them in the project environment before running this pipeline.",
      call. = FALSE
    )
  }
}

output_registry <- function(config) {
  c(
    amr_standardised_long_72 = file.path(
      config$processed_output_dir, "amr_standardised_long_72.csv"
    ),
    amr_flagged_low_confidence_72 = file.path(
      config$processed_output_dir, "amr_flagged_low_confidence_72.csv"
    ),
    amr_method_audit_72 = file.path(
      config$processed_output_dir, "amr_method_audit_72.csv"
    ),
    amr_functional_hits_with_contigs_72 = file.path(
      config$processed_output_dir, "amr_functional_hits_with_contigs_72.csv"
    ),
    amr_gene_presence_absence = file.path(
      config$table_output_dir, "amr_gene_presence_absence.csv"
    ),
    amr_shared_and_dataset_only_determinants = file.path(
      config$table_output_dir, "amr_shared_and_dataset_only_determinants.csv"
    ),
    amr_gene_prevalence_72 = file.path(
      config$table_output_dir, "amr_gene_prevalence_72.csv"
    ),
    amr_gene_class_mapping_72 = file.path(
      config$table_output_dir, "amr_gene_class_mapping_72.csv"
    ),
    amr_drug_class_summary_72 = file.path(
      config$table_output_dir, "amr_drug_class_summary_72.csv"
    ),
    aminoglycoside_determinants_72 = file.path(
      config$table_output_dir, "aminoglycoside_determinants_72.csv"
    ),
    aminoglycoside_profile_by_genome_72 = file.path(
      config$table_output_dir, "aminoglycoside_profile_by_genome_72.csv"
    ),
    amr_gene_fisher_tests_72 = file.path(
      config$table_output_dir, "amr_gene_fisher_tests_72.csv"
    ),
    amr_singletons_72 = file.path(
      config$table_output_dir, "amr_singletons_72.csv"
    ),
    amr_notable_prevalence_differences_72 = file.path(
      config$table_output_dir, "amr_notable_prevalence_differences_72.csv"
    ),
    amr_pipeline_qc_summary_72 = file.path(
      config$table_output_dir, "amr_pipeline_qc_summary_72.csv"
    ),
    amr_gene_prevalence_comparison_72 = file.path(
      config$figure_output_dir, "amr_gene_prevalence_comparison_72.png"
    ),
    amr_presence_absence_heatmap_selected_72 = file.path(
      config$figure_output_dir, "amr_presence_absence_heatmap_selected_72.png"
    ),
    amr_presence_absence_heatmap_clustered_72 = file.path(
      config$figure_output_dir, "amr_presence_absence_heatmap_clustered_72.png"
    ),
    amr_presence_absence_heatmap_clustering_audit_72 = file.path(
      config$figure_output_dir,
      "amr_presence_absence_heatmap_clustering_audit_72.tsv"
    ),
    amr_drug_class_prevalence_comparison_72 = file.path(
      config$figure_output_dir, "amr_drug_class_prevalence_comparison_72.png"
    ),
    amrfinder_comparative_analysis_72_log = file.path(
      config$log_output_dir, "amrfinder_comparative_analysis_72.log"
    )
  )
}

preflight_outputs <- function(paths, overwrite) {
  duplicated_paths <- unique(paths[duplicated(normalizePath(
    paths, mustWork = FALSE
  ))])
  if (length(duplicated_paths) > 0L) {
    stop(
      "Output registry contains duplicate destination(s): ",
      paste(duplicated_paths, collapse = ", "), call. = FALSE
    )
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !overwrite) {
    stop(
      "Refusing to overwrite existing output(s):\n- ",
      paste(existing, collapse = "\n- "),
      "\nReview them and rerun with --overwrite only if replacement is intended.",
      call. = FALSE
    )
  }
}

create_output_directories <- function(config) {
  dirs <- c(
    config$processed_output_dir, config$table_output_dir,
    config$figure_output_dir, config$log_output_dir
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  failed <- dirs[!dir.exists(dirs)]
  if (length(failed) > 0L) {
    stop("Could not create output directory: ", failed[[1]], call. = FALSE)
  }
}

log_message <- function(log_path, ..., echo = TRUE) {
  text <- paste0(..., collapse = "")
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ", text)
  cat(line, "\n", file = log_path, append = TRUE, sep = "")
  if (echo) message(text)
}

initialise_log <- function(log_path, config, overwrite) {
  cat("", file = log_path, append = FALSE)
  log_message(log_path, "AMRFinderPlus comparative pipeline started")
  log_message(log_path, "Overwrite enabled: ", overwrite)
  log_message(log_path, "R version: ", R.version.string)
  log_message(log_path, "Project root: ", config$project_root)
  log_message(log_path, "Reproductive AMR directory: ", config$reproductive_amr_dir)
  log_message(log_path, "Bacteraemia AMR directory: ", config$bacteraemia_amr_dir)
  log_message(log_path, "Metadata path: ", config$metadata_path)
  for (package in required_packages) {
    log_message(
      log_path, "Package ", package, ": ",
      as.character(utils::packageVersion(package)), echo = FALSE
    )
  }
}

parse_accession_from_filename <- function(path) {
  filename <- basename(path)
  match <- stringr::str_match(
    filename, "^(GC[AF]_[0-9]+\\.[0-9]+)\\.amrfinder\\.tsv$"
  )
  if (is.na(match[1, 2])) {
    stop(
      "Unexpected TSV filename: ", filename,
      ". Expected <GCA_or_GCF_accession>.amrfinder.tsv. ",
      "Group-combined TSVs are not valid primary inputs.",
      call. = FALSE
    )
  }
  match[1, 2]
}

discover_amrfinder_files <- function(config) {
  groups <- tibble::tribble(
    ~Source, ~input_dir, ~expected_files,
    "Reproductive", config$reproductive_amr_dir, 14L,
    "Bacteraemia", config$bacteraemia_amr_dir, 58L
  )

  manifests <- purrr::pmap(groups, function(Source, input_dir, expected_files) {
    if (!dir.exists(input_dir)) {
      stop("AMRFinderPlus input directory does not exist: ", input_dir, call. = FALSE)
    }
    all_tsv <- list.files(
      input_dir, pattern = "\\.tsv$", full.names = TRUE,
      recursive = FALSE, ignore.case = TRUE
    )
    if (length(all_tsv) != expected_files) {
      stop(
        "Expected exactly ", expected_files, " per-genome TSV files in ",
        input_dir, " but found ", length(all_tsv),
        ". Remove group-combined or unrelated TSVs from the configured input directory.",
        call. = FALSE
      )
    }
    tibble::tibble(
      Genome = vapply(all_tsv, parse_accession_from_filename, character(1)),
      Source = Source,
      Original_file = basename(all_tsv),
      File_path = normalizePath(all_tsv, mustWork = TRUE)
    )
  })

  manifest <- dplyr::bind_rows(manifests)
  if (anyDuplicated(manifest$File_path)) {
    stop("A raw TSV path appears more than once in the input manifest.", call. = FALSE)
  }
  if (anyDuplicated(manifest$Genome)) {
    duplicates <- unique(manifest$Genome[duplicated(manifest$Genome)])
    stop(
      "More than one AMRFinderPlus file was found for accession(s): ",
      paste(duplicates, collapse = ", "), call. = FALSE
    )
  }
  if (nrow(manifest) != 72L) {
    stop("Expected 72 unique per-genome files but found ", nrow(manifest), ".", call. = FALSE)
  }
  manifest
}

load_metadata <- function(path) {
  if (!file.exists(path)) {
    stop("Canonical metadata file does not exist: ", path, call. = FALSE)
  }
  metadata <- readr::read_tsv(
    path, show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  required <- c("assembly_accession", "dataset_category")
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0L) {
    stop(
      "Metadata is missing required column(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  metadata <- metadata |>
    dplyr::transmute(
      Genome = stringr::str_trim(.data$assembly_accession),
      Source = stringr::str_trim(.data$dataset_category)
    )
  if (nrow(metadata) != 72L) {
    stop("Expected 72 metadata rows but found ", nrow(metadata), ".", call. = FALSE)
  }
  if (any(is.na(metadata$Genome) | metadata$Genome == "")) {
    stop("Metadata contains a missing assembly accession.", call. = FALSE)
  }
  if (anyDuplicated(metadata$Genome)) {
    duplicates <- unique(metadata$Genome[duplicated(metadata$Genome)])
    stop("Duplicated metadata accession(s): ", paste(duplicates, collapse = ", "), call. = FALSE)
  }
  observed_groups <- sort(unique(metadata$Source))
  expected_groups <- sort(names(EXPECTED_GROUP_COUNTS))
  if (!identical(observed_groups, expected_groups)) {
    stop(
      "Metadata group values must be exactly Reproductive and Bacteraemia; found: ",
      paste(observed_groups, collapse = ", "), call. = FALSE
    )
  }
  counts <- table(metadata$Source)
  for (group in names(EXPECTED_GROUP_COUNTS)) {
    if (unname(counts[[group]]) != EXPECTED_GROUP_COUNTS[[group]]) {
      stop(
        "Expected ", EXPECTED_GROUP_COUNTS[[group]], " ", group,
        " metadata rows but found ", unname(counts[[group]]), ".", call. = FALSE
      )
    }
  }
  metadata
}

load_strain_metadata <- function(path, canonical_metadata) {
  # Read strain labels from the curated metadata and require an exact 72-genome join.
  if (!file.exists(path)) {
    stop("Curated strain metadata file does not exist: ", path, call. = FALSE)
  }
  strain_metadata <- readr::read_csv(
    path, show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  required <- c(
    "assembly_accession", "strain", "reproductive_bacteraemia_category"
  )
  missing <- setdiff(required, names(strain_metadata))
  if (length(missing) > 0L) {
    stop(
      "Curated strain metadata is missing required column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  strain_metadata <- strain_metadata |>
    dplyr::transmute(
      Genome = stringr::str_trim(.data$assembly_accession),
      Strain = stringr::str_trim(.data$strain),
      Source = stringr::str_trim(.data$reproductive_bacteraemia_category)
    )
  if (nrow(strain_metadata) != 72L || anyDuplicated(strain_metadata$Genome)) {
    stop(
      "Curated strain metadata must contain exactly 72 unique accessions.",
      call. = FALSE
    )
  }
  if (any(is.na(strain_metadata$Strain) | strain_metadata$Strain == "")) {
    stop("Every canonical genome must have a non-empty strain name.", call. = FALSE)
  }
  if (anyDuplicated(strain_metadata$Strain)) {
    duplicates <- unique(strain_metadata$Strain[duplicated(strain_metadata$Strain)])
    stop(
      "Strain names must be unique for heatmap axes; duplicated value(s): ",
      paste(duplicates, collapse = ", "), call. = FALSE
    )
  }
  joined <- dplyr::full_join(
    canonical_metadata, strain_metadata,
    by = "Genome", suffix = c("_canonical", "_curated")
  )
  if (nrow(joined) != 72L || any(is.na(joined$Source_canonical)) ||
      any(is.na(joined$Source_curated))) {
    stop(
      "Curated strain metadata does not match the exact canonical accession set.",
      call. = FALSE
    )
  }
  mismatched <- joined$Source_canonical != joined$Source_curated
  if (any(mismatched)) {
    stop(
      "Curated strain metadata source disagrees with the canonical manifest for: ",
      paste(joined$Genome[mismatched], collapse = ", "), call. = FALSE
    )
  }
  joined |>
    dplyr::transmute(
      Genome = .data$Genome,
      Source = .data$Source_canonical,
      Strain = .data$Strain
    ) |>
    dplyr::arrange(.data$Genome)
}

validate_genome_membership <- function(manifest, metadata) {
  missing_files <- setdiff(metadata$Genome, manifest$Genome)
  unexpected_files <- setdiff(manifest$Genome, metadata$Genome)
  if (length(missing_files) > 0L) {
    stop(
      "Metadata accession(s) without an AMRFinderPlus file: ",
      paste(missing_files, collapse = ", "), call. = FALSE
    )
  }
  if (length(unexpected_files) > 0L) {
    stop(
      "AMRFinderPlus file accession(s) absent from metadata: ",
      paste(unexpected_files, collapse = ", "), call. = FALSE
    )
  }
  membership <- dplyr::inner_join(
    manifest, metadata, by = "Genome", suffix = c("_file", "_metadata")
  )
  mismatched <- membership |>
    dplyr::filter(.data$Source_file != .data$Source_metadata)
  if (nrow(mismatched) > 0L) {
    stop(
      "Input-directory group disagrees with metadata for accession(s): ",
      paste(mismatched$Genome, collapse = ", "), call. = FALSE
    )
  }
  manifest |>
    dplyr::select(.data$Genome, .data$Source, .data$Original_file, .data$File_path) |>
    dplyr::arrange(factor(.data$Source, levels = c("Reproductive", "Bacteraemia")), .data$Genome)
}

validate_required_columns <- function(data, path) {
  if (!identical(names(data), EXPECTED_AMRFINDER_COLUMNS)) {
    missing <- setdiff(EXPECTED_AMRFINDER_COLUMNS, names(data))
    extra <- setdiff(names(data), EXPECTED_AMRFINDER_COLUMNS)
    stop(
      "AMRFinderPlus schema mismatch in ", path,
      ". Missing: ", ifelse(length(missing), paste(missing, collapse = ", "), "none"),
      "; extra/reordered columns: ",
      ifelse(length(extra), paste(extra, collapse = ", "), "none or order differs"),
      ". Expected the exact validated 23-column schema.", call. = FALSE
    )
  }
}

read_single_amrfinder_file <- function(Genome, Source, Original_file, File_path) {
  if (file.info(File_path)$size <= 0L) {
    stop("AMRFinderPlus TSV is empty: ", File_path, call. = FALSE)
  }
  data <- readr::read_tsv(
    File_path, show_col_types = FALSE, progress = FALSE,
    name_repair = "minimal", na = c("", "NA")
  )
  validate_required_columns(data, File_path)
  if (nrow(data) == 0L) {
    stop("AMRFinderPlus TSV contains a header but no records: ", File_path, call. = FALSE)
  }
  name_values <- as.character(data$Name)
  names_in_file <- unique(stats::na.omit(name_values))
  if (any(is.na(name_values)) || any(name_values != Genome) ||
      length(names_in_file) != 1L) {
    stop(
      "Filename accession ", Genome, " does not match the unique Name value in ",
      Original_file, ". Found: ", paste(names_in_file, collapse = ", "), call. = FALSE
    )
  }
  data |>
    dplyr::mutate(
      Genome = Genome,
      Source = Source,
      Original_file = Original_file,
      Raw_row_in_file = dplyr::row_number(),
      .before = 1
    )
}

load_amrfinder_data <- function(manifest) {
  purrr::pmap_dfr(manifest, read_single_amrfinder_file)
}

standardise_amrfinder <- function(raw_data) {
  raw_data |>
    dplyr::transmute(
      Genome = .data$Genome,
      Source = .data$Source,
      Raw_record_ID = paste(.data$Genome, .data$Raw_row_in_file, sep = ":"),
      Raw_row_in_file = .data$Raw_row_in_file,
      Gene = stringr::str_trim(as.character(.data$`Element symbol`)),
      Element_name = as.character(.data$`Element name`),
      Scope = as.character(.data$Scope),
      Type = dplyr::coalesce(
        stringr::str_to_upper(stringr::str_trim(as.character(.data$Type))),
        "MISSING"
      ),
      Subtype = as.character(.data$Subtype),
      Class = as.character(.data$Class),
      Subclass = as.character(.data$Subclass),
      Method = dplyr::coalesce(
        stringr::str_to_upper(stringr::str_trim(as.character(.data$Method))),
        "MISSING"
      ),
      Protein_ID = as.character(.data$`Protein id`),
      Contig_ID = as.character(.data$`Contig id`),
      Start = as.integer(.data$Start),
      Stop = as.integer(.data$Stop),
      Strand = as.character(.data$Strand),
      Target_length = as.numeric(.data$`Target length`),
      Reference_sequence_length = as.numeric(.data$`Reference sequence length`),
      Coverage_reference_percent = as.numeric(.data$`% Coverage of reference`),
      Identity_reference_percent = as.numeric(.data$`% Identity to reference`),
      Alignment_length = as.numeric(.data$`Alignment length`),
      Closest_reference_accession = as.character(.data$`Closest reference accession`),
      Closest_reference_name = as.character(.data$`Closest reference name`),
      HMM_accession = as.character(.data$`HMM accession`),
      HMM_description = as.character(.data$`HMM description`),
      Original_file = .data$Original_file
    ) |>
    dplyr::mutate(
      Class_for_analysis = dplyr::if_else(
        is.na(.data$Class) | stringr::str_trim(.data$Class) == "",
        "other/unclassified", .data$Class
      ),
      Subclass_for_analysis = dplyr::if_else(
        is.na(.data$Subclass) | stringr::str_trim(.data$Subclass) == "",
        "other/unclassified", .data$Subclass
      )
    ) |>
    dplyr::group_by(
      .data$Genome, .data$Gene, .data$Type, .data$Method, .data$Contig_ID,
      .data$Start, .data$Stop
    ) |>
    dplyr::mutate(Potential_duplicate_record = dplyr::n() > 1L) |>
    dplyr::ungroup()
}

validate_current_data_state <- function(data) {
  if (nrow(data) != EXPECTED_RAW_ROWS) {
    stop(
      "Expected ", EXPECTED_RAW_ROWS, " raw AMRFinderPlus rows but imported ",
      nrow(data), ". The raw source may have changed; review and update the ",
      "validated expectations before analysis.", call. = FALSE
    )
  }
  if (any(is.na(data$Gene) | data$Gene == "")) {
    stop("At least one raw record has a missing Element symbol/Gene.", call. = FALSE)
  }
  observed_methods <- sort(unique(data$Method))
  if (!identical(observed_methods, sort(EXPECTED_METHODS))) {
    warning(
      "Observed Method values differ from the reviewed current dataset. Expected: ",
      paste(sort(EXPECTED_METHODS), collapse = ", "), "; found: ",
      paste(observed_methods, collapse = ", "),
      ". Unrecognised methods will be excluded from functional presence and retained ",
      "in audit outputs for manual review.", call. = FALSE
    )
  }
  type_counts <- table(data$Type)
  observed_types <- sort(names(type_counts))
  if (!identical(observed_types, sort(names(EXPECTED_TYPE_COUNTS)))) {
    stop(
      "Observed Type values differ from the reviewed AMR/STRESS policy: ",
      paste(observed_types, collapse = ", "), call. = FALSE
    )
  }
  for (type in names(EXPECTED_TYPE_COUNTS)) {
    if (unname(type_counts[[type]]) != EXPECTED_TYPE_COUNTS[[type]]) {
      stop(
        "Expected ", EXPECTED_TYPE_COUNTS[[type]], " ", type, " rows but found ",
        unname(type_counts[[type]]),
        ". The raw source may have changed and requires review.", call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

classify_method_status <- function(data) {
  data |>
    dplyr::mutate(
      Functional_status = dplyr::case_when(
        .data$Type == "AMR" & .data$Method %in% ACCEPTED_METHODS ~
          "accepted_functional_amr_detection",
        .data$Type == "AMR" & .data$Method %in% EXCLUDED_METHODS ~
          "excluded_from_functional_presence",
        .data$Type == "STRESS" ~ "outside_functional_amr_analysis",
        TRUE ~ "manual_review_required"
      ),
      Coverage_below_80 = !is.na(.data$Coverage_reference_percent) &
        .data$Coverage_reference_percent < LOW_COVERAGE_THRESHOLD,
      Flag_reason = purrr::pmap_chr(
        list(
          .data$Type, .data$Method, .data$Coverage_below_80,
          .data$Potential_duplicate_record
        ),
        function(type, method, low_coverage, potential_duplicate) {
          flags <- character()
          if (method == "PARTIALX") flags <- c(flags, "PARTIALX")
          if (method == "INTERNAL_STOP") flags <- c(flags, "INTERNAL_STOP")
          if (!method %in% EXPECTED_METHODS) {
            flags <- c(flags, paste0("unrecognised_method:", method))
          }
          if (type == "STRESS") {
            flags <- c(flags, "STRESS_outside_functional_AMR_analysis")
          } else if (type != "AMR") {
            flags <- c(flags, paste0("Type_not_AMR:", type))
          }
          if (isTRUE(low_coverage)) flags <- c(flags, "coverage_below_80")
          if (isTRUE(potential_duplicate)) flags <- c(flags, "potential_duplicate_record")
          paste(flags, collapse = ";")
        }
      ),
      Functional_matrix_status = dplyr::case_when(
        .data$Type != "AMR" ~ "outside_functional_amr_analysis",
        .data$Method %in% EXCLUDED_METHODS ~ "excluded_from_functional_presence",
        !.data$Method %in% ACCEPTED_METHODS ~ "excluded_pending_method_review",
        .data$Coverage_below_80 ~ "included_in_functional_presence_but_coverage_flagged",
        TRUE ~ "included_in_functional_presence"
      )
    )
}

build_method_audit <- function(data) {
  data |>
    dplyr::group_by(.data$Method) |>
    dplyr::summarise(
      n_records = dplyr::n(),
      min_coverage = min(.data$Coverage_reference_percent, na.rm = TRUE),
      max_coverage = max(.data$Coverage_reference_percent, na.rm = TRUE),
      min_identity = min(.data$Identity_reference_percent, na.rm = TRUE),
      max_identity = max(.data$Identity_reference_percent, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      functional_status = dplyr::case_when(
        .data$Method %in% ACCEPTED_METHODS ~ "accepted_when_Type_is_AMR",
        .data$Method %in% EXCLUDED_METHODS ~ "excluded_from_functional_presence",
        TRUE ~ "manual_review_required"
      ),
      notes = dplyr::case_when(
        .data$Method == "EXACTX" ~ "Complete exact protein match",
        .data$Method == "BLASTX" ~ "High-coverage translated sequence match in current data",
        .data$Method == "POINTX" ~ "Reviewed point-mutation detection with complete current-data coverage",
        .data$Method == "PARTIALX" ~ "Partial hit retained in flagged audit output",
        .data$Method == "INTERNAL_STOP" ~ "Interrupted coding sequence retained in flagged audit output",
        TRUE ~ "Method requires manual review"
      )
    ) |>
    dplyr::arrange(match(.data$Method, EXPECTED_METHODS))
}

build_flagged_hits <- function(data) {
  data |>
    dplyr::filter(
      .data$Functional_status != "accepted_functional_amr_detection" |
        .data$Flag_reason != ""
    ) |>
    dplyr::select(
      .data$Raw_record_ID, .data$Raw_row_in_file,
      .data$Genome, .data$Source, .data$Gene, .data$Element_name,
      .data$Type, .data$Method,
      .data$Class, .data$Subclass, .data$Contig_ID,
      .data$Start, .data$Stop, .data$Strand,
      .data$Coverage_reference_percent, .data$Identity_reference_percent,
      .data$Target_length, .data$Reference_sequence_length, .data$Alignment_length,
      .data$Potential_duplicate_record,
      flag_reason = .data$Flag_reason,
      functional_matrix_status = .data$Functional_matrix_status,
      .data$Original_file
    ) |>
    dplyr::arrange(.data$Source, .data$Genome, .data$Gene, .data$Method)
}

accepted_functional_hits <- function(data) {
  # Apply the reviewed AMR Type and Method policy without a coverage exclusion.
  data |>
    dplyr::filter(.data$Type == "AMR", .data$Method %in% ACCEPTED_METHODS)
}

build_presence_absence <- function(functional_hits, metadata) {
  # Expand every metadata genome across accepted genes so zero-hit genomes remain.
  genes <- sort(unique(functional_hits$Gene))
  if (length(genes) == 0L) {
    stop("No accepted functional AMR genes remain after filtering.", call. = FALSE)
  }
  presence_long <- functional_hits |>
    dplyr::distinct(.data$Genome, .data$Gene) |>
    dplyr::mutate(Present = 1L)
  matrix <- tidyr::crossing(metadata, Gene = genes) |>
    dplyr::left_join(presence_long, by = c("Genome", "Gene")) |>
    dplyr::mutate(Present = dplyr::coalesce(.data$Present, 0L)) |>
    tidyr::pivot_wider(
      names_from = .data$Gene, values_from = .data$Present,
      values_fill = 0L, values_fn = max
    ) |>
    dplyr::arrange(
      factor(.data$Source, levels = c("Reproductive", "Bacteraemia")),
      .data$Genome
    )
  if (nrow(matrix) != 72L || anyDuplicated(matrix$Genome)) {
    stop("Presence/absence matrix must contain one row for each of 72 genomes.", call. = FALSE)
  }
  if (any(is.na(matrix$Source))) {
    stop("Presence/absence matrix contains a missing Source value.", call. = FALSE)
  }
  matrix_group_counts <- table(matrix$Source)
  observed_group_counts <- as.integer(
    matrix_group_counts[c("Reproductive", "Bacteraemia")]
  )
  required_group_counts <- as.integer(
    EXPECTED_GROUP_COUNTS[c("Reproductive", "Bacteraemia")]
  )
  if (!identical(observed_group_counts, required_group_counts)) {
    stop(
      "Presence/absence matrix must contain exactly 14 Reproductive and 58 ",
      "Bacteraemia genomes.", call. = FALSE
    )
  }
  matrix_values <- unlist(
    matrix[setdiff(names(matrix), c("Genome", "Source"))],
    use.names = FALSE
  )
  if (any(!matrix_values %in% c(0L, 1L))) {
    stop("Presence/absence matrix contains a value other than 0 or 1.", call. = FALSE)
  }
  matrix
}

presence_matrix_long <- function(matrix) {
  matrix |>
    tidyr::pivot_longer(
      cols = -c("Genome", "Source"),
      names_to = "Gene", values_to = "Present"
    )
}

summarise_gene_prevalence <- function(matrix) {
  presence_matrix_long(matrix) |>
    dplyr::group_by(.data$Gene) |>
    dplyr::summarise(
      Reproductive_count = sum(.data$Present[.data$Source == "Reproductive"]),
      Reproductive_total = sum(.data$Source == "Reproductive"),
      Bacteraemia_count = sum(.data$Present[.data$Source == "Bacteraemia"]),
      Bacteraemia_total = sum(.data$Source == "Bacteraemia"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Reproductive_percent = 100 * .data$Reproductive_count / .data$Reproductive_total,
      Bacteraemia_percent = 100 * .data$Bacteraemia_count / .data$Bacteraemia_total,
      Combined_count = .data$Reproductive_count + .data$Bacteraemia_count,
      Combined_total = .data$Reproductive_total + .data$Bacteraemia_total,
      Combined_percent = 100 * .data$Combined_count / .data$Combined_total,
      Prevalence_difference_percentage_points =
        .data$Reproductive_percent - .data$Bacteraemia_percent
    ) |>
    dplyr::arrange(dplyr::desc(.data$Combined_percent), .data$Gene)
}

format_prevalence_table <- function(prevalence) {
  prevalence |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(c(
        "Reproductive_percent", "Bacteraemia_percent", "Combined_percent",
        "Prevalence_difference_percentage_points"
      )),
      ~ round(.x, 2)
    ))
}

classify_shared_dataset_only <- function(prevalence) {
  prevalence |>
    dplyr::transmute(
      Gene = .data$Gene,
      Classification = dplyr::case_when(
        .data$Reproductive_count >= 1L & .data$Bacteraemia_count >= 1L ~ "Shared",
        .data$Reproductive_count >= 1L & .data$Bacteraemia_count == 0L ~
          "detected only in the reproductive-associated dataset (n=14)",
        .data$Reproductive_count == 0L & .data$Bacteraemia_count >= 1L ~
          "detected only in the bacteraemia-associated dataset (n=58)",
        TRUE ~ "not detected in either dataset"
      ),
      Reproductive_count = .data$Reproductive_count,
      Bacteraemia_count = .data$Bacteraemia_count,
      Total_count = .data$Combined_count
    ) |>
    dplyr::arrange(dplyr::desc(.data$Total_count), .data$Gene)
}

build_gene_class_mapping <- function(data) {
  mappings <- data |>
    dplyr::filter(.data$Type == "AMR") |>
    dplyr::count(
      .data$Gene, Class = .data$Class_for_analysis,
      Subclass = .data$Subclass_for_analysis,
      name = "n_records"
    )
  status <- mappings |>
    dplyr::group_by(.data$Gene) |>
    dplyr::summarise(n_mappings = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(mapping_status = dplyr::if_else(
      .data$n_mappings == 1L, "consistent", "review_required"
    ))
  mappings |>
    dplyr::left_join(status, by = "Gene") |>
    dplyr::select(
      .data$Gene, .data$Class, .data$Subclass, .data$n_records,
      .data$mapping_status
    ) |>
    dplyr::arrange(.data$Gene, .data$Class, .data$Subclass)
}

classify_broad_classes_one <- function(class_value, subclass_value) {
  class_missing <- is.na(class_value) || stringr::str_trim(class_value) == "" ||
    class_value == "other/unclassified"
  subclass_missing <- is.na(subclass_value) || stringr::str_trim(subclass_value) == "" ||
    subclass_value == "other/unclassified"
  if (class_missing && subclass_missing) return("other/unclassified")

  text <- stringr::str_to_lower(paste(class_value, subclass_value))
  patterns <- c(
    aminoglycoside = "aminoglycoside|gentamicin|amikacin|kanamycin|tobramycin|streptomycin",
    tetracycline = "tetracycline|tigecycline",
    macrolide = "macrolide|erythromycin",
    lincosamide = "lincosamide|lincomycin|clindamycin",
    streptogramin = "streptogramin",
    glycopeptide = "glycopeptide|vancomycin|teicoplanin"
  )
  matched <- names(patterns)[vapply(
    patterns, function(pattern) stringr::str_detect(text, pattern), logical(1)
  )]
  if (length(matched) == 0L) "other" else matched
}

add_broad_drug_classes <- function(data) {
  data |>
    dplyr::mutate(Broad_drug_class = purrr::map2(
      .data$Class_for_analysis, .data$Subclass_for_analysis,
      classify_broad_classes_one
    )) |>
    tidyr::unnest_longer("Broad_drug_class")
}

collapse_values <- function(values) {
  values <- sort(unique(stats::na.omit(as.character(values))))
  values <- values[values != ""]
  if (length(values) == 0L) "" else paste(values, collapse = "; ")
}

summarise_drug_classes <- function(functional_hits, metadata) {
  broad <- add_broad_drug_classes(functional_hits)
  required_classes <- c(
    "aminoglycoside", "tetracycline", "macrolide", "lincosamide",
    "streptogramin", "glycopeptide", "other", "other/unclassified"
  )
  genes <- broad |>
    dplyr::group_by(.data$Broad_drug_class) |>
    dplyr::summarise(genes_detected_in_class = collapse_values(.data$Gene), .groups = "drop")
  genome_presence <- broad |>
    dplyr::distinct(.data$Genome, .data$Broad_drug_class) |>
    dplyr::left_join(metadata, by = "Genome") |>
    dplyr::group_by(.data$Broad_drug_class) |>
    dplyr::summarise(
      reproductive_genomes_with_any_gene_in_class = sum(.data$Source == "Reproductive"),
      bacteraemia_genomes_with_any_gene_in_class = sum(.data$Source == "Bacteraemia"),
      combined_genomes_with_any_gene_in_class = dplyr::n_distinct(.data$Genome),
      .groups = "drop"
    )
  tibble::tibble(Broad_drug_class = required_classes) |>
    dplyr::left_join(genes, by = "Broad_drug_class") |>
    dplyr::left_join(genome_presence, by = "Broad_drug_class") |>
    dplyr::mutate(
      genes_detected_in_class = dplyr::coalesce(.data$genes_detected_in_class, ""),
      dplyr::across(
        dplyr::all_of(c(
          "reproductive_genomes_with_any_gene_in_class",
          "bacteraemia_genomes_with_any_gene_in_class",
          "combined_genomes_with_any_gene_in_class"
        )),
        ~ dplyr::coalesce(as.integer(.x), 0L)
      ),
      reproductive_percent_with_any_gene_in_class = round(
        100 * .data$reproductive_genomes_with_any_gene_in_class / 14, 2
      ),
      bacteraemia_percent_with_any_gene_in_class = round(
        100 * .data$bacteraemia_genomes_with_any_gene_in_class / 58, 2
      ),
      combined_percent_with_any_gene_in_class = round(
        100 * .data$combined_genomes_with_any_gene_in_class / 72, 2
      )
    )
}

build_aminoglycoside_outputs <- function(data, metadata) {
  aminoglycoside <- data |>
    dplyr::filter(.data$Type == "AMR") |>
    add_broad_drug_classes() |>
    dplyr::filter(.data$Broad_drug_class == "aminoglycoside")

  determinants <- aminoglycoside |>
    dplyr::transmute(
      .data$Raw_record_ID, .data$Raw_row_in_file,
      .data$Genome, .data$Source, .data$Gene, .data$Method,
      Class = .data$Class_for_analysis,
      Subclass = .data$Subclass_for_analysis,
      .data$Contig_ID, .data$Coverage_reference_percent,
      .data$Identity_reference_percent,
      functional_status = .data$Functional_status,
      flag_reason = .data$Flag_reason,
      .data$Original_file
    ) |>
    dplyr::arrange(.data$Source, .data$Genome, .data$Gene)

  accepted <- aminoglycoside |>
    dplyr::filter(.data$Functional_status == "accepted_functional_amr_detection") |>
    dplyr::group_by(.data$Genome) |>
    dplyr::summarise(
      accepted_aminoglycoside_genes = collapse_values(.data$Gene),
      aminoglycoside_gene_count = dplyr::n_distinct(.data$Gene),
      .groups = "drop"
    )
  flagged <- aminoglycoside |>
    dplyr::filter(.data$Flag_reason != "") |>
    dplyr::mutate(description = paste(.data$Gene, .data$Method, .data$Flag_reason, sep = "|")) |>
    dplyr::group_by(.data$Genome) |>
    dplyr::summarise(
      flagged_aminoglycoside_hits = collapse_values(.data$description),
      .groups = "drop"
    )
  profile <- metadata |>
    dplyr::left_join(accepted, by = "Genome") |>
    dplyr::left_join(flagged, by = "Genome") |>
    dplyr::mutate(
      accepted_aminoglycoside_genes = dplyr::coalesce(
        .data$accepted_aminoglycoside_genes, ""
      ),
      aminoglycoside_gene_count = dplyr::coalesce(
        .data$aminoglycoside_gene_count, 0L
      ),
      flagged_aminoglycoside_hits = dplyr::coalesce(
        .data$flagged_aminoglycoside_hits, ""
      )
    ) |>
    dplyr::arrange(
      factor(.data$Source, levels = c("Reproductive", "Bacteraemia")),
      .data$Genome
    )
  list(determinants = determinants, profile = profile)
}

run_gene_fisher_tests <- function(prevalence) {
  # Test matrix-derived gene counts only when at least two genomes carry the gene.
  tested <- prevalence |>
    dplyr::filter(.data$Combined_count >= 2L)
  results <- purrr::pmap_dfr(
    tested,
    function(Gene, Reproductive_count, Reproductive_total,
             Bacteraemia_count, Bacteraemia_total, ...) {
      contingency <- matrix(
        c(
          Reproductive_count, Reproductive_total - Reproductive_count,
          Bacteraemia_count, Bacteraemia_total - Bacteraemia_count
        ),
        nrow = 2L, byrow = TRUE,
        dimnames = list(
          Source = c("Reproductive", "Bacteraemia"),
          Detection = c("Present", "Absent")
        )
      )
      test <- stats::fisher.test(contingency, alternative = "two.sided")
      tibble::tibble(
        Gene = Gene,
        Odds_ratio = if (length(test$estimate) == 0L) NA_real_ else unname(test$estimate),
        Raw_p_value = test$p.value
      )
    }
  )
  tested |>
    dplyr::left_join(results, by = "Gene") |>
    dplyr::mutate(
      BH_adjusted_p_value = stats::p.adjust(.data$Raw_p_value, method = "BH"),
      Significant_q_lt_0_05 = .data$BH_adjusted_p_value < 0.05,
      Reproductive_percent = round(.data$Reproductive_percent, 2),
      Bacteraemia_percent = round(.data$Bacteraemia_percent, 2),
      Prevalence_difference_percentage_points = round(
        .data$Prevalence_difference_percentage_points, 2
      )
    ) |>
    dplyr::select(
      .data$Gene, .data$Reproductive_count, .data$Reproductive_total,
      .data$Reproductive_percent, .data$Bacteraemia_count,
      .data$Bacteraemia_total, .data$Bacteraemia_percent,
      .data$Prevalence_difference_percentage_points, .data$Odds_ratio,
      .data$Raw_p_value, .data$BH_adjusted_p_value,
      .data$Significant_q_lt_0_05
    ) |>
    dplyr::arrange(.data$Raw_p_value, .data$Gene)
}

build_singleton_table <- function(prevalence, functional_hits) {
  singleton_genes <- prevalence |>
    dplyr::filter(.data$Combined_count == 1L) |>
    dplyr::pull(.data$Gene)
  functional_hits |>
    dplyr::filter(.data$Gene %in% singleton_genes) |>
    dplyr::transmute(
      .data$Gene, .data$Genome, .data$Source, .data$Method,
      Class = .data$Class_for_analysis,
      Subclass = .data$Subclass_for_analysis,
      .data$Contig_ID,
      statistical_test_status =
        "insufficient occurrences for gene-level statistical testing"
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$Gene, .data$Genome, .data$Contig_ID)
}

build_notable_difference_table <- function(prevalence, fisher_results) {
  prevalence |>
    dplyr::filter(abs(.data$Prevalence_difference_percentage_points) >= 20) |>
    dplyr::select(
      .data$Gene, .data$Reproductive_percent, .data$Bacteraemia_percent,
      .data$Prevalence_difference_percentage_points
    ) |>
    dplyr::left_join(
      fisher_results |>
        dplyr::select(
          .data$Gene, .data$Odds_ratio, .data$Raw_p_value,
          .data$BH_adjusted_p_value
        ),
      by = "Gene"
    ) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(c(
        "Reproductive_percent", "Bacteraemia_percent",
        "Prevalence_difference_percentage_points"
      )),
      ~ round(.x, 2)
    )) |>
    dplyr::arrange(
      dplyr::desc(abs(.data$Prevalence_difference_percentage_points)),
      .data$Gene
    )
}

build_functional_hits_with_contigs <- function(functional_hits) {
  functional_hits |>
    dplyr::transmute(
      .data$Raw_record_ID, .data$Raw_row_in_file,
      .data$Genome, .data$Source, .data$Gene, .data$Method,
      Class = .data$Class_for_analysis,
      Subclass = .data$Subclass_for_analysis,
      .data$Contig_ID, .data$Coverage_reference_percent,
      .data$Identity_reference_percent, .data$Start, .data$Stop, .data$Strand
    ) |>
    dplyr::arrange(.data$Source, .data$Genome, .data$Contig_ID, .data$Start)
}

build_qc_summary <- function(data, matrix, prevalence, fisher_results, classifications) {
  method_counts <- table(factor(data$Method, levels = EXPECTED_METHODS))
  type_counts <- table(factor(data$Type, levels = names(EXPECTED_TYPE_COUNTS)))
  class_counts <- table(classifications$Classification)
  value_or_zero <- function(table_object, key) {
    if (key %in% names(table_object)) as.integer(table_object[[key]]) else 0L
  }
  shared_label <- "Shared"
  reproductive_label <- "detected only in the reproductive-associated dataset (n=14)"
  bacteraemia_label <- "detected only in the bacteraemia-associated dataset (n=58)"

  tibble::tibble(
    metric = c(
      "expected_genomes", "observed_genomes", "reproductive_genomes",
      "bacteraemia_genomes", "raw_rows", "amr_rows", "stress_rows",
      "accepted_functional_amr_rows", "exactx_rows", "blastx_rows",
      "pointx_rows", "partialx_rows", "internal_stop_rows",
      "unrecognised_method_rows", "coverage_below_80_rows",
      "potential_duplicate_rows", "genes_in_functional_matrix", "genes_tested",
      "singleton_genes", "shared_genes", "reproductive_dataset_only_genes",
      "bacteraemia_dataset_only_genes"
    ),
    value = c(
      72L, dplyr::n_distinct(matrix$Genome),
      sum(matrix$Source == "Reproductive"), sum(matrix$Source == "Bacteraemia"),
      nrow(data), as.integer(type_counts[["AMR"]]),
      as.integer(type_counts[["STRESS"]]),
      sum(data$Functional_status == "accepted_functional_amr_detection"),
      as.integer(method_counts[["EXACTX"]]), as.integer(method_counts[["BLASTX"]]),
      as.integer(method_counts[["POINTX"]]), as.integer(method_counts[["PARTIALX"]]),
      as.integer(method_counts[["INTERNAL_STOP"]]),
      sum(!data$Method %in% EXPECTED_METHODS), sum(data$Coverage_below_80),
      sum(data$Potential_duplicate_record),
      nrow(prevalence), nrow(fisher_results), sum(prevalence$Combined_count == 1L),
      value_or_zero(class_counts, shared_label),
      value_or_zero(class_counts, reproductive_label),
      value_or_zero(class_counts, bacteraemia_label)
    )
  )
}

plot_gene_prevalence <- function(prevalence) {
  plot_data <- prevalence |>
    dplyr::select(
      .data$Gene, Reproductive = .data$Reproductive_percent,
      Bacteraemia = .data$Bacteraemia_percent,
      .data$Combined_percent
    ) |>
    dplyr::slice_max(.data$Combined_percent, n = HEATMAP_TOP_N_GENES, with_ties = FALSE) |>
    tidyr::pivot_longer(
      cols = c("Reproductive", "Bacteraemia"),
      names_to = "Source", values_to = "Prevalence_percent"
    ) |>
    dplyr::mutate(Gene = stats::reorder(.data$Gene, .data$Combined_percent))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Gene, y = .data$Prevalence_percent, fill = .data$Source)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Prevalence of accepted functional AMR determinants",
      subtitle = "Top genes by combined prevalence",
      x = "AMR determinant", y = "Genomes with determinant (%)", fill = "Source"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

selected_heatmap_genes <- function(prevalence, functional_hits) {
  # Select the documented top-prevalence genes plus every aminoglycoside gene.
  top_genes <- prevalence |>
    dplyr::slice_max(.data$Combined_percent, n = HEATMAP_TOP_N_GENES, with_ties = FALSE) |>
    dplyr::pull(.data$Gene)
  functional_for_classes <- functional_hits
  if (!"Class_for_analysis" %in% names(functional_for_classes) &&
      "Class" %in% names(functional_for_classes)) {
    functional_for_classes$Class_for_analysis <- functional_for_classes$Class
  }
  if (!"Subclass_for_analysis" %in% names(functional_for_classes) &&
      "Subclass" %in% names(functional_for_classes)) {
    functional_for_classes$Subclass_for_analysis <- functional_for_classes$Subclass
  }
  required <- c("Class_for_analysis", "Subclass_for_analysis")
  missing <- setdiff(required, names(functional_for_classes))
  if (length(missing) > 0L) {
    stop(
      "Functional-hit table is missing heatmap class column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  aminoglycoside_genes <- functional_for_classes |>
    add_broad_drug_classes() |>
    dplyr::filter(.data$Broad_drug_class == "aminoglycoside") |>
    dplyr::pull(.data$Gene) |>
    unique()
  unique(c(top_genes, sort(aminoglycoside_genes)))
}

heatmap_detection_scale <- function() {
  # Use one explicit blue/red binary palette for every AMR heatmap variant.
  ggplot2::scale_fill_manual(
    values = c(Absent = "#2166AC", Present = "#B2182B"),
    drop = FALSE
  )
}

binary_jaccard_matrix <- function(data_matrix) {
  # Calculate Jaccard distances on raw binary calls, including a defined zero-union case.
  binary <- data_matrix != 0
  n_items <- nrow(binary)
  result <- matrix(
    0, nrow = n_items, ncol = n_items,
    dimnames = list(rownames(binary), rownames(binary))
  )
  if (n_items < 2L) return(result)
  for (i in seq_len(n_items - 1L)) {
    for (j in seq.int(i + 1L, n_items)) {
      union_count <- sum(binary[i, ] | binary[j, ])
      distance <- if (union_count == 0L) {
        0
      } else {
        sum(xor(binary[i, ], binary[j, ])) / union_count
      }
      result[i, j] <- distance
      result[j, i] <- distance
    }
  }
  result
}

distance_matrix_audit <- function(
    distance_matrix, dimension, accession_lookup = NULL, cluster_order) {
  # Store the complete symmetric matrix in long form with reproducible clustering metadata.
  audit <- as.data.frame(as.table(distance_matrix), stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    dplyr::rename(
      Item_1_key = Var1,
      Item_2_key = Var2,
      Distance = Freq
    ) |>
    dplyr::mutate(
      Heatmap_variant = "clustered",
      Dimension = dimension,
      Item_1 = .data$Item_1_key,
      Item_2 = .data$Item_2_key,
      Item_1_accession = NA_character_,
      Item_2_accession = NA_character_,
      Item_1_cluster_order = match(.data$Item_1_key, cluster_order),
      Item_2_cluster_order = match(.data$Item_2_key, cluster_order),
      Distance_method = "binary Jaccard; joint absences ignored; zero union defined as distance 0",
      Scaling_normalisation = "none; raw accepted-detection presence/absence values (0/1)",
      Linkage_method = "average (UPGMA)"
    )
  if (!is.null(accession_lookup)) {
    audit <- audit |>
      dplyr::mutate(
        Item_1 = unname(accession_lookup[.data$Item_1_key]),
        Item_2 = unname(accession_lookup[.data$Item_2_key]),
        Item_1_accession = .data$Item_1_key,
        Item_2_accession = .data$Item_2_key
      )
  }
  audit |>
    dplyr::select(dplyr::all_of(c(
      "Heatmap_variant", "Dimension", "Item_1", "Item_2",
      "Item_1_accession", "Item_2_accession", "Distance",
      "Item_1_cluster_order", "Item_2_cluster_order", "Distance_method",
      "Scaling_normalisation", "Linkage_method"
    )))
}

plot_presence_absence_heatmap <- function(
    matrix, prevalence, functional_hits, strain_metadata) {
  # Keep the original deterministic ordering while replacing accessions with strain labels.
  selected_genes <- selected_heatmap_genes(prevalence, functional_hits)
  gene_order <- prevalence |>
    dplyr::filter(.data$Gene %in% selected_genes) |>
    dplyr::arrange(.data$Combined_count, .data$Gene) |>
    dplyr::pull(.data$Gene)
  genome_order <- matrix |>
    dplyr::arrange(
      factor(.data$Source, levels = c("Bacteraemia", "Reproductive")),
      .data$Genome
    ) |>
    dplyr::pull(.data$Genome)
  strain_order <- strain_metadata$Strain[match(genome_order, strain_metadata$Genome)]
  plot_data <- presence_matrix_long(matrix) |>
    dplyr::filter(.data$Gene %in% selected_genes) |>
    dplyr::left_join(strain_metadata, by = c("Genome", "Source")) |>
    dplyr::mutate(
      Gene = factor(.data$Gene, levels = gene_order),
      Strain = factor(.data$Strain, levels = strain_order),
      Present = factor(.data$Present, levels = c(0, 1), labels = c("Absent", "Present"))
    )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Gene, y = .data$Strain, fill = .data$Present)
  ) +
    ggplot2::geom_tile() +
    heatmap_detection_scale() +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Source), scales = "free_y", space = "free_y"
    ) +
    ggplot2::labs(
      title = "Accepted functional AMR determinant presence/absence",
      subtitle = "Top 30 genes by prevalence plus all accepted aminoglycoside genes",
      x = "AMR determinant", y = "Strain", fill = "Detection"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

build_clustered_presence_absence_heatmap <- function(
    matrix, prevalence, functional_hits, strain_metadata) {
  # Cluster both genomes and genes using average-linkage Jaccard distances.
  selected_genes <- selected_heatmap_genes(prevalence, functional_hits)
  binary <- matrix |>
    dplyr::arrange(.data$Genome) |>
    dplyr::select("Genome", dplyr::all_of(selected_genes))
  numeric_matrix <- as.matrix(binary[, selected_genes, drop = FALSE])
  storage.mode(numeric_matrix) <- "numeric"
  rownames(numeric_matrix) <- binary$Genome

  row_distance <- binary_jaccard_matrix(numeric_matrix)
  column_distance <- binary_jaccard_matrix(t(numeric_matrix))
  row_cluster <- stats::hclust(stats::as.dist(row_distance), method = "average")
  column_cluster <- stats::hclust(
    stats::as.dist(column_distance), method = "average"
  )
  genome_order <- rownames(numeric_matrix)[row_cluster$order]
  gene_order <- colnames(numeric_matrix)[column_cluster$order]
  strain_lookup <- stats::setNames(
    strain_metadata$Strain, strain_metadata$Genome
  )
  strain_order <- unname(strain_lookup[genome_order])

  plot_data <- presence_matrix_long(matrix) |>
    dplyr::filter(.data$Gene %in% selected_genes) |>
    dplyr::left_join(strain_metadata, by = c("Genome", "Source")) |>
    dplyr::mutate(
      Gene = factor(.data$Gene, levels = gene_order),
      Strain = factor(.data$Strain, levels = rev(strain_order)),
      Present = factor(
        .data$Present, levels = c(0, 1), labels = c("Absent", "Present")
      )
    )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Gene, y = .data$Strain, fill = .data$Present)
  ) +
    ggplot2::geom_tile() +
    heatmap_detection_scale() +
    ggplot2::labs(
      title = "Clustered accepted functional AMR determinant presence/absence",
      subtitle = paste(
        "Rows and columns: binary Jaccard distance, average linkage;",
        "raw 0/1 calls without scaling"
      ),
      x = "AMR determinant", y = "Strain", fill = "Detection"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  unclustered_policy <- tibble::tibble(
    Heatmap_variant = "unclustered",
    Dimension = "rows_and_columns",
    Item_1 = NA_character_, Item_2 = NA_character_,
    Item_1_accession = NA_character_, Item_2_accession = NA_character_,
    Distance = NA_real_,
    Item_1_cluster_order = NA_integer_, Item_2_cluster_order = NA_integer_,
    Distance_method = "not applied; rows ordered by source/accession and columns by prevalence/gene",
    Scaling_normalisation = "none; raw accepted-detection presence/absence values (0/1)",
    Linkage_method = "not applied"
  )
  row_audit <- distance_matrix_audit(
    row_distance, "rows (strains)", strain_lookup, genome_order
  )
  column_audit <- distance_matrix_audit(
    column_distance, "columns (AMR determinants)", NULL, gene_order
  )
  list(
    plot = plot,
    audit = dplyr::bind_rows(unclustered_policy, row_audit, column_audit)
  )
}

plot_drug_class_prevalence <- function(drug_class_summary) {
  drug_class_summary |>
    dplyr::select(
      .data$Broad_drug_class,
      Reproductive = .data$reproductive_percent_with_any_gene_in_class,
      Bacteraemia = .data$bacteraemia_percent_with_any_gene_in_class
    ) |>
    tidyr::pivot_longer(
      cols = c("Reproductive", "Bacteraemia"),
      names_to = "Source", values_to = "Prevalence_percent"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$Broad_drug_class, y = .data$Prevalence_percent,
        fill = .data$Source
      )
    ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "Genome-level prevalence of accepted AMR drug classes",
      x = "Broad drug class", y = "Genomes with at least one determinant (%)",
      fill = "Source"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

safe_write_csv <- function(data, path, overwrite) {
  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  readr::write_csv(data, path, na = "")
}

safe_write_tsv <- function(data, path, overwrite) {
  # Apply the same overwrite guard to the clustering-distance audit TSV.
  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  readr::write_tsv(data, path, na = "")
}

safe_save_plot <- function(plot, path, overwrite, width, height) {
  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing figure: ", path, call. = FALSE)
  }
  ggplot2::ggsave(
    filename = path, plot = plot, width = width, height = height,
    dpi = 300, bg = "white"
  )
}

write_outputs <- function(objects, plots, paths, overwrite) {
  csv_names <- setdiff(names(paths), c(
    "amr_gene_prevalence_comparison_72",
    "amr_presence_absence_heatmap_selected_72",
    "amr_presence_absence_heatmap_clustered_72",
    "amr_presence_absence_heatmap_clustering_audit_72",
    "amr_drug_class_prevalence_comparison_72",
    "amrfinder_comparative_analysis_72_log"
  ))
  missing_objects <- setdiff(csv_names, names(objects))
  if (length(missing_objects) > 0L) {
    stop(
      "No data object was registered for output(s): ",
      paste(missing_objects, collapse = ", "), call. = FALSE
    )
  }
  for (name in csv_names) {
    safe_write_csv(objects[[name]], paths[[name]], overwrite)
  }
  safe_write_tsv(
    objects$amr_presence_absence_heatmap_clustering_audit_72,
    paths[["amr_presence_absence_heatmap_clustering_audit_72"]], overwrite
  )
  safe_save_plot(
    plots$gene_prevalence, paths[["amr_gene_prevalence_comparison_72"]],
    overwrite, width = 10, height = 9
  )
  safe_save_plot(
    plots$presence_absence,
    paths[["amr_presence_absence_heatmap_selected_72"]],
    overwrite, width = 14, height = 16
  )
  safe_save_plot(
    plots$presence_absence_clustered,
    paths[["amr_presence_absence_heatmap_clustered_72"]],
    overwrite, width = 14, height = 16
  )
  safe_save_plot(
    plots$drug_class, paths[["amr_drug_class_prevalence_comparison_72"]],
    overwrite, width = 10, height = 7
  )
}

run_heatmaps_only <- function(args, config, paths) {
  # Regenerate only heatmap artifacts from frozen validated AMR tables.
  heatmap_paths <- paths[c(
    "amr_presence_absence_heatmap_selected_72",
    "amr_presence_absence_heatmap_clustered_72",
    "amr_presence_absence_heatmap_clustering_audit_72"
  )]
  preflight_outputs(heatmap_paths, args$overwrite)
  dir.create(config$figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  required_inputs <- c(
    matrix = paths[["amr_gene_presence_absence"]],
    prevalence = paths[["amr_gene_prevalence_72"]],
    functional_hits = file.path(
      config$processed_output_dir, "amr_functional_hits_with_contigs_72.csv"
    )
  )
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs) > 0L) {
    stop(
      "Heatmaps-only mode is missing validated upstream input(s): ",
      paste(missing_inputs, collapse = ", "), call. = FALSE
    )
  }
  matrix <- readr::read_csv(
    required_inputs[["matrix"]], show_col_types = FALSE, progress = FALSE
  )
  prevalence <- readr::read_csv(
    required_inputs[["prevalence"]], show_col_types = FALSE, progress = FALSE
  )
  functional <- readr::read_csv(
    required_inputs[["functional_hits"]], show_col_types = FALSE, progress = FALSE
  )
  metadata <- load_metadata(config$metadata_path)
  strain_metadata <- load_strain_metadata(config$strain_metadata_path, metadata)
  if (nrow(matrix) != 72L || anyDuplicated(matrix$Genome) ||
      !setequal(matrix$Genome, metadata$Genome)) {
    stop(
      "Frozen AMR matrix must contain exactly the 72 canonical genomes.",
      call. = FALSE
    )
  }
  clustered <- build_clustered_presence_absence_heatmap(
    matrix, prevalence, functional, strain_metadata
  )
  unclustered <- plot_presence_absence_heatmap(
    matrix, prevalence, functional, strain_metadata
  )
  safe_save_plot(
    unclustered, heatmap_paths[["amr_presence_absence_heatmap_selected_72"]],
    args$overwrite, width = 14, height = 16
  )
  safe_save_plot(
    clustered$plot,
    heatmap_paths[["amr_presence_absence_heatmap_clustered_72"]],
    args$overwrite, width = 14, height = 16
  )
  safe_write_tsv(
    clustered$audit,
    heatmap_paths[["amr_presence_absence_heatmap_clustering_audit_72"]],
    args$overwrite
  )
  message("Heatmap-only regeneration completed successfully.")
}

main <- function() {
  # Parse options and preflight every destination before creating runtime files.
  args <- parse_args()
  check_packages(required_packages)
  config <- get_config()
  paths <- output_registry(config)
  if (args$heatmaps_only) {
    run_heatmaps_only(args, config, paths)
    return(invisible(NULL))
  }
  preflight_outputs(paths, args$overwrite)
  create_output_directories(config)
  log_path <- paths[["amrfinder_comparative_analysis_72_log"]]
  initialise_log(log_path, config, args$overwrite)

  # Discover per-genome TSVs and validate their one-to-one metadata membership.
  manifest <- discover_amrfinder_files(config)
  log_message(log_path, "Per-genome files discovered: ", nrow(manifest))
  metadata <- load_metadata(config$metadata_path)
  strain_metadata <- load_strain_metadata(
    config$strain_metadata_path, metadata
  )
  manifest <- validate_genome_membership(manifest, metadata)
  log_message(log_path, "Metadata and file membership validation passed")

  # Import all raw rows, standardise the exact schema and attach audit statuses.
  raw <- load_amrfinder_data(manifest)
  data <- standardise_amrfinder(raw)
  validate_current_data_state(data)
  data <- classify_method_status(data)
  log_message(log_path, "Raw rows imported and validated: ", nrow(data))
  log_message(
    log_path, "Observed Method values: ",
    paste(sort(unique(data$Method)), collapse = ", ")
  )

  # Build QC/audit objects before constructing the accepted functional matrix.
  method_audit <- build_method_audit(data)
  flagged <- build_flagged_hits(data)
  functional <- accepted_functional_hits(data)
  matrix <- build_presence_absence(functional, metadata)

  # Derive prevalence, class, aminoglycoside and statistical summaries from the matrix policy.
  prevalence <- summarise_gene_prevalence(matrix)
  classifications <- classify_shared_dataset_only(prevalence)
  class_mapping <- build_gene_class_mapping(data)
  drug_classes <- summarise_drug_classes(functional, metadata)
  aminoglycoside <- build_aminoglycoside_outputs(data, metadata)
  fisher <- run_gene_fisher_tests(prevalence)
  singletons <- build_singleton_table(prevalence, functional)
  notable <- build_notable_difference_table(prevalence, fisher)
  contig_hits <- build_functional_hits_with_contigs(functional)
  qc <- build_qc_summary(data, matrix, prevalence, fisher, classifications)

  # Derive clustered orders and the complete distance audit before registration.
  clustered_heatmap <- build_clustered_presence_absence_heatmap(
    matrix, prevalence, functional, strain_metadata
  )

  # Register every table explicitly so missing or duplicated destinations fail safely.
  objects <- list(
    amr_standardised_long_72 = data,
    amr_flagged_low_confidence_72 = flagged,
    amr_method_audit_72 = method_audit,
    amr_functional_hits_with_contigs_72 = contig_hits,
    amr_gene_presence_absence = matrix,
    amr_shared_and_dataset_only_determinants = classifications,
    amr_gene_prevalence_72 = format_prevalence_table(prevalence),
    amr_gene_class_mapping_72 = class_mapping,
    amr_drug_class_summary_72 = drug_classes,
    aminoglycoside_determinants_72 = aminoglycoside$determinants,
    aminoglycoside_profile_by_genome_72 = aminoglycoside$profile,
    amr_gene_fisher_tests_72 = fisher,
    amr_singletons_72 = singletons,
    amr_notable_prevalence_differences_72 = notable,
    amr_pipeline_qc_summary_72 = qc,
    amr_presence_absence_heatmap_clustering_audit_72 = clustered_heatmap$audit
  )

  # Build the documented plots from validated in-memory analysis objects.
  plots <- list(
    gene_prevalence = plot_gene_prevalence(prevalence),
    presence_absence = plot_presence_absence_heatmap(
      matrix, prevalence, functional, strain_metadata
    ),
    presence_absence_clustered = clustered_heatmap$plot,
    drug_class = plot_drug_class_prevalence(drug_classes)
  )

  # Write all registered outputs only after every transformation succeeds.
  write_outputs(objects, plots, paths, args$overwrite)

  # Finish the audit log with result counts and full R session provenance.
  log_message(log_path, "Accepted functional AMR rows: ", nrow(functional))
  log_message(log_path, "Genes in functional matrix: ", nrow(prevalence))
  log_message(log_path, "Genes tested with Fisher's exact test: ", nrow(fisher))
  log_message(log_path, "All registered outputs written successfully")
  cat(
    paste(capture.output(utils::sessionInfo()), collapse = "\n"), "\n",
    file = log_path, append = TRUE, sep = ""
  )

  message("\nAnalysis complete. Statistical interpretation should prioritise ",
          "prevalence, direction and effect size alongside p-values and ",
          "FDR-adjusted values because the source groups are small and unbalanced ",
          "(n=14 versus n=58). Review odds ratios and prevalence differences even ",
          "when q >= 0.05, especially differences of at least 20 percentage points.")
}

if (sys.nframe() == 0L) {
  main()
}
