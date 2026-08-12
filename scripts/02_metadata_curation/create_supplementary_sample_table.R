#!/usr/bin/env Rscript

# Validate the selected sample catalog and create analysis metadata plus the
# manuscript-ready Supplementary Table S1. Uses base R only.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
cli <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% cli

# --- Canonical inputs and generated outputs ---
source_file <- file.path(project_root, "data", "metadata", "sample_metadata_source_72.tsv")
accession_file <- file.path(project_root, "data", "accession_lists", "selected_72_accessions.tsv")
existing_curated <- file.path(project_root, "data", "metadata", "curated_metadata_72_genomes.csv")
curated_output <- existing_curated
supplementary_output <- file.path(
  project_root, "manuscript", "supplementary", "table_s1_sample_metadata.tsv"
)
qc_output <- file.path(project_root, "results", "tables", "metadata_qc", "sample_identifier_qc.tsv")

# --- Read public metadata and the selected-accession contract ---
required_inputs <- c(source_file, accession_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Required input(s) missing:\n", paste(missing_inputs, collapse = "\n"))
}

source <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE,
                     na.strings = c("", "NA"), quote = "\"")
selected <- read.delim(accession_file, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"), quote = "\"")

# --- Validate assembly, BioSample, strain and study identifiers ---
required <- c(
  "assembly_accession", "biosample_accession", "bioproject_accession",
  "strain", "dataset_category", "study_title"
)
missing_columns <- setdiff(required, names(source))
if (length(missing_columns) > 0L) {
  stop("Sample catalog lacks required columns: ", paste(missing_columns, collapse = ", "))
}

trim <- function(x) trimws(as.character(x))
source$assembly_accession <- trim(source$assembly_accession)
source$biosample_accession <- trim(source$biosample_accession)
source$bioproject_accession <- trim(source$bioproject_accession)
source$strain <- trim(source$strain)
source$dataset_category <- trim(source$dataset_category)
selected$assembly_accession <- trim(selected$assembly_accession)
selected$dataset_category <- trim(selected$dataset_category)

issues <- character()
if (nrow(source) != 72L) issues <- c(issues, paste0("sample catalog row count is ", nrow(source), ", expected 72"))
if (anyDuplicated(source$assembly_accession)) issues <- c(issues, "assembly accessions are not unique")
if (anyDuplicated(source$biosample_accession)) issues <- c(issues, "BioSample accessions are not unique")
if (anyDuplicated(source$strain)) issues <- c(issues, "strain names are not unique")
if (any(!grepl("^GC[AF]_[0-9]+\\.[0-9]+$", source$assembly_accession))) issues <- c(issues, "invalid assembly accession format")
if (any(!grepl("^SAM[NED][A-Z]?[0-9]+$", source$biosample_accession))) issues <- c(issues, "invalid BioSample accession format")
if (any(!grepl("^PRJ(NA|EB|DB)[0-9]+$", source$bioproject_accession))) issues <- c(issues, "invalid BioProject accession format")
if (any(is.na(source$strain) | source$strain == "")) issues <- c(issues, "missing strain name")

source_key <- paste(source$assembly_accession, source$dataset_category, sep = "|")
selected_key <- paste(selected$assembly_accession, selected$dataset_category, sep = "|")
if (!setequal(source_key, selected_key)) {
  issues <- c(issues, "sample catalog does not exactly match the selected accession/category manifest")
}
group_counts <- table(source$dataset_category)
if (!identical(as.integer(group_counts[c("Bacteraemia", "Reproductive")]), c(58L, 14L))) {
  issues <- c(issues, "dataset counts are not 58 Bacteraemia and 14 Reproductive")
}

if (length(issues) > 0L) stop(paste(issues, collapse = "\n"))

# --- Preserve existing analysis fields when rebuilding metadata ---
amr_columns <- c(
  "assembly_accession", "amr_gene_count", "hlgr_gentamicin_marker_status",
  "hlgr_gentamicin_markers"
)
if (file.exists(existing_curated)) {
  old <- read.csv(existing_curated, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
  available <- intersect(amr_columns, names(old))
  old <- old[, available, drop = FALSE]
  source <- merge(source, old, by = "assembly_accession", all.x = TRUE, sort = FALSE)
  source <- source[match(selected$assembly_accession, source$assembly_accession), , drop = FALSE]
}

# --- Build analysis and manuscript views from one source ---
curated_columns <- intersect(c(
  "assembly_accession", "biosample_accession", "bioproject_accession", "strain",
  "dataset_category", "source_group", "sample_type", "country",
  "geographic_location", "collection_year", "host", "study_title",
  "amr_gene_count", "hlgr_gentamicin_marker_status", "hlgr_gentamicin_markers"
), names(source))
curated <- source[, curated_columns, drop = FALSE]
names(curated)[names(curated) == "dataset_category"] <- "reproductive_bacteraemia_category"
if ("collection_year" %in% names(curated)) names(curated)[names(curated) == "collection_year"] <- "year"

source$sample_number <- seq_len(nrow(source))
source$dataset_id <- ifelse(
  source$dataset_category == "Reproductive", "reproductive_14", "bacteraemia_58"
)
supplementary_columns <- intersect(c(
  "sample_number", "dataset_id", "dataset_category", "strain", "biosample_accession",
  "assembly_accession", "bioproject_accession", "study_title", "assembly_level",
  "source_group", "sample_type", "host", "geographic_location", "country",
  "collection_year", "sequencing_technology", "chromosome_count",
  "plasmid_status", "plasmid_count", "associated_paper", "pathology_context"
), names(source))
supplementary <- source[, supplementary_columns, drop = FALSE]

# --- Write outputs only after every validation passes ---
outputs <- c(curated_output, supplementary_output, qc_output)
existing <- outputs[file.exists(outputs)]
if (length(existing) > 0L && !overwrite) {
  stop("Refusing to overwrite existing output(s); rerun after review with --overwrite:\n",
       paste(existing, collapse = "\n"))
}
for (directory in unique(dirname(outputs))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

write.csv(curated, curated_output, row.names = FALSE, na = "")
write.table(supplementary, supplementary_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
qc <- data.frame(
  check = c("samples", "unique_assemblies", "unique_biosamples", "unique_strains",
            "reproductive_samples", "bacteraemia_samples", "manifest_exact_match"),
  value = c(nrow(source), length(unique(source$assembly_accession)),
            length(unique(source$biosample_accession)), length(unique(source$strain)),
            unname(group_counts[["Reproductive"]]), unname(group_counts[["Bacteraemia"]]), TRUE)
)
write.table(qc, qc_output, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Validated 72 assembly/BioSample/strain mappings.\n")
cat("Curated metadata: ", curated_output, "\n", sep = "")
cat("Supplementary Table S1: ", supplementary_output, "\n", sep = "")
cat("Identifier QC: ", qc_output, "\n", sep = "")
