#!/usr/bin/env Rscript

# Assign formal seven-locus PubMLST sequence types to the canonical 72 genomes.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

LOCI <- c("gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL")
EXPECTED_GROUPS <- c(Reproductive = 14L, Bacteraemia = 58L)
BLAST_COLUMNS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore", "qseq", "sseq")

# Resolve repository-relative inputs and outputs.
get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  derived <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  derived
}

# Parse the pinned reference directory and overwrite option.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  ref_index <- match("--reference-dir", args)
  if (is.na(ref_index) || ref_index == length(args)) stop("Usage: Rscript ... --reference-dir PATH [--overwrite]", call. = FALSE)
  reference <- args[[ref_index + 1L]]
  consumed <- c(ref_index, ref_index + 1L, which(args == "--overwrite"))
  if (length(setdiff(seq_along(args), consumed))) stop("Unknown command-line argument", call. = FALSE)
  list(reference = normalizePath(reference, mustWork = TRUE), overwrite = "--overwrite" %in% args)
}

# Refuse mixed-snapshot outputs unless overwrite is explicit.
preflight_outputs <- function(paths, overwrite) {
  if (anyDuplicated(normalizePath(paths, mustWork = FALSE))) stop("Duplicate output path in registry", call. = FALSE)
  existing <- paths[file.exists(paths)]
  if (length(existing) && !overwrite) stop("Refusing to overwrite: ", paste(existing, collapse = ", "), call. = FALSE)
  walk(unique(dirname(paths)), dir.create, recursive = TRUE, showWarnings = FALSE)
}

# Validate 72 versioned accessions and one canonical FNA per manifest row.
build_manifest <- function(root) {
  metadata <- read_tsv(file.path(root, "data", "accession_lists", "selected_72_accessions.tsv"), col_types = cols(.default = col_character())) |>
    transmute(Genome = assembly_accession, Source = dataset_category)
  if (nrow(metadata) != 72L || n_distinct(metadata$Genome) != 72L) stop("Manifest must contain 72 unique accessions", call. = FALSE)
  observed_groups <- table(metadata$Source)
  if (!setequal(names(observed_groups), names(EXPECTED_GROUPS)) ||
      any(as.integer(observed_groups[names(EXPECTED_GROUPS)]) != unname(EXPECTED_GROUPS))) {
    stop("Manifest does not contain the expected 14/58 groups", call. = FALSE)
  }
  dirs <- c(
    Reproductive = file.path(root, "local_archive", "large_outputs", "Efaecalis_14_AMRFinder", "02_genomes_fna"),
    Bacteraemia = file.path(root, "local_archive", "large_outputs", "Efaecalis_58_AMRFinder", "02_genomes_fna")
  )
  files <- imap_dfr(dirs, function(path, source) {
    fna <- sort(list.files(path, pattern = "\\.fna$", full.names = TRUE))
    if (length(fna) != EXPECTED_GROUPS[[source]]) stop("Unexpected FNA count in ", path, call. = FALSE)
    tibble(Genome = str_remove(basename(fna), "\\.fna$"), Source_from_path = source, FNA_path = normalizePath(fna))
  })
  if (any(!str_detect(files$Genome, "^GC[AF]_\\d+\\.\\d+$")) || n_distinct(files$Genome) != 72L) stop("FNA filenames must encode 72 versioned accessions", call. = FALSE)
  joined <- full_join(metadata, files, by = "Genome")
  if (nrow(joined) != 72L || anyNA(joined$FNA_path) || any(joined$Source != joined$Source_from_path)) stop("FNA/metadata join failed", call. = FALSE)
  joined |> select(-Source_from_path) |> arrange(Genome)
}

# Parse the numeric identifier from each official allele header.
parse_allele <- function(header) {
  value <- str_match(header, "(?:^|_)(\\d+)$")[, 2]
  if (anyNA(value)) stop("Unparseable PubMLST allele header: ", paste(header[is.na(value)], collapse = ", "), call. = FALSE)
  value
}

# Execute one allele-library versus genome BLASTN comparison.
run_blast <- function(allele_fasta, genome_fna) {
  output <- tempfile(fileext = ".tsv")
  on.exit(unlink(output), add = TRUE)
  args <- c("-task", "blastn", "-dust", "no", "-query", shQuote(allele_fasta), "-subject", shQuote(genome_fna),
            "-evalue", "1e-20", "-perc_identity", "80", "-qcov_hsp_perc", "80", "-max_target_seqs", "20",
            "-outfmt", shQuote(paste("6", paste(BLAST_COLUMNS, collapse = " "))), "-out", shQuote(output))
  status <- system2("blastn", args = args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) stop("blastn failed for ", genome_fna, ": ", paste(status, collapse = " "), call. = FALSE)
  if (!file.size(output)) return(tibble())
  read_tsv(output, col_names = BLAST_COLUMNS, col_types = cols(.default = col_character()), show_col_types = FALSE) |>
    mutate(
      across(c(pident, evalue, bitscore), as.numeric),
      across(c(length, mismatch, gapopen, gaps, qstart, qend, sstart, send, qlen, slen), as.integer),
      Locus_from_query = str_extract(qseqid, "^[^_]+"), Official_allele = parse_allele(qseqid), Allele_coverage_percent = 100 * length / qlen,
      Contig = sseqid, Start = pmin(sstart, send), Stop = pmax(sstart, send),
      Strand = if_else(sstart <= send, "+", "-"), Location_key = paste(Contig, ((Start + Stop) %/% 2L) %/% 100L, sep = ":")
    )
}

# Create a stable SHA-256 token for an observed novel sequence.
sequence_token <- function(sequence) {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  writeLines(sequence, path, useBytes = TRUE)
  checksum <- system2("sha256sum", shQuote(path), stdout = TRUE)
  paste0("NOVEL_", substr(word(checksum[[1]], 1), 1, 12))
}

# Apply the exact known-allele and explicit unresolved-state policy.
classify_locus <- function(genome, source, locus, hits, fna_path) {
  exact <- hits |> filter(pident == 100, length == qlen, gaps == 0)
  exact_distinct <- exact |> distinct(Official_allele, Location_key, .keep_all = TRUE)
  exact_alleles <- sort(unique(exact_distinct$Official_allele))
  exact_locations <- n_distinct(exact_distinct$Location_key)
  accepted <- token <- ""
  if (length(exact_alleles) > 1L) {
    status <- "Multiple_exact_alleles"; reason <- paste("Distinct exact official alleles:", paste(exact_alleles, collapse = ";"))
  } else if (length(exact_alleles) == 1L && exact_locations > 1L) {
    status <- "Multiple_locus_copies"; accepted <- token <- exact_alleles[[1]]; reason <- paste("Same exact allele at", exact_locations, "locations")
  } else if (length(exact_alleles) == 1L) {
    status <- "Exact_known_allele"; accepted <- token <- exact_alleles[[1]]; reason <- "One unambiguous exact full-length official allele"
  } else {
    full <- hits |> filter(Allele_coverage_percent >= 99.999, pident >= 95)
    if (nrow(full)) {
      best <- full |> filter(bitscore >= max(bitscore) * 0.999)
      if (n_distinct(best$Location_key) == 1L) {
        chosen <- best |> arrange(desc(pident), gaps, as.integer(Official_allele)) |> slice(1)
        observed <- str_replace_all(toupper(chosen$sseq[[1]]), "[^ACGT]", "")
        if (nchar(observed)) {
          status <- "Novel_allele_candidate"; token <- sequence_token(observed); reason <- "One complete high-similarity locus sequence has no exact official allele"
        } else {
          status <- "Unresolved"; reason <- "Complete hit did not yield an unambiguous observed sequence"
        }
      } else {
        status <- "Unresolved"; reason <- paste("Best complete non-exact hits occur at", n_distinct(best$Location_key), "locations")
      }
    } else if (any(hits$Allele_coverage_percent >= 80 & hits$pident >= 80)) {
      status <- "Non_exact_best_match"; reason <- "Only incomplete or non-qualifying nearest official allele matches"
    } else if (nrow(hits)) {
      status <- "Unresolved"; reason <- "Hits exist below the minimum audit threshold"
    } else {
      status <- "Missing_locus"; reason <- "No PubMLST allele-reference hit detected"
    }
  }
  best <- if (nrow(hits)) hits |> arrange(desc(bitscore), desc(pident)) |> slice(1) else tibble()
  summary <- tibble(
    Genome = genome, Source = source, Locus = locus, Assignment_status = status,
    Accepted_allele = accepted, Dataset_profile_token = token,
    Exact_alleles_detected = paste(exact_alleles, collapse = ";"), Exact_location_count = exact_locations,
    Best_reference_allele = if (nrow(best)) best$Official_allele else NA_character_,
    Best_identity_percent = if (nrow(best)) best$pident else NA_real_,
    Best_allele_coverage_percent = if (nrow(best)) best$Allele_coverage_percent else NA_real_,
    Best_contig = if (nrow(best)) best$Contig else NA_character_, Best_start = if (nrow(best)) best$Start else NA_integer_,
    Best_stop = if (nrow(best)) best$Stop else NA_integer_, Best_strand = if (nrow(best)) best$Strand else NA_character_,
    Status_reason = reason, FNA_path = fna_path
  )
  retained_exact <- hits |> filter(pident == 100 & length == qlen & gaps == 0)
  retained_near <- hits |> filter(!(pident == 100 & length == qlen & gaps == 0)) |>
    arrange(desc(bitscore), desc(pident), desc(Allele_coverage_percent)) |> distinct(Official_allele, Location_key, .keep_all = TRUE) |> slice_head(n = 20)
  retained <- bind_rows(retained_exact, retained_near) |> distinct(qseqid, Contig, Start, Stop, .keep_all = TRUE) |>
    transmute(Genome = genome, Source = source, Locus = locus, Official_allele, Contig, Start, Stop, Strand,
              Identity_percent = pident, Allele_coverage_percent, Alignment_length = length, Mismatches = mismatch,
              Gaps = gaps, Evalue = evalue, Bitscore = bitscore,
              Is_exact_full_length = pident == 100 & length == qlen & gaps == 0)
  list(summary = summary, audit = retained)
}

# Match exact seven-locus combinations to official profiles only.
build_assignments <- function(locus_qc, manifest, profiles, grouping_column) {
  wide <- locus_qc |>
    select(Genome, Locus, Accepted_allele, Dataset_profile_token, Assignment_status) |>
    pivot_wider(names_from = Locus, values_from = c(Accepted_allele, Dataset_profile_token, Assignment_status))
  map_dfr(manifest$Genome, function(genome) {
    row <- wide |> filter(Genome == genome)
    alleles <- set_names(map_chr(LOCI, ~ row[[paste0("Accepted_allele_", .x)]]), LOCI)
    tokens <- set_names(map_chr(LOCI, ~ row[[paste0("Dataset_profile_token_", .x)]]), LOCI)
    statuses <- set_names(map_chr(LOCI, ~ row[[paste0("Assignment_status_", .x)]]), LOCI)
    formal_complete <- all(!is.na(alleles) & alleles != "")
    matched <- if (formal_complete) {
      profile_matches <- Reduce(`&`, map2(profiles[LOCI], alleles, ~ .x == .y))
      profiles[profile_matches, , drop = FALSE]
    } else profiles[0, , drop = FALSE]
    if (nrow(matched) == 1L) {
      mlst_status <- "Assigned"; st <- matched$ST[[1]]; reason <- "Exact seven-locus profile matched one official PubMLST ST"
      cc <- if (!is.null(grouping_column)) matched[[grouping_column]][[1]] else NA_character_
    } else if (formal_complete) {
      mlst_status <- "Unassigned_profile_not_in_PubMLST"; st <- "Unassigned"; reason <- "Exact known-allele combination absent from pinned profile table"; cc <- NA_character_
    } else {
      mlst_status <- "Unassigned_incomplete_profile"; st <- "Unassigned"
      reason <- paste(paste0(names(statuses)[alleles == "" | is.na(alleles)], ":", statuses[alleles == "" | is.na(alleles)]), collapse = "; "); cc <- NA_character_
    }
    tibble(
      Genome = genome, Source = manifest$Source[match(genome, manifest$Genome)],
      !!!as.list(alleles),
      Allele_profile = paste0(LOCI, ifelse(alleles == "" | is.na(alleles), "UNRESOLVED", alleles), collapse = "-"),
      Dataset_profile = paste0(LOCI, ifelse(tokens == "" | is.na(tokens), "UNRESOLVED", tokens), collapse = "-"),
      Dataset_profile_complete = all(!is.na(tokens) & tokens != ""), MLST_status = mlst_status,
      MLST_status_reason = reason, ST = st, PubMLST_CC_if_officially_available = cc
    )
  })
}

# Run all 504 comparisons and promote outputs only after validation.
main <- function() {
  args <- parse_args(); root <- get_project_root(); reference <- args$reference
  if (!nzchar(Sys.which("blastn"))) stop("blastn is unavailable in PATH", call. = FALSE)
  profiles <- read_tsv(file.path(reference, "profiles.tsv"), col_types = cols(.default = col_character()), na = character())
  grouping <- names(profiles)[tolower(names(profiles)) %in% c("clonal_complex", "cc")]
  grouping <- if (length(grouping)) grouping[[1]] else NULL
  if (anyDuplicated(profiles$ST) || anyDuplicated(profiles[LOCI])) stop("Pinned profile table is duplicated", call. = FALSE)
  manifest <- build_manifest(root)
  output_dir <- file.path(root, "data", "processed", "mlst")
  outputs <- file.path(output_dir, c("mlst_allele_hit_audit_72.tsv", "mlst_locus_assignment_qc_72.csv", "formal_mlst_assignments_72.csv", "mlst_assignment_qc_summary_72.csv"))
  preflight_outputs(outputs, args$overwrite)

  # Combine the seven official query libraries so each genome is searched only once.
  combined_alleles <- tempfile(fileext = ".fasta")
  on.exit(unlink(combined_alleles), add = TRUE)
  writeLines(unlist(map(file.path(reference, paste0(LOCI, ".alleles.fasta")), readLines, warn = FALSE)), combined_alleles)

  # Run one combined BLASTN search per genome, then apply locus-specific classification.
  results <- pmap(manifest, function(Genome, Source, FNA_path) {
    message("MLST BLASTN: ", Genome)
    genome_hits <- run_blast(combined_alleles, FNA_path)
    map(LOCI, function(locus) classify_locus(Genome, Source, locus, filter(genome_hits, Locus_from_query == locus), FNA_path))
  }) |> flatten()
  locus_qc <- map_dfr(results, "summary"); audit <- map_dfr(results, "audit")
  if (nrow(locus_qc) != 504L || anyDuplicated(locus_qc[c("Genome", "Locus")])) stop("Expected 504 unique genome-locus summaries", call. = FALSE)
  assignments <- build_assignments(locus_qc, manifest, profiles, grouping)
  if (nrow(assignments) != 72L || n_distinct(assignments$Genome) != 72L) stop("Expected 72 unique assignment rows", call. = FALSE)
  qc <- bind_rows(
    tibble(Metric = c("genomes", "reproductive_genomes", "bacteraemia_genomes", "genome_locus_records", "formal_st_assigned", "formal_st_unassigned", "complete_dataset_profiles", "pubmlst_curator_grouping_column"),
           Value = as.character(c(72, 14, 58, 504, sum(assignments$MLST_status == "Assigned"), sum(assignments$MLST_status != "Assigned"), sum(assignments$Dataset_profile_complete), if (is.null(grouping)) "not_present" else grouping))),
    count(locus_qc, Assignment_status) |> transmute(Metric = paste0("locus_status_", Assignment_status), Value = as.character(n))
  )
  write_tsv(audit, outputs[[1]]); write_csv(locus_qc, outputs[[2]]); write_csv(assignments, outputs[[3]]); write_csv(qc, outputs[[4]])
  print(qc, n = Inf)
}

# Run only during a deliberate Rscript invocation.
if (!interactive() && Sys.getenv("EFAECALIS_SKIP_MAIN", "0") != "1") main()
