#!/usr/bin/env Rscript

# Retrieve and pin the official Enterococcus faecalis PubMLST scheme resources.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

LOCI <- c("gdh", "gyd", "pstS", "gki", "aroE", "xpt", "yqiL")
BASE_URL <- "https://rest.pubmlst.org/db/pubmlst_efaecalis_seqdef"

# Resolve the project root from this script or an explicit environment override.
get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  derived <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  derived
}

# Parse the snapshot date and overwrite safeguard.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  date_arg <- grep("^--date=", args, value = TRUE)
  snapshot_date <- if (length(date_arg)) sub("^--date=", "", date_arg[[1]]) else format(Sys.Date(), "%Y-%m-%d")
  if (!str_detect(snapshot_date, "^\\d{4}-\\d{2}-\\d{2}$")) stop("--date must use YYYY-MM-DD", call. = FALSE)
  unknown <- args[!args %in% c("--overwrite") & !str_detect(args, "^--date=")]
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  list(date = snapshot_date, overwrite = "--overwrite" %in% args)
}

# Compute a SHA-256 checksum with the system utility already used by the project.
sha256_file <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(output)) stop("sha256sum failed for ", path, call. = FALSE)
  str_split(output[[1]], "\\s+", simplify = TRUE)[1, 1]
}

# Parse FASTA records without requiring another R package.
validate_alleles <- function(path, locus) {
  lines <- readLines(path, warn = FALSE)
  headers <- sub("^>", "", lines[str_detect(lines, "^>")]) |> word(1)
  if (!length(headers)) stop("Empty allele FASTA for ", locus, call. = FALSE)
  ids <- str_match(headers, "(?:^|_)(\\d+)$")[, 2]
  if (anyNA(ids) || anyDuplicated(ids)) stop("Invalid or duplicated allele identifiers for ", locus, call. = FALSE)
  starts <- which(str_detect(lines, "^>"))
  ends <- c(starts[-1] - 1L, length(lines))
  sequences <- map2_chr(starts, ends, ~ paste0(lines[(.x + 1L):.y], collapse = "") |> toupper())
  if (any(!str_detect(sequences, "^[ACGT]+$")) || anyDuplicated(sequences)) {
    stop("Invalid or duplicated official allele sequences for ", locus, call. = FALSE)
  }
  length(headers)
}

# Download all resources into a temporary directory before promotion.
main <- function() {
  config <- parse_args()
  root <- get_project_root()
  destination <- file.path(root, "references", "mlst", "pubmlst", paste0("efaecalis_", config$date))
  if (dir.exists(destination) && !config$overwrite) stop("Snapshot exists; refusing to overwrite: ", destination, call. = FALSE)
  if (dir.exists(destination) && config$overwrite) unlink(destination, recursive = TRUE)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("efaecalis_pubmlst_", tmpdir = dirname(destination))
  dir.create(temporary)
  on.exit(if (dir.exists(temporary)) unlink(temporary, recursive = TRUE), add = TRUE)

  resources <- map_dfr(LOCI, function(locus) {
    url <- paste0(BASE_URL, "/loci/", locus, "/alleles_fasta")
    filename <- paste0(locus, ".alleles.fasta")
    path <- file.path(temporary, filename)
    download.file(url, path, mode = "wb", quiet = TRUE, method = "libcurl")
    tibble(Resource = locus, File = filename, Endpoint = url, Records = validate_alleles(path, locus))
  })
  profile_url <- paste0(BASE_URL, "/schemes/1/profiles_csv")
  profile_path <- file.path(temporary, "profiles.tsv")
  download.file(profile_url, profile_path, mode = "wb", quiet = TRUE, method = "libcurl")
  profiles <- read_tsv(profile_path, col_types = cols(.default = col_character()), na = character())
  required <- c("ST", LOCI)
  if (!all(required %in% names(profiles))) stop("Profile response lacks: ", paste(setdiff(required, names(profiles)), collapse = ", "), call. = FALSE)
  if (!nrow(profiles) || anyDuplicated(profiles$ST) || anyDuplicated(profiles[LOCI])) stop("Official profile table is empty or duplicated", call. = FALSE)
  resources <- bind_rows(resources, tibble(Resource = "profiles", File = "profiles.tsv", Endpoint = profile_url, Records = nrow(profiles))) |>
    mutate(Retrieved_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), SHA256 = map_chr(file.path(temporary, File), sha256_file))
  write_tsv(resources, file.path(temporary, "provenance.tsv"))
  writeLines(paste(resources$SHA256, resources$File, sep = "  "), file.path(temporary, "SHA256SUMS"))
  schema <- list(
    scheme = "Enterococcus faecalis seven-locus MLST",
    locus_order = LOCI,
    profile_columns = names(profiles),
    curator_grouping_columns = names(profiles)[tolower(names(profiles)) %in% c("clonal_complex", "cc")]
  )
  dput(schema, file.path(temporary, "schema.dput"))
  if (!file.rename(temporary, destination)) stop("Could not promote validated snapshot to ", destination, call. = FALSE)
  message("Pinned PubMLST snapshot: ", destination)
  print(resources)
}

# Run only when invoked deliberately with Rscript.
if (!interactive() && Sys.getenv("EFAECALIS_SKIP_MAIN", "0") != "1") main()
