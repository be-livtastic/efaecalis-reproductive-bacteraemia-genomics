#!/usr/bin/env Rscript

# Create the writable directory skeleton used by the reproducible workflow.
# Run from anywhere; the repository root is resolved from this script.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

# --- Required writable directories ---
directories <- file.path(project_root, c(
  "data/raw",
  "data/interim",
  "data/processed",
  "analysis",
  "results/figures",
  "results/tables/annotation_qc",
  "results/tables/metadata_qc",
  "manuscript/supplementary"
))

# --- Create missing paths safely ---
for (directory in directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

cat("Project directories are ready under:\n", project_root, "\n", sep = "")
