#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

get_project_root <- function() {
  override <- Sys.getenv("EFAECALIS_PROJECT_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Run this script via Rscript so the project root can be inferred.", call. = FALSE)
  }
  script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
}

van_columns <- c(
  "vanA", "vanB", "vanH-A", "vanH-B", "vanR-A", "vanR-B",
  "vanS-A", "vanS-B", "vanW-B", "vanX-A", "vanX-B",
  "vanY-A", "vanY-B", "vanZ-A"
)

main <- function() {
  root <- get_project_root()
  metadata_path <- file.path(root, "data", "metadata", "sample_metadata_source_72.tsv")
  plasmid_path <- file.path(root, "data", "metadata", "plasmid_status_by_genome_72.tsv")
  amr_path <- file.path(root, "results", "tables", "amr", "amr_gene_presence_absence.csv")
  out_dir <- file.path(root, "results", "figures", "amr")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- read_tsv(metadata_path, na = c("", "NA"), show_col_types = FALSE) |>
    transmute(
      assembly_accession = assembly_accession,
      strain = strain,
      dataset_category = dataset_category,
      source_group = if_else(is.na(source_group), "Unknown", source_group)
    )

  plasmid <- read_tsv(plasmid_path, na = c("", "NA"), show_col_types = FALSE) |>
    transmute(
      assembly_accession = assembly_accession,
      strain = strain,
      dataset_category = dataset_category,
      plasmid_status = if_else(plasmid_status == "Y", "Present", "Absent")
    )

  amr <- read_csv(amr_path, show_col_types = FALSE)

  if (!all(c("Genome", van_columns) %in% names(amr))) {
    stop("Expected vancomycin-associated AMR columns are missing.", call. = FALSE)
  }

  plot_df <- metadata |>
    left_join(plasmid, by = c("assembly_accession", "strain", "dataset_category"), relationship = "many-to-one") |>
    left_join(
      amr |>
        transmute(
          Genome,
          van_any = if_else(
            rowSums(across(all_of(van_columns), ~ as.integer(.x) > 0), na.rm = TRUE) > 0,
            "Van genes present",
            "No van genes"
          )
        ),
      by = c("assembly_accession" = "Genome")
    ) |>
    mutate(
      strain = factor(
        strain,
        levels = unique(strain[order(dataset_category, source_group, strain)])
      ),
      dataset_category = factor(dataset_category, levels = c("Reproductive", "Bacteraemia")),
      plasmid_status = coalesce(plasmid_status, "Absent"),
      source_group = if_else(is.na(source_group), "Unknown", source_group),
      source_label = case_when(
        source_group == "Semen sample" ~ "Semen",
        source_group == "Vaginal swab" ~ "Vaginal",
        source_group == "Clinical laboratory" ~ "Clinic",
        TRUE ~ str_replace_all(source_group, " ", " ")
      )
    )

  plot_long <- plot_df |>
    select(strain, dataset_category, plasmid_status, source_group, van_any) |>
    mutate(
      dataset_category = if_else(is.na(dataset_category), "Unknown", as.character(dataset_category)),
      source_group = if_else(is.na(source_group), "Unknown", source_group)
    ) |>
    pivot_longer(
      cols = c(dataset_category, plasmid_status, source_group, van_any),
      names_to = "feature",
      values_to = "value"
    ) |>
    mutate(
      feature = factor(
        feature,
        levels = c("dataset_category", "plasmid_status", "source_group", "van_any"),
        labels = c("Dataset", "Plasmid", "Source", "Van genes")
      ),
      value_label = case_when(
        feature == "Dataset" ~ value,
        feature == "Plasmid" ~ if_else(value == "Present", "Plasmid", "No plasmid"),
        feature == "Source" ~ case_when(
          value == "Semen sample" ~ "Semen",
          value == "Vaginal swab" ~ "Vaginal",
          value == "Clinical laboratory" ~ "Clinic",
          TRUE ~ value
        ),
        feature == "Van genes" ~ if_else(value == "Van genes present", "Van", "No van")
      )
    )

  fill_palette <- c(
    "Reproductive" = "#4C78A8",
    "Bacteraemia" = "#F58518",
    "Present" = "#2E7D32",
    "Absent" = "#D1D5DB",
    "Semen sample" = "#66C2A5",
    "Vaginal swab" = "#FC8D62",
    "Clinical laboratory" = "#8DA0CB",
    "Unknown" = "#BDBDBD",
    "Van genes present" = "#7B1FA2",
    "No van genes" = "#E5E7EB"
  )

  plot <- ggplot(plot_long, aes(x = feature, y = strain)) +
    geom_tile(aes(fill = value), colour = "white", linewidth = 0.15) +
    geom_text(aes(label = value_label), size = 2.5) +
    scale_fill_manual(values = fill_palette, guide = "none") +
    labs(
      title = "Strain-level plasmid status, source and vancomycin-associated genes",
      subtitle = "Only observed metadata and AMR presence/absence values were used; no inferred plasmid-gene mappings were added.",
      x = NULL,
      y = "Strain"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 5.5),
      panel.grid = element_blank(),
      plot.margin = margin(0.35, 0.6, 0.35, 0.2, "cm")
    ) +
    coord_fixed(ratio = 0.38)

  out_path <- file.path(out_dir, "strain_plasmid_source_van_72.png")
  ggsave(out_path, plot, width = 13.5, height = 18, dpi = 300, bg = "white")
  message("Saved plot to ", out_path)
}

if (!interactive()) main()
