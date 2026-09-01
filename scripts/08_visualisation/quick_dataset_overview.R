# =============================================================================
# Quick overview of the E. faecalis phylogeny dataset
# =============================================================================

# -----------------------------
# 1. Install and load packages
# -----------------------------
required_packages <- c(
  "tidyverse",
  "janitor",
  "scales"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

invisible(lapply(required_packages, library, character.only = TRUE))

# -----------------------------
# 2. Portable file paths
# -----------------------------
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath(
    file.path(getwd(), "scripts", "08_visualisation",
              "quick_dataset_overview.R"),
    mustWork = TRUE
  )
}

# EFAECALIS_PROJECT_ROOT supports non-standard installations without storing a
# personal Windows, WSL or Linux path in version control.
PROJECT_DIR <- Sys.getenv(
  "EFAECALIS_PROJECT_ROOT",
  unset = normalizePath(
    file.path(dirname(script_path), "..", ".."),
    mustWork = TRUE
  )
)
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
ALLOW_OVERWRITE <- tolower(
  Sys.getenv("EFAECALIS_ALLOW_OVERWRITE", unset = "false")
) %in% c("1", "true", "yes")

INPUT_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "metadata",
  "curated_metadata_72_genomes.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "figures",
  "dataset_overview"
)

if (dir.exists(OUTPUT_DIR) &&
    length(list.files(OUTPUT_DIR, all.files = TRUE, no.. = TRUE)) > 0 &&
    !ALLOW_OVERWRITE) {
  stop(
    "Output directory is not empty: ", OUTPUT_DIR, "\n",
    "Set EFAECALIS_ALLOW_OVERWRITE=true only after reviewing existing files."
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Metadata file was not found:\n",
    INPUT_FILE,
    "\nCheck the project-folder spelling and filename."
  )
}

# -----------------------------
# 3. Helper functions
# -----------------------------
first_existing_column <- function(data, candidates) {
  found <- candidates[candidates %in% names(data)]

  if (length(found) == 0) {
    return(NULL)
  }

  found[[1]]
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

normalise_binary_status <- function(x) {
  value <- stringr::str_to_lower(clean_text(x))

  dplyr::case_when(
    is.na(value) ~ "Unknown",

    stringr::str_detect(
      value,
      "absent|negative|no|false|not detected"
    ) ~ "Absent",

    stringr::str_detect(
      value,
      "present|positive|yes|true|detected|hlgr|gentamicin"
    ) ~ "Present",

    TRUE ~ stringr::str_to_title(value)
  )
}

save_plot <- function(plot, filename, width = 9, height = 6) {
  ggsave(
    filename = file.path(OUTPUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

# -----------------------------
# 4. Import and standardise data
# -----------------------------
raw_data <- readr::read_csv(
  INPUT_FILE,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "Unknown", "unknown")
) |>
  janitor::clean_names()

assembly_col <- first_existing_column(
  raw_data,
  c(
    "assembly_accession",
    "assembly",
    "genome_id",
    "assembly_key"
  )
)

category_col <- first_existing_column(
  raw_data,
  c(
    "dataset_category",
    "category",
    "source_category",
    "dataset",
    "group"
  )
)

year_col <- first_existing_column(
  raw_data,
  c(
    "collection_year",
    "year",
    "isolation_year",
    "sample_year"
  )
)

source_col <- first_existing_column(
  raw_data,
  c(
    "source_group",
    "isolation_source",
    "sample_source",
    "source",
    "sample_type"
  )
)

amr_count_col <- first_existing_column(
  raw_data,
  c(
    "amr_gene_count",
    "functional_amr_gene_count",
    "unique_amr_gene_count",
    "functional_amr_hit_count",
    "amr_count"
  )
)

hlgr_col <- first_existing_column(
  raw_data,
  c(
    "hlgr_gentamicin_marker_status",
    "hlgr_status",
    "gentamicin_marker_status",
    "gentamicin_status",
    "hlgr_detected",
    "gentamicin_detected"
  )
)

data <- tibble(
  assembly_accession = if (!is.null(assembly_col)) {
    clean_text(raw_data[[assembly_col]])
  } else {
    paste0("Genome_", seq_len(nrow(raw_data)))
  },

  dataset_category = if (!is.null(category_col)) {
    clean_text(raw_data[[category_col]])
  } else {
    NA_character_
  },

  collection_year = if (!is.null(year_col)) {
    suppressWarnings(as.integer(raw_data[[year_col]]))
  } else {
    NA_integer_
  },

  isolation_source = if (!is.null(source_col)) {
    clean_text(raw_data[[source_col]])
  } else {
    NA_character_
  },

  amr_gene_count = if (!is.null(amr_count_col)) {
    suppressWarnings(as.numeric(raw_data[[amr_count_col]]))
  } else {
    NA_real_
  },

  hlgr_gentamicin_status = if (!is.null(hlgr_col)) {
    normalise_binary_status(raw_data[[hlgr_col]])
  } else {
    "Unknown"
  }
)

data <- data |>
  mutate(
    dataset_category = case_when(
      str_detect(str_to_lower(dataset_category), "repro") ~ "Reproductive",
      str_detect(str_to_lower(dataset_category), "bacter|blood") ~ "Bacteraemia",
      TRUE ~ dataset_category
    ),
    dataset_category = replace_na(dataset_category, "Unknown"),
    isolation_source = replace_na(isolation_source, "Unknown"),
    hlgr_gentamicin_status = replace_na(
      hlgr_gentamicin_status,
      "Unknown"
    )
  )

readr::write_csv(
  data,
  file.path(OUTPUT_DIR, "dataset_overview_cleaned.csv")
)

# -----------------------------
# 5. Dataset summary table
# -----------------------------
valid_years <- data$collection_year[
  !is.na(data$collection_year)
]

valid_amr <- data$amr_gene_count[
  !is.na(data$amr_gene_count)
]

summary_table <- tibble(
  metric = c(
    "Total genomes",
    "Reproductive genomes",
    "Bacteraemia genomes",
    "Other or unknown category genomes",
    "Earliest collection year",
    "Latest collection year",
    "Genomes with collection year",
    "Genomes with AMR gene count",
    "Distinct AMR gene-count values",
    "Mean AMR genes per genome",
    "Median AMR genes per genome",
    "Minimum AMR genes per genome",
    "Maximum AMR genes per genome",
    "HLGR/gentamicin marker present",
    "HLGR/gentamicin marker absent",
    "HLGR/gentamicin status unknown"
  ),

  value = c(
    nrow(data),

    sum(
      data$dataset_category == "Reproductive",
      na.rm = TRUE
    ),

    sum(
      data$dataset_category == "Bacteraemia",
      na.rm = TRUE
    ),

    sum(
      !data$dataset_category %in%
        c("Reproductive", "Bacteraemia"),
      na.rm = TRUE
    ),

    if (length(valid_years) > 0) {
      min(valid_years)
    } else {
      NA
    },

    if (length(valid_years) > 0) {
      max(valid_years)
    } else {
      NA
    },

    sum(!is.na(data$collection_year)),

    sum(!is.na(data$amr_gene_count)),

    if (length(valid_amr) > 0) {
      dplyr::n_distinct(valid_amr)
    } else {
      NA
    },

    if (length(valid_amr) > 0) {
      round(mean(valid_amr), 2)
    } else {
      NA
    },

    if (length(valid_amr) > 0) {
      round(median(valid_amr), 2)
    } else {
      NA
    },

    if (length(valid_amr) > 0) {
      min(valid_amr)
    } else {
      NA
    },

    if (length(valid_amr) > 0) {
      max(valid_amr)
    } else {
      NA
    },

    sum(
      data$hlgr_gentamicin_status == "Present",
      na.rm = TRUE
    ),

    sum(
      data$hlgr_gentamicin_status == "Absent",
      na.rm = TRUE
    ),

    sum(
      data$hlgr_gentamicin_status == "Unknown",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  summary_table,
  file.path(OUTPUT_DIR, "dataset_summary.csv")
)

# -----------------------------
# 6. AMR gene-count frequency table
# -----------------------------
amr_frequency <- data |>
  filter(!is.na(amr_gene_count)) |>
  count(
    amr_gene_count,
    dataset_category,
    name = "number_of_genomes"
  ) |>
  arrange(amr_gene_count, dataset_category)

readr::write_csv(
  amr_frequency,
  file.path(OUTPUT_DIR, "amr_gene_count_frequency.csv")
)

# -----------------------------
# 7. Graph 1: genomes by dataset category
# -----------------------------
category_counts <- data |>
  count(dataset_category, name = "genome_count") |>
  arrange(genome_count)

plot_1 <- ggplot(
  category_counts,
  aes(
    x = reorder(dataset_category, genome_count),
    y = genome_count,
    fill = dataset_category
  )
) +
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = genome_count),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Reproductive" = "#7B2CBF",
      "Bacteraemia" = "#D62828",
      "Unknown" = "#6C757D"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    breaks = scales::pretty_breaks()
  ) +
  labs(
    title = "Genomes by dataset category",
    subtitle = "Number of E. faecalis genomes in each comparison group",
    x = NULL,
    y = "Number of genomes"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  plot_1,
  "01_genomes_by_dataset_category.png"
)

# -----------------------------
# 8. Graph 2: collection-year distribution
# -----------------------------
if (sum(!is.na(data$collection_year)) > 0) {
  plot_2 <- data |>
    filter(!is.na(collection_year)) |>
    ggplot(
      aes(
        x = collection_year,
        fill = dataset_category
      )
    ) +
    geom_histogram(
      binwidth = 1,
      boundary = 0,
      closed = "left",
      position = "stack"
    ) +
    scale_fill_manual(
      values = c(
        "Reproductive" = "#7B2CBF",
        "Bacteraemia" = "#F77F00",
        "Unknown" = "#6C757D"
      )
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks()
    ) +
    labs(
      title = "Collection-year distribution",
      subtitle = "Number of genomes collected in each year",
      x = "Collection year",
      y = "Number of genomes",
      fill = "Dataset category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )

  save_plot(
    plot_2,
    "02_collection_year_distribution.png"
  )
}

# -----------------------------
# 9. Graph 3: ranked AMR gene counts
# -----------------------------
if (sum(!is.na(data$amr_gene_count)) > 0) {
  ranked_amr <- data |>
    filter(!is.na(amr_gene_count)) |>
    arrange(dataset_category, amr_gene_count, assembly_accession) |>
    group_by(dataset_category) |>
    mutate(genome_rank = row_number()) |>
    ungroup()

  plot_3 <- ggplot(
    ranked_amr,
    aes(
      x = genome_rank,
      y = amr_gene_count,
      colour = dataset_category
    )
  ) +
    geom_segment(
      aes(
        xend = genome_rank,
        y = 0,
        yend = amr_gene_count
      ),
      linewidth = 0.7,
      alpha = 0.65
    ) +
    geom_point(
      size = 2.8
    ) +
    facet_wrap(
      ~ dataset_category,
      scales = "free_x"
    ) +
    scale_colour_manual(
      values = c(
        "Reproductive" = "#4361EE",
        "Bacteraemia" = "#E63946",
        "Unknown" = "#6C757D"
      )
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(),
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(
      title = "Ranked AMR gene counts by genome",
      subtitle = paste(
        "Each point represents one genome, ordered from lowest",
        "to highest AMR gene count within each dataset"
      ),
      x = "Genome rank within dataset",
      y = "AMR gene count per genome",
      colour = "Dataset category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )

  save_plot(
    plot_3,
    "03_ranked_amr_gene_counts.png",
    width = 11,
    height = 6
  )
}

# -----------------------------
# 10. Graph 4: AMR burden by category
# -----------------------------
if (
  sum(!is.na(data$amr_gene_count)) > 0 &&
  n_distinct(data$dataset_category) > 1
) {
  plot_4 <- data |>
    filter(!is.na(amr_gene_count)) |>
    ggplot(
      aes(
        x = dataset_category,
        y = amr_gene_count,
        fill = dataset_category
      )
    ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.55,
      width = 0.55
    ) +
    geom_jitter(
      aes(colour = dataset_category),
      width = 0.14,
      alpha = 0.75,
      size = 2.3,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "Reproductive" = "#2A9D8F",
        "Bacteraemia" = "#E9C46A",
        "Unknown" = "#6C757D"
      )
    ) +
    scale_colour_manual(
      values = c(
        "Reproductive" = "#006D77",
        "Bacteraemia" = "#BC6C25",
        "Unknown" = "#495057"
      )
    ) +
    labs(
      title = "AMR burden by dataset category",
      subtitle = paste(
        "Functional AMR gene counts in reproductive",
        "and bacteraemia genomes"
      ),
      x = NULL,
      y = "AMR gene count per genome"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none"
    )

  save_plot(
    plot_4,
    "04_amr_burden_by_dataset_category.png"
  )
}

# -----------------------------
# 11. Graph 5: HLGR/gentamicin marker status
# -----------------------------
hlgr_counts <- data |>
  count(
    hlgr_gentamicin_status,
    name = "genome_count"
  ) |>
  arrange(genome_count)

plot_5 <- ggplot(
  hlgr_counts,
  aes(
    x = reorder(
      hlgr_gentamicin_status,
      genome_count
    ),
    y = genome_count,
    fill = hlgr_gentamicin_status
  )
) +
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = genome_count),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Present" = "#D00000",
      "Absent" = "#2A9D8F",
      "Unknown" = "#ADB5BD"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    breaks = scales::pretty_breaks()
  ) +
  labs(
    title = "HLGR/gentamicin marker status",
    subtitle = paste(
      "Genome counts by high-level gentamicin",
      "resistance marker status"
    ),
    x = NULL,
    y = "Number of genomes"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  plot_5,
  "05_hlgr_gentamicin_status.png"
)

# -----------------------------
# 12. Graph 6: isolation-source distribution
# -----------------------------
source_counts <- data |>
  count(
    isolation_source,
    name = "genome_count"
  ) |>
  arrange(genome_count)

plot_6 <- ggplot(
  source_counts,
  aes(
    x = reorder(
      isolation_source,
      genome_count
    ),
    y = genome_count,
    fill = isolation_source
  )
) +
  geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = genome_count),
    hjust = -0.15,
    size = 3.5
  ) +
  coord_flip(clip = "off") +
  scale_fill_viridis_d(
    option = "C",
    begin = 0.12,
    end = 0.9
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = scales::pretty_breaks()
  ) +
  labs(
    title = "Isolation-source distribution",
    subtitle = "Types of samples represented in the dataset",
    x = NULL,
    y = "Number of genomes"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  plot_6,
  "06_isolation_source_distribution.png",
  width = 10,
  height = max(
    6,
    0.35 * nrow(source_counts) + 2
  )
)

# -----------------------------
# 13. Graph 7: AMR gene count by year
# -----------------------------
if (
  sum(!is.na(data$collection_year)) > 0 &&
  sum(!is.na(data$amr_gene_count)) > 0
) {
  plot_7 <- data |>
    filter(
      !is.na(collection_year),
      !is.na(amr_gene_count)
    ) |>
    ggplot(
      aes(
        x = collection_year,
        y = amr_gene_count,
        colour = dataset_category,
        shape = dataset_category
      )
    ) +
    geom_point(
      size = 3,
      alpha = 0.78
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.85
    ) +
    scale_colour_manual(
      values = c(
        "Reproductive" = "#7209B7",
        "Bacteraemia" = "#F72585",
        "Unknown" = "#6C757D"
      )
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks()
    ) +
    labs(
      title = "AMR gene count by collection year",
      subtitle = paste(
        "Exploratory relationship between collection year",
        "and AMR burden"
      ),
      x = "Collection year",
      y = "AMR gene count per genome",
      colour = "Dataset category",
      shape = "Dataset category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )

  save_plot(
    plot_7,
    "07_amr_gene_count_by_year.png"
  )
}

# -----------------------------
# 14. Console summary
# -----------------------------
cat("\n============================================================\n")
cat("QUICK DATASET OVERVIEW COMPLETED\n")
cat("============================================================\n")
cat("Input file:\n", INPUT_FILE, "\n\n")
cat("Output folder:\n", OUTPUT_DIR, "\n\n")

cat(
  "Total genomes:",
  nrow(data),
  "\n"
)

cat(
  "Reproductive genomes:",
  sum(
    data$dataset_category == "Reproductive",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Bacteraemia genomes:",
  sum(
    data$dataset_category == "Bacteraemia",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Collection years available:",
  sum(!is.na(data$collection_year)),
  "\n"
)

cat(
  "AMR gene counts available:",
  sum(!is.na(data$amr_gene_count)),
  "\n"
)

cat(
  "Distinct AMR gene-count values:",
  dplyr::n_distinct(
    data$amr_gene_count[
      !is.na(data$amr_gene_count)
    ]
  ),
  "\n"
)

cat(
  "HLGR/gentamicin marker present:",
  sum(
    data$hlgr_gentamicin_status == "Present",
    na.rm = TRUE
  ),
  "\n"
)

if (
  length(valid_amr) > 0 &&
  dplyr::n_distinct(valid_amr) == 1
) {
  cat(
    "\nNOTE: Every genome has the same AMR gene count (",
    unique(valid_amr),
    ").\n",
    "The ranked and boxplot figures will therefore appear flat.\n",
    "Check that the metadata contains a per-genome AMR gene-count column,\n",
    "rather than a binary AMR-present field.\n",
    sep = ""
  )
}

cat(
  "\nCreated summary tables and all graphs supported by the available data.\n"
)

cat("============================================================\n")
