#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape); library(treeio); library(ggtree); library(ggplot2)
  library(dplyr); library(readr); library(tidyr); library(phangorn)
})

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing argument: ", flag)
  args[[i + 1]]
}
tree_path <- value("--tree")
metadata_path <- value("--metadata")
output_dir <- value("--output-dir")
tables_dir <- value("--tables-dir")
expected_tips <- as.integer(value("--expected-tips"))
if (!file.exists(tree_path) || !file.exists(metadata_path)) stop("Tree or metadata input is missing")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "captions"), recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_path)
metadata <- read_csv(metadata_path, show_col_types = FALSE)
required <- c("assembly_accession", "reproductive_bacteraemia_category")
if (!all(required %in% names(metadata))) stop("Metadata lacks required columns: ", paste(setdiff(required, names(metadata)), collapse=", "))
if (anyDuplicated(metadata$assembly_accession)) stop("Duplicate metadata assembly accessions")
tree_metadata <- metadata %>% filter(assembly_accession %in% tree$tip.label)
join_qc <- full_join(tibble(assembly_accession = tree$tip.label, in_tree = TRUE),
                    tree_metadata %>% mutate(in_metadata = TRUE), by = "assembly_accession") %>%
  mutate(in_tree = replace_na(in_tree, FALSE), in_metadata = replace_na(in_metadata, FALSE),
         qc_status = if_else(in_tree & in_metadata, "MATCH", "UNMATCHED"))
write_tsv(join_qc, file.path(tables_dir, "tree_metadata_join_qc.tsv"))
if (length(tree$tip.label) != expected_tips ||
    sum(join_qc$qc_status == "MATCH") != expected_tips || nrow(join_qc) != expected_tips)
  stop("Tree/metadata join does not match expected tips: ", expected_tips)

plot_data <- tree_metadata %>% rename(label = assembly_accession)
rooted <- midpoint(tree)
base_plot <- function(x, layout = "rectangular") {
  ggtree(x, layout = layout) %<+% plot_data +
    geom_tippoint(aes(colour = reproductive_bacteraemia_category), size = 2) +
    scale_colour_manual(values = c(Reproductive = "#0072B2", Bacteraemia = "#D55E00"),
                        na.value = "grey50", name = "Dataset group") +
    theme_tree2() + theme(legend.position = "right")
}
save_plot <- function(plot, stem, width = 11, height = 9) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width=width, height=height, dpi=300)
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width=width, height=height)
}
save_plot(base_plot(tree), "combined_rectangular_unrooted")
save_plot(base_plot(rooted), "combined_rectangular_midpoint_rooted")
save_plot(base_plot(tree, "fan"), "combined_fan", 10, 10)
for (group in c("Reproductive", "Bacteraemia")) {
  keep <- tree_metadata %>% filter(reproductive_bacteraemia_category == group) %>% pull(assembly_accession)
  pruned <- keep.tip(rooted, keep)
  save_plot(base_plot(pruned), paste0(tolower(group), "_only_pruned"), 10, ifelse(group == "Reproductive", 6, 10))
}
distances <- tibble(label = tree$tip.label, root_to_tip = node.depth.edgelength(rooted)[seq_along(tree$tip.label)])
diagnostic <- ggplot(distances, aes(x = reorder(label, root_to_tip), y = root_to_tip)) +
  geom_col(fill = "#666666") + coord_flip() + labs(x = "Assembly accession", y = "Midpoint-rooted display distance") +
  theme_minimal(base_size = 9)
save_plot(diagnostic, "branch_length_diagnostic", 9, 12)
write_tsv(join_qc, file.path(tables_dir, "tip_to_metadata_qc.tsv"))
caption <- paste(
  "Nine-locus concatenated housekeeping-gene maximum-likelihood phylogeny of 72 genomes.",
  "Colours distinguish dataset categories. Midpoint rooting is used only for display;",
  "the inferred unrooted topology is retained. Genomic proximity does not establish",
  "pathogenicity, transmission, or epidemiological linkage."
)
writeLines(caption, file.path(output_dir, "captions", "phylogeny_9_locus_caption.txt"))
message("Figures generated; tree tips and metadata matched exactly ", expected_tips, ":", expected_tips, ".")
