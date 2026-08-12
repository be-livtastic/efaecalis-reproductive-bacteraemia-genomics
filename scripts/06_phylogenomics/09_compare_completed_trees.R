#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing argument: ", flag)
  args[[i + 1]]
}
metadata_path <- value("--metadata")
output_dir <- value("--output-dir")
tree_specs <- strsplit(value("--trees"), ",", fixed = TRUE)[[1]]
if (length(tree_specs) != 3) stop("Exactly three name=treefile specifications are required")
tree_paths <- setNames(sub("^[^=]+=", "", tree_specs), sub("=.*$", "", tree_specs))
if (any(!file.exists(tree_paths))) stop("Missing tree: ", paste(tree_paths[!file.exists(tree_paths)], collapse = ", "))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
group <- setNames(meta$reproductive_bacteraemia_category, meta$assembly_accession)
trees <- lapply(tree_paths, read.tree)
shared <- Reduce(intersect, lapply(trees, function(x) x$tip.label))

descendant_tips <- function(tree, node) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  out <- integer()
  for (child in children) out <- c(out, if (child <= Ntip(tree)) child else descendant_tips(tree, child))
  out
}
support_pair <- function(label) {
  if (is.null(label) || is.na(label) || label == "") return(c(NA_real_, NA_real_))
  z <- strsplit(label, "/", fixed = TRUE)[[1]]
  if (length(z) == 1) c(NA_real_, as.numeric(z)) else as.numeric(z[1:2])
}
canonical_split <- function(a, universe) {
  b <- setdiff(universe, a)
  a <- sort(a); b <- sort(b)
  if (length(a) < length(b)) return(paste(a, collapse = ";"))
  if (length(b) < length(a)) return(paste(b, collapse = ";"))
  min(paste(a, collapse = ";"), paste(b, collapse = ";"))
}
splits_for_tree <- function(tree, tree_name) {
  rows <- list()
  for (node in seq_len(tree$Nnode) + Ntip(tree)) {
    tips <- intersect(tree$tip.label[descendant_tips(tree, node)], shared)
    other <- setdiff(shared, tips)
    if (length(tips) < 2 || length(other) < 2) next
    sp <- support_pair(tree$node.label[node - Ntip(tree)])
    key <- canonical_split(tips, shared)
    small <- strsplit(key, ";", fixed = TRUE)[[1]]
    rows[[length(rows) + 1]] <- data.frame(
      tree = tree_name, split_key = key, smaller_side_size = length(small),
      sh_alrt = sp[1], ufboot = sp[2], supported = !is.na(sp[1]) && sp[1] >= 80 && sp[2] >= 95,
      reproductive_n = sum(group[small] == "Reproductive", na.rm = TRUE),
      bacteraemia_n = sum(group[small] == "Bacteraemia", na.rm = TRUE), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
splits <- do.call(rbind, Map(splits_for_tree, trees, names(trees)))
major <- splits[splits$supported & splits$smaller_side_size >= 4, ]
presence <- table(major$split_key, factor(major$tree, levels = names(trees)))
major$shared_by_n_trees <- rowSums(presence > 0)[match(major$split_key, rownames(presence))]
major <- major[order(-major$shared_by_n_trees, -major$smaller_side_size, major$tree), ]
write.table(major, file.path(output_dir, "supported_major_clades.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

placement_rows <- list()
branch_rows <- list()
for (nm in names(trees)) {
  tr <- trees[[nm]]
  d <- cophenetic.phylo(tr)
  terminal <- setNames(tr$edge.length[match(seq_len(Ntip(tr)), tr$edge[, 2])], tr$tip.label)
  med <- median(terminal); madv <- mad(terminal, constant = 1)
  for (tip in tr$tip.label) {
    branch_rows[[length(branch_rows) + 1]] <- data.frame(
      tree = nm, assembly_accession = tip, group = group[tip], terminal_branch = terminal[tip],
      terminal_over_median = terminal[tip] / med,
      robust_z = if (madv > 0) (terminal[tip] - med) / madv else NA_real_, stringsAsFactors = FALSE)
  }
  supported <- splits[splits$tree == nm & splits$supported, ]
  full_supported_sides <- list()
  for (node in seq_len(tr$Nnode) + Ntip(tr)) {
    sp <- support_pair(tr$node.label[node - Ntip(tr)])
    if (is.na(sp[1]) || sp[1] < 80 || sp[2] < 95) next
    side <- tr$tip.label[descendant_tips(tr, node)]
    if (length(side) >= 2) full_supported_sides[[length(full_supported_sides) + 1]] <- side
    complement <- setdiff(tr$tip.label, side)
    if (length(complement) >= 2) full_supported_sides[[length(full_supported_sides) + 1]] <- complement
  }
  reps <- intersect(names(group)[group == "Reproductive"], tr$tip.label)
  for (tip in reps) {
    candidates <- full_supported_sides[vapply(full_supported_sides, function(x) tip %in% x, logical(1))]
    local <- if (length(candidates)) candidates[[which.min(vapply(candidates, length, integer(1)))]] else character()
    distances <- d[tip, setdiff(tr$tip.label, tip)]
    nearest <- names(which.min(distances))
    placement_rows[[length(placement_rows) + 1]] <- data.frame(
      tree = nm, assembly_accession = tip, nearest_tip = nearest,
      nearest_tip_group = group[nearest], nearest_distance = min(distances),
      smallest_supported_clade_size = if (length(local)) length(local) else NA_integer_,
      reproductive_in_supported_clade = if (length(local)) sum(group[local] == "Reproductive") else NA_integer_,
      bacteraemia_in_supported_clade = if (length(local)) sum(group[local] == "Bacteraemia") else NA_integer_,
      mean_distance_to_other_reproductive = mean(d[tip, setdiff(reps, tip)]),
      mean_distance_to_bacteraemia = mean(d[tip, intersect(names(group)[group == "Bacteraemia"], tr$tip.label)]),
      stringsAsFactors = FALSE)
  }
}
placements <- do.call(rbind, placement_rows)
branches <- do.call(rbind, branch_rows)
write.table(placements, file.path(output_dir, "reproductive_isolate_placements.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(branches, file.path(output_dir, "terminal_branch_lengths.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

tree_summary <- do.call(rbind, lapply(names(trees), function(nm) {
  tr <- trees[[nm]]; b <- branches[branches$tree == nm, ]
  data.frame(tree = nm, taxa = Ntip(tr), total_tree_length = sum(tr$edge.length),
    zero_terminal_branches = sum(b$terminal_branch == 0),
    terminal_branch_95th_percentile = unname(quantile(b$terminal_branch, 0.95)),
    maximum_terminal_branch = max(b$terminal_branch),
    supported_internal_splits = sum(splits$tree == nm & splits$supported),
    supported_major_splits = sum(major$tree == nm), stringsAsFactors = FALSE)
}))
write.table(tree_summary, file.path(output_dir, "tree_comparison_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("Compared", length(trees), "trees; shared taxa:", length(shared), "\n")
print(tree_summary, row.names = FALSE)
cat("Supported major split keys present in all three trees:", sum(rowSums(presence > 0) == 3), "\n")
