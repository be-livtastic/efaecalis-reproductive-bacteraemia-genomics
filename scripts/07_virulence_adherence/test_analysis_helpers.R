#!/usr/bin/env Rscript

# Load helper functions without creating any analysis output.
Sys.setenv(EFAECALIS_SKIP_MAIN = "1")
script_arg <- commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
source(file.path(dirname(script_path), "analyse_virulence_adherence.R"))

# Check conservative Ebp propagation and accepted-only completeness.
stopifnot(composite_ebp_status("accepted_present", "accepted_present", "accepted_present") == "accepted_present")
stopifnot(composite_ebp_status("accepted_present", "ambiguous_multiple_hit", "accepted_present") == "ambiguous_multiple_hit")
stopifnot(composite_ebp_status("accepted_present", "review_required", "accepted_present") == "review_required")
stopifnot(composite_ebp_status("accepted_present", "flagged_partial", "accepted_present") == "flagged_partial")
stopifnot(composite_ebp_status("accepted_present", "not_detected", "accepted_present") == "not_detected")

# Check that BH correction uses the base-R implementation and preserves ordering.
raw <- c(0.001, 0.02, 0.3)
adjusted <- p.adjust(raw, method = "BH")
stopifnot(isTRUE(all.equal(adjusted, c(0.003, 0.03, 0.3))))
cat("Virulence R helper tests passed\n")
