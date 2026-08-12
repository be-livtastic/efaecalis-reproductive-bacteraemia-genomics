#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) stop("Missing ", flag); args[[i+1]] }
scores <- read.delim(value("--scores"), check.names = FALSE)
variance <- read.delim(value("--variance"), check.names = FALSE)
out <- value("--output-dir"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
cols <- c(Reproductive = "#0072B2", Bacteraemia = "#D55E00")
plot_pair <- function(x, y, stem) {
  png(file.path(out, paste0(stem, ".png")), width=1800, height=1400, res=180)
  par(mar=c(5,5,1,1)); g <- scores$reproductive_bacteraemia_category
  plot(scores[[x]], scores[[y]], pch=19, col=cols[g], xlab=sprintf("%s (%.2f%%)",x,100*variance$explained_variance_proportion[variance$component==x]), ylab=sprintf("%s (%.2f%%)",y,100*variance$explained_variance_proportion[variance$component==y]))
  legend("topright", names(cols), col=cols, pch=19, bty="n")
  reps <- which(g == "Reproductive"); text(scores[[x]][reps], scores[[y]][reps], labels=scores$assembly_accession[reps], pos=3, cex=.48)
  dev.off()
  pdf(file.path(out, paste0(stem, ".pdf")), width=10, height=8)
  plot(scores[[x]], scores[[y]], pch=19, col=cols[g], xlab=sprintf("%s (%.2f%%)",x,100*variance$explained_variance_proportion[variance$component==x]), ylab=sprintf("%s (%.2f%%)",y,100*variance$explained_variance_proportion[variance$component==y]))
  legend("topright", names(cols), col=cols, pch=19, bty="n"); text(scores[[x]][reps], scores[[y]][reps], labels=scores$assembly_accession[reps], pos=3, cex=.48); dev.off()
}
plot_pair("PC1", "PC2", "pca_pc1_pc2")
plot_pair("PC1", "PC3", "pca_pc1_pc3")
message("PCA plots written for ", nrow(scores), " genomes")
