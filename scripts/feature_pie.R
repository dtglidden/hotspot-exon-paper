#!/usr/bin/env Rscript
## Plots a pie chart of the feature types used in the machine learning models
suppressPackageStartupMessages({
  library(RColorBrewer)
})

features <- as.data.frame(matrix(c(
  "hek_ws", "Exon",
  "hek_wu", "Exon",
  "mutation_base_change", "Mutation",
  "w5score", "Motif",
  "w3score", "Motif",
  "m5score", "Mutation",
  "m3score", "Mutation",
  "mwdif_5score", "Mutation",
  "mwdif_3score", "Mutation",
  "ss5usage", "Transcript",
  "ss3usage", "Transcript",
  "chasin_ese_density", "Exon",
  "chasin_ess_density", "Exon"
), ncol=2, byrow=T, dimnames=list(c(), c("feature", "level"))))

fLevels <- factor(features$level, c("Mutation", "Motif", "Exon", "Transcript"))
tFLevels <- table(fLevels)
labels <- sprintf("%s (%d)", names(tFLevels), tFLevels)
myColors <- brewer.pal(length(tFLevels), "Set2")
pdf(file.path("..", "plots", "feature_pie.pdf"))
pie(tFLevels, labels=labels, col=myColors, cex=2)
dev.off()
