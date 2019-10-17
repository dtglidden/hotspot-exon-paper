#!/usr/bin/env Rscript
## Plots a pie chart of the feature types used in the machine learning models
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(stringi)
})

features <- read.csv(file.path("..", "data", "feature_table.csv"), stringsAsFactors=F)

fLevels <- factor(stri_trans_totitle(features$level), c("Mutation", "Motif", "Exon", "Transcript"))
tFLevels <- table(fLevels)
labels <- sprintf("%s (%d)", names(tFLevels), tFLevels)
myColors <- brewer.pal(length(tFLevels), "Set2")
pdf(file.path("..", "plots", "feature_pie.pdf"), width=9)
pie(tFLevels, labels=labels, col=myColors, cex=2)
dev.off()
