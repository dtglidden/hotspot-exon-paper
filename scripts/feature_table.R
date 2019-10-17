#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gridExtra)
})

## Make these variables to reduce spelling errors
mutation <- "mutation"
motif <- "motif"
exon <- "exon"
transcript <- "transcript"

featureDf <- as.data.frame(matrix(c(
  "mutation_base_change", mutation, "The actual SNV (e.g G->T)",
  "isWtPaired", mutation, "RNA-RNA hybridization state at the position of the mutation for the WT base",
  "isMutPaired", mutation, "RNA-RNA hybridization state at the position of the mutation for the Mut base",
  "mwdif_3score", motif, "Mut 3'SS score - WT 3'SS score",
  "mwdif_5score", motif, "Mut 5'SS score - WT 5'SS score",
  "avgDeltaEIScores", motif, "Average change in hexamer EI for all hexamers affected by the mutation",
  "avgDeltaA3SSScores", motif, "Average change in Shendure A3SS score for all hexamers affected by the mutation",
  "avgDeltaA5SSScores", motif, "Average change in Shendure A5SS score for all hexamers affected by the mutation",
  "hek_ws", exon, "MaPSy WT in vivo spliced counts",
  "hek_wu", exon, "MaPSy WT in vivo unspliced counts",
  "w3score", exon, "WT 3'SS Maxent score",
  "w5score", exon, "WT 5'SS Maxent score",
  "ss3usage", exon, "3'SS usage in HEK 293T cells",
  "ss5usage", exon, "5'SS usage in HEK 293T cells",
  "ei", exon, "Average EI score for all hexamers in the exon",
  "nExons", transcript, "Number of exons in the canonical transcript",
  "geneLength", transcript, "Full length of the gene"),
  ncol=3, byrow=T, dimnames=list(c(), c("feature", "level", "description"))), stringsAsFactors=F)
write.csv(featureDf, file.path("..", "data", "feature_table.csv"), quote=F, row.names=F)

pdf(file.path("..", "plots", "feature_table.pdf"), width=9.5, height=5)
grid.table(featureDf)
dev.off()
