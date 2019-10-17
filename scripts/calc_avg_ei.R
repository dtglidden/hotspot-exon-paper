#!/usr/bin/env Rscript
## Calculate the average EI score in all exons
## Hexamers that don't have 0 score are not included ('none' in the hexamer table)

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
})

exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
system.time( exs$ei <- AvgEI(exs) )
saveRDS(exs, file.path("..", "data", "exonic_avg_ei.rds"))
