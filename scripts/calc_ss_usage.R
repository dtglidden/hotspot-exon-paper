#!/usr/bin/env Rscript
## Calculates splice site usage for each exon based on RNA-seq data
## NB: Because sequencing experiments vary in their coverage, not all exons
## may have their splice sites scored, and some might only have one splice site scored

source(file.path("..", "lib", "feature_gen.R"))
seqFile <- file.path(Sys.getenv("HOME"), "scratch", "amiloride", "star", "ctrl", "SJ.out.tab")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exs <- exons(txdb)
usage <- SSUsage2(seqFile, exs)
saveRDS(usage, file.path("..", "data", "ss_usage.rds"))
