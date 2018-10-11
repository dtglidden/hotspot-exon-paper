#!/usr/bin/env Rscript
## Calculate Maxent scores for all exons

source(file.path("..", "lib", "feature_gen.R"))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exs <- exons(txdb)
ss5 <- Score5SS(exs)
ss3 <- Score3SS(exs)
saveRDS(list(ss5=ss5, ss3=ss3), file.path("..", "data", "wt_maxent.rds"))
