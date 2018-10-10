#!/usr/bin/env Rscript
## Calculate the Chasin ESE and ESS densities in all exons

source(file.path("..", "lib", "feature_gen.R"))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exs <- exons(txdb)
chasin_ese_density <- ChasinCisDen(exs)
chasin_ess_density <- ChasinCisDen(exs, getESEs=F)
exs$chasin_ese_density <- chasin_ese_density
exs$chasin_ess_density <- chasin_ess_density
saveRDS(exs, file.path("..", "data", "chasin_cis_element_densities.rds"))
