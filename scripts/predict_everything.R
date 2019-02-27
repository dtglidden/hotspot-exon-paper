#!/usr/bin/env Rscript
## Predict the effect of all exonic mutations on splicing
## Rank exons by how many predicted mutations affect splicing
## NB: You had to remove hek_wu and hek_ws from the GBM model
## since these cannot be calculated from the genomic sequence

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  library(sqldf)
  library(gbm)
  library(VariantAnnotation)
})


gr <- readRDS(file.path("..", "data", "allMutsScored.rds"))
gbmModel <- readRDS(file.path("..", "data", "ml_model_no_hek.rds"))

test <- as.data.frame(mcols(gr))[, c(
  "mutation_base_change",
  "w5score",
  "w3score",
  "m5score",
  "m3score",
  "mwdif_5score",
  "mwdif_3score",
  "ss5usage",
  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

## 'response' type scales values to the 0-1 interval
## (i.e. the probability of affecting splicing or not)
gr$splicingProb <- predict(gbmModel$model, test, n.trees=gbmModel$n.trees, type="response")
saveRDS("gr", file.path("..", "data", "allMutsPredicted.rds"))

vcf <- VCF(gr,
           fixed=DataFrame(REF=DNAStringSet(gr$ref),
                           ALT=DNAStringSet(gr$alt)),
           info=DataFrame(SP=gr$splicingProb),
           collapsed=F)
names(geno(vcf)) <- character(0)
writeVcf(vcf, file.path("..", "data", "allMutsPredicted.vcf"))
