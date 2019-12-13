#!/usr/bin/env Rscript
## Predict the effect of all exonic mutations on splicing
## Rank exons by how many predicted mutations affect splicing
## NB: You had to remove hek_wu and hek_ws from the GBM model
## since these cannot be calculated from the genomic sequence

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
  library(sqldf)
  library(dplyr)
  library(gbm)
  library(VariantAnnotation)
})


gr <- readRDS(file.path("..", "data", "allMutsScored.rds"))
gbmModel <- readRDS(file.path("..", "data", "ml_model_no_hek.rds"))

featureDf <- read.csv(file.path("..", "data", "feature_table.csv"), stringsAsFactors=F) %>%
  filter(!feature %in% c("hek_wu", "hek_ws"))
test <- as.data.frame(mcols(gr))[, featureDf$feature]

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
