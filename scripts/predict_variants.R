#!/usr/bin/env Rscript
## Provides functions to predict exon skipping in a list of variants

## You got the liftOver chain file from:
## rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz .
## You have to gunzip it

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(rtracklayer)
  library(gbm)
  source(file.path("..", "lib", "feature_gen.R"))
  source(file.path("..", "lib", "feature_query.R"))
  source(file.path("..", "lib", "granges_util.R"))
})

## Loads a VCF file into a GRanges object
## vcf: a VCF object
## gnome: a string representing the genome of the VCF
ImportVcf <- function(vcf, gnome="hg19", unknownStrand=F) {
  si <- Seqinfo(genome=gnome)
  gr <- rowRanges(vcf)
  ## Explicitly setting the seqlevels first can avoid weird errors
  seqlevels(gr) <- seqlevels(si)
  seqinfo(gr) <- si
  #if (!unknownStrand) {
  #  strand(gr) <- "+"
  #}
  return(gr)
}

## Prepares imported data for prediction
## gr: a GRanges object of imported variants
## exs: a GRanges object of exons with calculated features
PrepareData <- function(gr, exs) {
  gn <- genome(gr)[[1]]
  if (gn == "hg38") {
    chain <- import.chain("hg38ToHg19.over.chain")
    gr <- unlist(liftOver(gr, chain))
  }
#  gr <- AssociateVars(gr)
  df <- mergeByOverlaps(gr, exs)
  testExs <- df$exs
  testExs$ref <- as.character(df$REF)
  testExs$alt <- as.character(unlist(df$ALT))
  testExs$pos <- start(df$gr) - start(df$exs) + 1
  return(CalcMutFeatures(testExs))
}

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

saveRDS(gr, file.path("..", "data", "hek293_vars_predicted.rds"))

vcf <- VCF(gr,
           fixed=DataFrame(REF=DNAStringSet(gr$ref),
                           ALT=DNAStringSet(gr$alt)),
           info=DataFrame(SP=gr$splicingProb),
           collapsed=F)
names(geno(vcf)) <- character(0)
writeVcf(vcf, file.path("..", "data", "hek293VarsPredicted.vcf"))
