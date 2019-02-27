#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
})

nCores <- detectCores()
#exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

bases <- c("A", "T", "C", "G")
## Returns a list of vectors where the list index
## corresponds to the nucleotide position relative to 3'SS
## and the vector holds the 3 possible mutations
AllMutsInExon <- function(ex, containsFullSpliceSites=F,
                          genome=BSgenome.Hsapiens.UCSC.hg19) {
  ## Get the WT sequence of the whole exon
  if (containsFullSpliceSites) {
    exStart <- start(ex) + 20
    exEnd <- end(ex) - 6
    exLoci <- exStart:exEnd
    exLen <- exEnd-exStart + 1
    positions <- 1:exLen
    seq <- substr(getSeq(genome, ex, as.character=T), 21, exLen + 20)
  } else {
    exStart <- start(ex)
    exEnd <- end(ex)
    exLoci <- exStart:exEnd
    exLen <- exEnd-exStart + 1
    positions <- 1:exLen
    seq <- getSeq(genome, ex, as.character=T)
  }
  return(unlist(GRangesList(mclapply(positions, function(i) {
    gr <- ex
    gr$ref <- substr(seq, i, i)
    gr$pos <- i
    mutGr <- rep(gr, 3)
    mutGr$alt <- bases[!bases %in% gr$ref]
    mutGr
  }))))
}

## Save the mutations for future reference
saveRDS(unlist(GRangesList(mclapply(exsWithSSUsage, AllMutsInExon))),
        file.path("..", "data", "allMuts.rds"))


## Calculate mutation features
allMuts <- readRDS(file.path("..", "data", "allMuts.rds"))
allMuts <- CalcMutFeatures(allMuts)
saveRDS(allMuts, file.path("..", "data", "allMutsCalced.rds"))
