## This file has functions for calculating genomic features (i.e. the don't depend on a supplied mutation)
library(GenomicRanges, quietly=T)
library(BSgenome.Hsapiens.UCSC.hg19, quietly=T)

## Calls a MAXENT Perl script on a list of sequences
## seq: A character vector of splice site sequences
## script: The appropriate perl script to compute the 5'SS or 5'SS scores (either "score5.pl" or "score3.pl)
## Returns: A dataframe with the sequences and their associated splice site scores
MaxentPerl <- function(seqs, script="score5.pl") {
  wd <- getwd()
  setwd(file.path("..", "maxent"))
  scores <- system2("perl", script, stdout=T, input=seqs)
  df <- data.frame(do.call(rbind, strsplit(scores, "\t", fixed=T)), stringsAsFactors=F)
  colnames(df) <- c("seq", "score")
  df$score <- as.numeric(df$score)
  setwd(wd)
  return(df)
}

## Calculate the MAXENT score for the 5' splice sites of exons
## exons: A GRanges object that contains the genomic coordinates for the exons to be scored
## Returns: A dataframe with the sequences and their associated splice site scores
Score5SS <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  # The 5'SS is the last 3 nucleotides of the exon + the next 6 intronic nucleotides
  # Add the 6 downstream intronic bases to the coordinates
  sites <- flank(exons, 6, start=F)
  # Then make sure the total width of the range is the length of the 5'SS (starting from the last nucleotide in the site)
  sites <- resize(sites, 9, fix="end")

  seqs <- getSeq(genome, sites, as.character=T)

  return(MaxentPerl(seqs))
}

## Calculate the MAXENT score for the 3' splice sites of exons
## exons: A GRanges object that contains the genomic coordinates for the exons to be scored
## Returns: A dataframe with the sequences and their associated splice site scores
Score3SS <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  # The 3'SS is 20 nucleotides of the upstream intron + the first 3 nucleotides of the exon
  # Add the 20 upstream intronic bases to the coordinates
  sites <- flank(exons, 20)
  # Then make sure the total width of the range is the length of the 3'SS (starting from the first nucleotide in the site)
  sites <- resize(sites, 23)

  seqs <- getSeq(genome, sites, as.character=T)

  return(MaxentPerl(seqs, script="score3.pl"))
}
