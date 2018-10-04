## This file has functions for calculating genomic features (i.e. the don't depend on a supplied mutation)
library(GenomicRanges, quietly=T)
library(BSgenome.Hsapiens.UCSC.hg19, quietly=T)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rPython, quietly=T)

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

## Calculate the 5'SS or 3'SS usage for all splice sites in an RNA-seq study aligned with STAR
## file: The file holding all splice junctions (SJ.out.tab)
## From the RNAStar manual (column descriptions):
## column 1: chromosome
## column 2: first base of the intron (1-based)
## column 3: last base of the intron (1-based)
## column 4: strand (0: undefined, 1: +, 2: -)
## column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
## AT/AC, 6: GT/AT
## column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
## column 7: number of uniquely mapping reads crossing the junction
## column 8: number of multi-mapping reads crossing the junction
## column 9: maximum spliced alignment overhang
## Returns: A 2-element list, where the first element has all exons with associated 5'SS usage,
## and the second element has all exons with associated 3'SS usage
## We must split up the exons this way because we might not get enough sequencing coverage to
## calculate the SS usage for both sites in all exons
SSUsage <- function(file, exs) {
  python.load("ss_usage.py")
  ## TODO Could get the nrows to pass to read.csv
  tempfile <- python.call("get_usage", file)
  df <- read.csv(tempfile, colClasses=c("annotated" = "logical", "intron_location" = "character"))
  file.remove(tempfile)
  levels(df$strand) <- c("-", "+", "*") # The "*" strand is originally loaded as "0"
  gr <- GRanges(seqnames=df$chrom,
                ranges=IRanges(start=df$intron_start, end=df$intron_stop),
                strand=df$strand,
                intron_motif=df$intron_motif,
                annotated=df$annotated,
                psi5=df$psi5,
                psi3=df$psi3,
                seqinfo=seqinfo(exs))
  ## Ensuring that we use annotated splice sites prevents the case when we have an intron that goes between splice sites
  ## of exons of two overlapping genes
  gr <- gr[mcols(gr)$annotated]
  gr <- gr[strand(gr) != "*"]

  ## maxgap=1 allows introns adjacent to exons to be considered 'overlapping'
  ## Essentially gets junctions that go to at least one canonical splice site
  op <- findOverlapPairs(exs, gr, maxgap=1)
  ss5 <- op[abs(end(first(op)) - start(second(op))) == 1]
  ss3 <- op[abs(start(first(op)) - end(second(op))) == 1]

  ## Determine which exon IDs are duplicated, so we can check those for the best psi5/psi3 values
  ## This should speed up the endoapply step by only looking at splice sites with multiple options
  ss5Dupes <- duplicated(first(ss5)$exon_id) | duplicated(first(ss5)$exon_id, fromLast=T)
  ss3Dupes <- duplicated(first(ss3)$exon_id) | duplicated(first(ss3)$exon_id, fromLast=T)

  ## Split the splice site exon/intron pairs into those that have duplicate exons and those that don't
  ss5ByDupes <- split(ss5, ss5Dupes)
  ss3ByDupes <- split(ss3, ss3Dupes)

  ## Split the duplicated exons by their exon IDs
  ss5ByDupesAndExID <- split(second(ss5ByDupes$"TRUE"), first(ss5ByDupes$"TRUE")$exon_id)
  ss3ByDupesAndExID <- split(second(ss3ByDupes$"TRUE"), first(ss3ByDupes$"TRUE")$exon_id)

  ## SLOW: Select the exons from each set of duplicates that has the highest psi value for its relevant splice site
  ## Then merge the GRangesList back into a GRanges object, because we don't need separation by exon ID anymore
  goodSS5Introns <- unlist(endoapply(ss5ByDupesAndExID, function(gr) gr[which(gr$psi5 == max(gr$psi5))[[1]]]))
  goodSS3Introns <- unlist(endoapply(ss3ByDupesAndExID, function(gr) gr[which(gr$psi3 == max(gr$psi3))[[1]]]))

  ## Get the list of exons after filtering out duplicates
  ## The ordering of exons is the same as the introns we previously selected because the 'split' function
  ## does not scramble the exon IDs
  goodSS5Exons <- first(ss5ByDupes$"TRUE")[!duplicated(first(ss5ByDupes$"TRUE")$exon_id)]
  goodSS3Exons <- first(ss3ByDupes$"TRUE")[!duplicated(first(ss3ByDupes$"TRUE")$exon_id)]

  ## Attach the SS usage values to the exons
  ## This is defined as the psi value opposite to the splice site in question
  ## For example, for SS5 usage, this is the psi3 value for the intron with the highest psi5 value
  goodSS5Exons$usage <- goodSS5Introns$psi3
  goodSS3Exons$usage <- goodSS3Introns$psi5

  return(list(ss5=goodSS5Exons, ss3=goodSS3Exons))
}
