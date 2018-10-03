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
  introns <- gr[unique(subjectHits(findOverlaps(exs, gr, maxgap=1)))]
  ## Get introns that don't overlap any exons. Makes sure the splice site that was not matched in the previous command
  ## is not cryptic (going into the exon). Using annotated splice sites should ensure not cryptic splicing into the intron
  introns <- introns[!introns %over% exs]
  return(introns)
}
