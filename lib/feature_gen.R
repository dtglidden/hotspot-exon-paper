## This file has functions for calculating genomic features (i.e. the don't depend on a supplied mutation)
library(GenomicRanges, quietly=T)
library(BSgenome.Hsapiens.UCSC.hg19, quietly=T)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rPython, quietly=T)
library(sqldf, quietly=T)
source(file.path("..", "lib", "feature_query.R"))

## Calls a MAXENT Perl script on a list of sequences
## seq: A character vector of splice site sequences
## script: The appropriate perl script to compute the 5'SS or 5'SS scores (either "score5.pl" or "score3.pl)
## Returns: A dataframe with the sequences and their associated splice site scores
MaxentPerl <- function(seqs, script="score5.pl") {
  wd <- getwd()
  setwd(file.path("..", "maxent"))
  scores <- system2("perl", script, stdout=T, input=seqs$seq)
  df <- data.frame(do.call(rbind, strsplit(scores, "\t", fixed=T)), stringsAsFactors=F)
  colnames(df) <- c("seq", "score")
  df$exon_id <- seqs$exon_id
  df$score <- as.numeric(df$score)
  setwd(wd)
  return(df)
}

## Get the 5'SS sequences for a list of exons
## exons: A GRanges object of exons
GetSS5Seq <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  # The 5'SS is the last 3 nucleotides of the exon + the next 6 intronic nucleotides
  # Add the 6 downstream intronic bases to the coordinates
  sites <- flank(exons, 6, start=F)

  # Then make sure the total width of the range is the length of the 5'SS (starting from the last nucleotide in the site)
  sites <- resize(sites, 9, fix="end")

  # Necessary when the resizing functions above create GRanges that are outside of their chromosomes
  sites <- sites[end(sites) <= seqlengths(exs)[as.character(seqnames(sites))] & start(sites) > 0]

  return(list(seq=getSeq(genome, sites, as.character=T), exon_id=sites$exon_id))
}

## Get the 3'SS sequences for a list of exons
## exons: A GRanges object of exons
GetSS3Seq <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  # The 3'SS is 20 nucleotides of the upstream intron + the first 3 nucleotides of the exon
  # Add the 20 upstream intronic bases to the coordinates
  sites <- flank(exons, 20)

  # Then make sure the total width of the range is the length of the 3'SS (starting from the first nucleotide in the site)
  sites <- resize(sites, 23)

  # Necessary when the resizing functions above create GRanges that are outside of their chromosomes
  sites <- sites[end(sites) <= seqlengths(exs)[as.character(seqnames(sites))] & start(sites) > 0]

  return(list(seq=getSeq(genome, sites, as.character=T), exon_id=sites$exon_id))
}

## Calculate the MAXENT score for the 5' splice sites of exons
## exons: A GRanges object that contains the genomic coordinates for the exons to be scored
## Returns: A dataframe with the sequences and their associated splice site scores
Score5SS <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  return(MaxentPerl(GetSS5Seq(exons, genome)))
}

## Calculate the MAXENT score for the 3' splice sites of exons
## exons: A GRanges object that contains the genomic coordinates for the exons to be scored
## Returns: A dataframe with the sequences and their associated splice site scores
Score3SS <- function(exons, genome=BSgenome.Hsapiens.UCSC.hg19) {
  return(MaxentPerl(GetSS3Seq(exons, genome), script="score3.pl"))
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
  python.load(file.path("..", "lib", "ss_usage.py"))
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

## Gets SS usage data based on another method
## files: A character vector of paths to splice junctions files
## exs: A GRanges object of exons to have their splice site usages added
## cellLine: A string representing the cell line associated with the splicing data
SSUsage2 <- function(files, exs, cellLine="hek293", gtf=file.path("..", "data", "hg19.gtf")) {
  script <- file.path("..", "scripts", "splice_site_usage.py")
  fileSpec <- sprintf("%s:%s", cellLine, paste(files, collapse=","))
  outDir <- tempdir()
  print(outDir)
  system2("python", args=c(script, gtf, fileSpec, outDir))
  outFiles <- file.path(outDir, dir(outDir))
  ss5Usage <- read.delim(outFiles[[2]])
  ss3Usage <- read.delim(outFiles[[1]])
  ss5UsageGr <- GRanges(seqnames=ss5Usage$Chrom,
                        ranges=IRanges(start=ss5Usage$Site, end=ss5Usage$Site),
                        strand=ss5Usage$Strand,
                        seqinfo=Seqinfo(genome="hg19"),
                        usage=ss5Usage$Usage)
  ss3UsageGr <- GRanges(seqnames=ss3Usage$Chrom,
                        ranges=IRanges(start=ss3Usage$Site, end=ss3Usage$Site),
                        strand=ss3Usage$Strand,
                        seqinfo=Seqinfo(genome="hg19"),
                        usage=ss3Usage$Usage)
  file.remove(outFiles)

  ## Associate exons with splice site data
  ss5Merge <- mergeByOverlaps(ss5UsageGr, exs, type="end")
  ss3Merge <- mergeByOverlaps(ss3UsageGr, exs, type="start")

  ## Pull out GRanges and add usage data to it
  ss5UsageGr <- ss5Merge$exs
  ss3UsageGr <- ss3Merge$exs
  ss5UsageGr$usage <- ss5Merge$usage
  ss3UsageGr$usage <- ss3Merge$usage

  return(list(ss5=ss5UsageGr, ss3=ss3UsageGr))
}

## Tile through a string (window size n) and increment if a substring occurs in an environment (hashtable)
countSubstringsOfN <- function(x, n, e) {
  nc <- nchar(x)
  return(sum(substring(x, 1:(nc - n + 1), n:nc) %in% e))
}

## Calculate the density of Chasin exonic cis-elements in exons
## exons: A GRanges object of exons
## genome: A BSgenome object that holds sequence info
## getESEs: A logical determining if we should calculate the ESE density (vs. the ESS density)
ChasinCisDen <- function(exons, getESEs=T, genome=BSgenome.Hsapiens.UCSC.hg19) {
  elems <- read.table(file.path("..", "data", ifelse(getESEs, "chasin_enhancers.txt", "chasin_silencers.txt")),
                      stringsAsFactors=F)$V1
  exLn <- width(exons)

  countSubsV <- Vectorize(countSubstringsOfN, "x", USE.NAMES=F)
  total <- countSubsV(getSeq(genome, exons, as.character=T), 6, elems)
  return(total / exLn)
}

## Calculates MaPSy-style log2 splicing efficiency
## spliced: Numeric vector of read counts for spliced species
## unspliced: Numeric vector of read counts for unspliced species
MapsySplicingEfficiency <- function(spliced, unspliced) {
  splicedPseudoCounts <- spliced + 1
  splicedSum <- sum(splicedPseudoCounts)
  splicedEfficiency <- splicedPseudoCounts / splicedSum

  unsplicedPseudoCounts <- unspliced + 1
  unsplicedSum <- sum(unsplicedPseudoCounts)
  unsplicedEfficiency <- unsplicedPseudoCounts / unsplicedSum

  return(log2(splicedEfficiency/unsplicedEfficiency))
}

## Calculates allelic skew (a 1.5-fold change defines a mutation as an ESM)
AllelicSkew <- function(ws, wu, ms, mu) {
  wuPseudo <- wu + 1
  muPseudo <- mu + 1
  wsPseudo <- ws + 1
  msPseudo <- ms + 1
  return(log2((msPseudo/muPseudo) / (wsPseudo/wuPseudo)))
}


#######################
## Mutation Features ##
#######################

## Load a Mutation file (parses text supplied to web form)
## Output dataframe columns:
## chr: chromosome name
## pos: 1-based mutation position
## strand: DNA strand ("+" or "-")
## ref: reference allele
## alt: alternative allele
ParseMutText <- function(mutFile) {
  return(read.table(mutFile, sep=" ",
                    col.names=c("chr", "pos", "strand", "ref", "alt"),
                    colClasses=c("ref"="character", "alt"="character")))
}

## Get the exons associated with the mutation data (includes splice sites in the ranges)
## mutations: A dataframe of mutation info (should be parsed from ParseMutText of ParseMutVcf)
## exs: A GRanges object of all exons in a genome
## Returns a GRanges object with the ranges of the exons that overlap the mutations
## Passes along the mutation info in the metadata columns
Muts2Exons <- function(mutations, exs=exons(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
  ## Cover the intronic portion of the 3'SS
  exsF <- flank(exs, 20)

  ## Resize to cover through the intronic portion of the 5'SS
  exs <- resize(exsF, width(exs) + 26)

  mutsGr <- GRanges(mutations$chr, IRanges(mutations$pos, width=1), mutations$strand, seqinfo=seqinfo(exs),
                    ref=mutations$ref, alt=mutations$alt)
  mutsAndExs <- mergeByOverlaps(mutsGr, exs)
  mutsAndExs <- mutsAndExs[!duplicated(mutsAndExs$mutsGr), ]
  goodExons <- mutsAndExs$exs
  goodExons$pos <- start(mutsAndExs$mutsGr)
  goodExons$ref <- mutsAndExs$ref
  goodExons$alt <- mutsAndExs$alt
  return(goodExons)
}

## Get sequences with a mutation applied (includes complete 5' and 3'SS)
## mutExons A GRanges object with exon ranges, and metadata holding mutation info
## relative to the start of the sequence (1-based)
GetMutSeq <- function(mutExons, genome=BSgenome.Hsapiens.UCSC.hg19, hasSpliceSites=T, relativePos=F) {
  if (!hasSpliceSites) {
    start(mutExons) <- ifelse(strand(mutExons) == "+", start(mutExons) - 20, start(mutExons) - 6)
    end(mutExons) <- ifelse(strand(mutExons) == "+", end(mutExons) + 6, end(mutExons) + 20)
    mutExons$pos <- mutExons$pos + 20
  }
  wtSeqs <- getSeq(genome, mutExons, as.character=T)
  exStart <- ifelse(strand(mutExons) == "-", end(mutExons), start(mutExons))
  if (relativePos) {
    mutLoc <- mutExons$pos
  } else {
    mutLoc <- ifelse(strand(mutExons) == "-", exStart - mutExons$pos, mutExons$pos - exStart)
  }
  mutSeqs <- paste0(substring(wtSeqs, 1, mutLoc),
                    mutExons$alt,
                    substring(wtSeqs, mutLoc + 2, nchar(wtSeqs)))
  return(mutSeqs)
}

## Gets 5'SS Maxent scores for exons with associated mutations
## mutations: A dataframe containing mutation info (chr, pos, strand, ref, alt)
## genome: A BSgenome object to extract the sequence information
MutScoreSS5 <- function(exs=Muts2Exons(mutations, exons(TxDb.Hsapiens.UCSC.hg19.knownGene)),
                        genome=BSgenome.Hsapiens.UCSC.hg19) {
  mutSeqs <- GetMutSeq(exs)
  seqLens <- nchar(mutSeqs)
  exs$seq <- substring(mutSeqs, seqLens - 8, seqLens)
  return(MaxentPerl(exs))
}

## Gets 3'SS Maxent scores for exons with associated mutations
## mutations: A dataframe containing mutation info (chr, pos, strand, ref, alt)
## genome: A BSgenome object to extract the sequence information
MutScoreSS3 <- function(exs=Muts2Exons(mutations, exons(TxDb.Hsapiens.UCSC.hg19.knownGene)),
                        genome=BSgenome.Hsapiens.UCSC.hg19) {
  mutSeqs <- GetMutSeq(exs)
  exs$seq <- substring(mutSeqs, 1, 23)
  return(MaxentPerl(exs, script="score3.pl"))
}

## Gets the absolute difference between mut/wt maxent 5'SS scores
MutWtDiffSS5 <- function(mutations, exs=exons(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
  goodExons <- Muts2Exons(mutations, exs)
  wtScores <- QuerySS5Scores(goodExons$exon_id)
  mutScores <- MutScoreSS5(goodExons)
  df <- merge(wtScores, mutScores, by="exon_id")
  return(data.frame(exon_id=df$exon_id, mwDiffScore=abs(df$score - df$ss5score)))
}

## Gets the absolute difference between mut/wt maxent 3'SS scores
MutWtDiffSS3 <- function(mutations, exs=exons(TxDb.Hsapiens.UCSC.hg19.knownGene)) {
  goodExons <- Muts2Exons(mutations, exs)
  wtScores <- QuerySS3Scores(goodExons$exon_id)
  mutScores <- MutScoreSS3(goodExons)
  df <- merge(wtScores, mutScores, by="exon_id")
  return(data.frame(exon_id=df$exon_id, mwDiffScore=abs(df$score - df$ss3score)))
}

## Calculate all mutation features for a list of mutations
## muts: A dataframe of mutation info in the same form as the output of ParseMutText
## Only returns mutations in exons that have splice site usage data
CalcMutFeatures <- function(exs, genome=BSgenome.Hsapiens.UCSC.hg19) {
#  exs <- Muts2Exons(muts)
#  exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

  ## Have to merge dataframes via sqldf because the regular merge function affects the ordering.
  ## This way, the exon IDs are in the same order as those in 'gr'
#  grMcols <- as.data.frame(mcols(exs))
#  usageMcols <- as.data.frame(mcols(exsWithSSUsage))
#  mcols(exs) <- sqldf(paste("SELECT grMcols.*, ss5usage, ss3usage, chasin_ese_density, chasin_ess_density",
#                            "FROM grMcols LEFT JOIN usageMcols",
#                            "ON grMcols.exon_id=usageMcols.exon_id;"))
#  exs <- exs[complete.cases(mcols(exs))]

  ## The actual mutation
  exs$mutation_base_change <- paste0(exs$ref, "->", exs$alt)

  ## Get wt seq and scores so we can generate the mutant ones
#  exs$wtSeq <- Rle(getSeq(genome, exs, as.character=T))
#  wtSS5Scores <- QuerySS5Scores(exs$exon_id)
#  wtSS3Scores <- QuerySS3Scores(exs$exon_id)

  ## Get mut sequence and scores
  mutSeqs <- GetMutSeq(exs, hasSpliceSites=F, relativePos=T)
  ## 5'SS
  seqLens <- nchar(mutSeqs)
  exs$seq <- substring(mutSeqs, seqLens - 8, seqLens)
  system.time(mutSS5Scores <- MaxentPerl(exs))
  ## 3'SS
  exs$seq <- substring(mutSeqs, 1, 23)
  system.time(mutSS3Scores <- MaxentPerl(exs, script="score3.pl"))

  exs$seq <- NULL

#  dfSS5 <- merge(wtSS5Scores, mutSS5Scores, by="exon_id")
#  dfSS3 <- merge(wtSS3Scores, mutSS3Scores, by="exon_id")
  exs$w5score <- exs$ss5score # for ML
  exs$m5score <- mutSS5Scores$score
  exs$mwdif_5score <- exs$w5score - exs$m5score
  exs$w3score <- exs$ss3score # for ML
  exs$m3score <- mutSS3Scores$score
  exs$mwdif_3score <- exs$w3score - exs$m3score

  return(exs)
}
