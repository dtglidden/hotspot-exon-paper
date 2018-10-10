#!/usr/bin/env Rscript
## Calculates the differences b/w wt/mut Maxent scores in HGMD alleles
## Only the SSMs outside of the canonical GT/AG dinucleotides
source(file.path("..", "lib", "feature_gen.R"))

hgmd <- read.csv(file.path("/", "home", "datasets", "HGMD_2016", "HGMD_Advanced_Substitutions.csv"), stringsAsFactors=F)
wtPat <- "([[:alpha:]]+)\\[([[:alpha:]])/[[:alpha:]]\\]([[:alpha:]]+)"
## The Sequence context contains the [wt/mut] and 30 nt of flanking sequence on each side
hgmd$wtSequenceContext <- sub(wtPat, "\\1\\2\\3", hgmd$sequence_context_hg19)

## Remove all entries that don't contain a splice site in their sequence context
hgmd <- hgmd[!grepl("^[[:upper:]]+$", hgmd$wtSequenceContext), ]

## Give the indices of the 5'SS in the sequence context (first intronic base), or -1 if there isn't one
hgmd$ss5Locs <- regexpr("(?<=[[:upper:]])[[:lower:]]", hgmd$wtSequenceContext, perl=T)
## Checks if the mutation is within the bounds of the 5'SS, but not the dinucleotide
hgmd$isSS5Mut <- (31 - hgmd$ss5Locs) %in% c(-3:-1, 2:5)

## Give the indices of the 3'SS in the sequence context (first intronic base), or -1 if there isn't one
hgmd$ss3Locs <- regexpr("[[:lower:]](?=[[:upper:]])", hgmd$wtSequenceContext, perl=T)
## Checks if the mutation is within the bounds of the 3'SS, but not the dinucleotide
hgmd$isSS3Mut <- (31 - hgmd$ss3Locs) %in% c(-19:-2, 1:3)

## Keep only the SSMs that aren't in the dinucleotides
hgmd <- hgmd[hgmd$isSS5Mut | hgmd$isSS3Mut, ]

## There are some weird entries in HGMD, so we should remove any intron motif that isn't the canonical GT/AG
intronPat <- "[[:upper:]]gt|ag[[:upper:]]"
hgmd <- hgmd[grepl(intronPat, hgmd$wtSequenceContext), ]
## TODO: Make sure the mutation isn't at the intron-exon boundary (intron-side), causing a lowercase letter to be uppercase mistakenly

## Extract the splice sites and get their Maxent scores
## First, split the dataframe by the type of splice site
hgmd.s <- split(hgmd, hgmd$isSS5Mut)
names(hgmd.s) <- vapply(names(hgmd.s), function(n) if (n == "TRUE") "ss5" else "ss3", "")

## Regex pattern to extract the mutant allele
mutPat <- "([[:alpha:]]+)\\[[[:alpha:]]/([[:alpha:]])\\]([[:alpha:]]+)"

## Calculate the wt/mut 5'SS scores
ss5wtSeq <- toupper(substring(hgmd.s$ss5$wtSequenceContext, hgmd.s$ss5$ss5Locs - 3, hgmd.s$ss5$ss5Locs + 5))
ss5wtScores <- MaxentPerl(ss5wtSeq)
hgmd.s$ss5$mutSequenceContext <- sub(mutPat, "\\1\\2\\3", hgmd.s$ss5$sequence_context_hg19)
ss5mutSeq <- toupper(substring(hgmd.s$ss5$mutSequenceContext, hgmd.s$ss5$ss5Locs - 3, hgmd.s$ss5$ss5Locs + 5))
ss5mutScores <- MaxentPerl(ss5mutSeq)

## Calculate the wt/mut 3'SS scores
ss3wtSeq <- toupper(substring(hgmd.s$ss3$wtSequenceContext, hgmd.s$ss3$ss3Locs - 19, hgmd.s$ss3$ss3Locs + 3))
ss3wtScores <- MaxentPerl(ss3wtSeq, script="score3.pl")
hgmd.s$ss3$mutSequenceContext <- sub(mutPat, "\\1\\2\\3", hgmd.s$ss3$sequence_context_hg19)
ss3mutSeq <- toupper(substring(hgmd.s$ss3$mutSequenceContext, hgmd.s$ss3$ss3Locs - 19, hgmd.s$ss3$ss3Locs + 3))
ss3mutScores <- MaxentPerl(ss3mutSeq, script="score3.pl")

## Get all exons, with associated splice sites scores and usages
#usage <- SSUsage("~/scratch/amiloride/star/ctrl/SJ.out.tab", exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
usage <- readRDS(file.path("..", "data", "ss_usage.rds"))

## Get the splice site sequences for all exons that have splice site usage data
allSS5Scores <- Score5SS(usage$ss5)
allSS3Scores <- Score3SS(usage$ss3)
usage$ss5$ssSeq <- allSS5Scores$seq
usage$ss5$ssScore <- allSS5Scores$score
usage$ss3$ssSeq <- allSS3Scores$seq
usage$ss3$ssScore <- allSS3Scores$score

## Get all exons that have normally have the splice site sequences of the HGMD variants, and those that don't
ss5ExonsWithSSV <- split(usage$ss5, usage$ss5$ssSeq %in% ss5mutSeq)
ss3ExonsWithSSV <- split(usage$ss3, usage$ss3$ssSeq %in% ss3mutSeq)

## Shows the ranges of maxent scores b/w the HGMD mut 5'SS sequences that are found in wt genomic exons
## and those that aren't found in any wt genomic exon
## NB: There are no mut 3'SS HGMD sequences found in wt genomic exons. This is probably because the 3'SS
## is 23 bases instead of 9, and therefore dramatically less likely to match another wt 3'SS
boxplot(list(SSV_in_genome=ss5ExonsWithSSV$"TRUE"$ssScore, SSV_not_in_genome=ss5ExonsWithSSV$"FALSE"$ssScore), ylab="5'SS MAXENT")
