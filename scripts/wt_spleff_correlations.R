#!/usr/bin/env Rscript
## Generates a correlation heatmap for the following features:
## MaPSy WT Splicing Efficiency
## SS3Usage (psi3), SS5Usage (psi5) with nearest possible flanking exons
## SS5 and SS3 Maxent scores
## SS5score - flanking SS3score; SS3score - flanking SS5score
source("feature_gen.R")
library(corrplot)

GetPsis <- function(fIntrons, fExons) {
  exonsWithPsi5DF <- mergeByOverlaps(fExons, flank(fIntrons, 1))
  exonsWithPsi5 <- exonsWithPsi5DF$fExons
  mcols(exonsWithPsi5)$psi5 <- exonsWithPsi5DF$psi5
  mcols(exonsWithPsi5)$psi3flank <- exonsWithPsi5DF$psi3
  exonsWithPsi3DF <- mergeByOverlaps(fExons, flank(fIntrons, 1, start=F))
  exonsWithPsi3 <- exonsWithPsi3DF$fExons
  mcols(exonsWithPsi3)$psi3 <- exonsWithPsi3DF$psi3
  mcols(exonsWithPsi3)$psi5flank <- exonsWithPsi3DF$psi5
  return(merge(exonsWithPsi5, exonsWithPsi3))
}

exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
mutData <- read.delim(file.path("..", "data", "input_ML_complete_table.txt"))
sjFile <- file.path("~", "scratch", "amiloride", "star", "ctrl", "SJ.out.tab")
mapsyExons <- GRanges(seqnames=mutData$chrom_exon,
                      ranges=IRanges(start=mutData$start_exon, end=mutData$end_exon),
                      strand=mutData$strand,
                      seqinfo=Seqinfo(genome="hg19"),
                      vivo_ref_spliced=mutData$vivo_ref_spliced,
                      vivo_ref_unspliced=mutData$vivo_ref_unspliced,
                      ss5seq=mutData$splice_site_5_ref,
                      ss5score=mutData$splice_site_5_ref_score,
                      ss3seq=mutData$splice_site_3_ref,
                      ss3score=mutData$splice_site_3_ref_score)
## Remove duplicates because this dataset has entries for each mutation, which may be in the same exon
## We only care about the WT exons
mapsyExons <- mapsyExons[!duplicated(ranges(mapsyExons))]

## Get WT Splicing Efficiency
wsCounts <- mcols(mapsyExons)$vivo_ref_spliced + 1 #pseudocounts
wsSum <- sum(wsCounts)
wsEff <- wsCounts / wsSum

wuCounts <- mcols(mapsyExons)$vivo_ref_unspliced + 1 #pseudocounts
wuSum <- sum(wuCounts)
wuEff <- wuCounts / wuSum

mcols(mapsyExons)$splEff <- log2(wsEff/wuEff)

introns <- SSUsage(sjFile, exs)

goodExons <- GetPsis(introns, mapsyExons)
goodExons$psi3flank_minus_psi5 <- goodExons$psi3flank - goodExons$psi5
goodExons$psi3flank_minus_psi5_abs <- abs(goodExons$psi3flank - goodExons$psi5)
goodExons$psi5flank_minus_psi3 <- goodExons$psi5flank - goodExons$psi3
goodExons$psi5flank_minus_psi3_abs <- abs(goodExons$psi5flank - goodExons$psi3)
goodExons$ssScoreSum <- goodExons$ss5score + goodExons$ss3score
goodExons$ssScoreDiff <- goodExons$ss5score - goodExons$ss3score
goodExons$ssScoreAbsDiff <- abs(goodExons$ss5score - goodExons$ss3score)

allExons <- GetPsis(introns, exs)
## Which exons in the total list of exons are preceded by our good exons (MaPSy exons)
## I.e. get the downstream exons and their associated splice site usage info
dwnExons <- allExons[precede(goodExons, allExons)]
## Get the Maxent scores of the 5'SSs in these exons (competing splice sites)
goodExons$ss5score_competing <- Score5SS(dwnExons)$score
goodExons$ss5score_competing_diff <- goodExons$ss5score - goodExons$ss5score_competing

## Now get the upstream exons
upExons <- allExons[follow(goodExons, allExons)]
## Get the Maxent scores of the 3'SSs in these exons (competing splice sites)
goodExons$ss3score_competing <- Score3SS(upExons)$score
goodExons$ss3score_competing_diff <- goodExons$ss3score - goodExons$ss3score_competing

## Sum of absolute differences (SAD) b/w competing splice site scores
goodExons$maxent_competing_sad <- abs(goodExons$ss5score_competing_diff) + abs(goodExons$ss3score_competing_diff)

## Plot the correlation matrix
corMat <- cor(as.data.frame(mcols(goodExons)[, c("splEff", "ss5score", "ss3score", "ssScoreSum", "ssScoreDiff", "ssScoreAbsDiff", "ss5score_competing", "ss5score_competing_diff", "ss3score_competing", "ss3score_competing_diff", "maxent_competing_sad", "psi5", "psi3", "psi5flank", "psi3flank", "psi5flank_minus_psi3", "psi5flank_minus_psi3_abs", "psi3flank_minus_psi5", "psi3flank_minus_psi5_abs")]))

colPal <- colorRampPalette(rev(c("#7F0000", "red", "#FF7F00", "yellow", "white",
        "cyan", "#007FFF", "blue", "#00007F")))
#pdf("feature_cor.pdf")
corrplot(corMat, method="ellipse", col=colPal(200))
#dev.off()
