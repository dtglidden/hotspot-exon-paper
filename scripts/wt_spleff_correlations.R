#!/usr/bin/env Rscript
## Generates a correlation heatmap for the following features:
## MaPSy WT Splicing Efficiency
## SS3Usage (psi3), SS5Usage (psi5) with nearest possible flanking exons
## SS5 and SS3 Maxent scores
## SS5score - flanking SS3score; SS3score - flanking SS5score
source(file.path("..", "lib", "feature_query.R"))
source(file.path("..", "lib", "feature_gen.R"))
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
                      vivo_ref_unspliced=mutData$vivo_ref_unspliced)
## Remove duplicates because this dataset has entries for each mutation, which may be in the same exon
## We only care about the WT exons
mapsyExons <- mapsyExons[!duplicated(ranges(mapsyExons))]

ov <- findOverlaps(mapsyExons, exs, type="equal")
ids <- mcols(exs)[subjectHits(ov), "exon_id"]
mapsyExons <- mapsyExons[queryHits(ov)]
mapsyExons$exon_id <- ids

## Get WT Splicing Efficiency
wsCounts <- mcols(mapsyExons)$vivo_ref_spliced + 1 #pseudocounts
wsSum <- sum(wsCounts)
wsEff <- wsCounts / wsSum

wuCounts <- mcols(mapsyExons)$vivo_ref_unspliced + 1 #pseudocounts
wuSum <- sum(wuCounts)
wuEff <- wuCounts / wuSum

mcols(mapsyExons)$splEff <- log2(wsEff/wuEff)


exonData <- QueryExonsWithSSUsage()
mapsyExons <- AddMcols(mapsyExons, exonData)
mapsyExons$ssScoreSum <- mapsyExons$ss5score + mapsyExons$ss3score
mapsyExons$ssScoreDiff <- mapsyExons$ss5score - mapsyExons$ss3score
mapsyExons$ssScoreAbsDiff <- abs(mapsyExons$ss5score - mapsyExons$ss3score)
mapsyExons$ssUsageSum <- mapsyExons$ss5usage + mapsyExons$ss3usage
mapsyExons$ssUsageDiff <- mapsyExons$ss5usage - mapsyExons$ss3usage
mapsyExons$ssUsageAbsDiff <- abs(mapsyExons$ss5usage - mapsyExons$ss3usage)


## Plot the correlation matrix
corMat <- cor(as.data.frame(mcols(mapsyExons)[, c(
  "splEff",
  "ss5score",
  "ss3score",
  "ssScoreSum",
  "ssScoreDiff",
  "ssScoreAbsDiff",
  "ss5usage",
  "ss3usage",
  "ssUsageSum",
  "ssUsageDiff",
  "ssUsageAbsDiff",
  "chasin_ese_density",
  "chasin_ess_density"
)]))

colPal <- colorRampPalette(rev(c("#7F0000", "red", "#FF7F00", "yellow", "white",
        "cyan", "#007FFF", "blue", "#00007F")))
pdf(file.path("..", "lib", "feature_cor.pdf"))
corrplot(corMat, method="ellipse", col=colPal(200))
dev.off()
