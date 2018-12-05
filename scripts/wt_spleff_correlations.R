#!/usr/bin/env Rscript
## Generates a correlation heatmap for the following features:
## MaPSy WT Splicing Efficiency
## SS3Usage (psi3), SS5Usage (psi5) with nearest possible flanking exons
## SS5 and SS3 Maxent scores
## SS5score - flanking SS3score; SS3score - flanking SS5score
suppressPackageStartupMessages({
source(file.path("..", "lib", "feature_query.R"))
source(file.path("..", "lib", "feature_gen.R"))
library(corrplot)
})

exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
mutData <- read.delim(file.path("..", "data", "input_ML_complete_table.txt"))
sjFile <- file.path("~", "scratch", "amiloride", "star", "ctrl", "SJ.out.tab")
mapsyExons <- GRanges(seqnames=mutData$chrom_exon,
                      ranges=IRanges(start=mutData$start_exon, end=mutData$end_exon),
                      strand=mutData$strand,
                      seqinfo=Seqinfo(genome="hg19"),
                      vivo_ref_spliced=mutData$vivo_ref_spliced,
                      vivo_ref_unspliced=mutData$vivo_ref_unspliced,
                      vivo_alt_spliced=mutData$vivo_alt_spliced,
                      vivo_alt_unspliced=mutData$vivo_alt_unspliced)
## Remove duplicates because this dataset has entries for each mutation, which may be in the same exon
## We only care about the WT exons
mapsyExons <- mapsyExons[!duplicated(ranges(mapsyExons))]

ov <- findOverlaps(mapsyExons, exs, type="equal")
ids <- mcols(exs)[subjectHits(ov), "exon_id"]
mapsyExons <- mapsyExons[queryHits(ov)]
mapsyExons$exon_id <- ids

mapsyExons$splEff <- MapsySplicingEfficiency(mapsyExons$vivo_ref_spliced, mapsyExons$vivo_ref_unspliced)
mapsyExons$alSkew <- AllelicSkew(
  mapsyExons$vivo_ref_spliced,
  mapsyExons$vivo_ref_unspliced,
  mapsyExons$vivo_alt_spliced,
  mapsyExons$vivo_alt_unspliced)


exonData <- QueryExonsWithSSUsage()
mapsyExons <- AddMcols(mapsyExons, exonData)
mapsyExons$ssScoreSum <- mapsyExons$ss5score + mapsyExons$ss3score
mapsyExons$ssScoreDiff <- mapsyExons$ss5score - mapsyExons$ss3score
mapsyExons$ssScoreAbsDiff <- abs(mapsyExons$ss5score - mapsyExons$ss3score)
mapsyExons$ssUsageSum <- mapsyExons$ss5usage + mapsyExons$ss3usage
mapsyExons$ssUsageDiff <- mapsyExons$ss5usage - mapsyExons$ss3usage
mapsyExons$ssUsageAbsDiff <- abs(mapsyExons$ss5usage - mapsyExons$ss3usage)


## Plot the correlation matrix
featureNames <- c(
  "Allelic Skew"="alSkew",
  "Splicing Efficiency"="splEff",
  "WT Spliced Counts"="vivo_ref_spliced",
  "WT Unspliced Counts"="vivo_ref_unspliced",
  "5'SS Maxent Score"="ss5score",
  "3'SS Maxent Score"="ss3score",
  "Maxent Score Sum"="ssScoreSum",
  "Maxent Score Diff"="ssScoreDiff",
  "Maxent Score Abs Diff"="ssScoreAbsDiff",
  "5'SS Usage"="ss5usage",
  "3'SS Usage"="ss3usage",
  "SS Usage Sum"="ssUsageSum",
  "SS Usage Diff"="ssUsageDiff",
  "SS Usage Abs Diff"="ssUsageAbsDiff",
  "ESE Density"="chasin_ese_density",
  "ESS Density"="chasin_ess_density"
)

corMat <- cor(as.data.frame(mcols(mapsyExons)[, featureNames]))
rownames(corMat) <- names(featureNames)
colnames(corMat) <- names(featureNames)

colPal <- colorRampPalette(c(
  "#00007F",
  "blue",
  "#007FFF",
  "cyan",
  "white",
  "yellow",
  "#FF7F00",
  "red",
  "#7F0000"
))

pdf(file.path("..", "plots", "feature_cor.pdf"))
corrplot(corMat, method="color", col=colPal(200))
dev.off()
