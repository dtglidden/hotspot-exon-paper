#!/usr/bin/env Rscript
## Generates a correlation heatmap for the following features:
## MaPSy WT Splicing Efficiency
## SS3Usage (psi3), SS5Usage (psi5) with nearest possible flanking exons
## SS5 and SS3 Maxent scores
## SS5score - flanking SS3score; SS3score - flanking SS5score
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_query.R"), chdir=T)
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
  library(corrplot)
})

mapsyExons <- readRDS(file.path("..", "data", "mapsy_features_gr.rds"))
mapsyExons <- mapsyExons[!duplicated(ranges(mapsyExons))]

mapsyExons$splEff <- MapsySplicingEfficiency(mapsyExons$hek_ws, mapsyExons$hek_wu)
mapsyExons$alSkew <- AllelicSkew(
  mapsyExons$hek_ws,
  mapsyExons$hek_wu,
  mapsyExons$hek_ms,
  mapsyExons$hek_mu)

exonData <- QueryExonsWithSSUsage()
mapsyExons <- AddMcols(mapsyExons, exonData)


## Plot the correlation matrix
featureNames <- c(
  "Allelic Imbalance"="alSkew",
  "Splicing Efficiency"="splEff",
  "WT Spliced Counts"="hek_ws",
  "WT Unspliced Counts"="hek_wu",
  "5'SS Maxent Score"="ss5score",
  "3'SS Maxent Score"="ss3score",
  "5'SS Usage"="ss5usage",
  "3'SS Usage"="ss3usage",
  "Avg EI"="ei",
  "Avg Delta EI"="avgDeltaEIScores",
  "Avg Delta A3SS"="avgDeltaA3SSScores",
  "Avg Delta A5SS"="avgDeltaA5SSScores",
  "# Exons"="nExons",
  "Gene Length"="geneLength"
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
