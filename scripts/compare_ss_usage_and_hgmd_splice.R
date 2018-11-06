#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  library(sqldf)
  library(ggplot2)
  library(dplyr)
})

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

## Get HGMD variants with associated SS usage values
hgmdGr <- Df2Gr(read.delim(file.path("..", "data", "hgmd_splice.txt"), colClasses=c("acc_num"="character")))

hgmdMerge <- mergeByOverlaps(exsWithSSUsage, hgmdGr)
hgmdMerge <- hgmdMerge[!duplicated(hgmdMerge$hgmdGr), ]

## For each exon in hgmdMerge, get all of the exons in the transcripts associated them
## List index corresponds to exon_id
txByExon <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="exon")
txByExon <- txByExon[hgmdMerge$exon_id]
topTxs <- sapply(txByExon, function(x) x[1]$tx_name)
exonsByTx <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="tx", use.names=T)
## Filter for the transcripts that have an exon in txByExon. Have them ordered by the exons in txByExon
exonsByTxFilt <- exonsByTx[names(exonsByTx) %in% topTxs]
grlLengths <- sapply(exonsByTxFilt, length)
## Get list of exon_ids for each transcript (faster to iterate through than GRangesList
exonsByTxList <- lapply(exonsByTxFilt, function(x) x$exon_id)
## Get indices in exonsByTxList that correspond to the exons in the order of hgmdMerge$exon_id
indices <- sapply(hgmdMerge$exon_id, function(x) which(sapply(exonsByTxList, function(y) x %in% y))[[1]])
organizedExonsBySSVExon <- exonsByTxFilt[indices]
organizedExonsWithSSUsage <- lapply(organizedExonsBySSVExon, function(x) mergeByOverlaps(x, exsWithSSUsage))
avgSS5Usage <- sapply(organizedExonsWithSSUsage, function(x) mean(x$ss5usage))
avgSS3Usage <- sapply(organizedExonsWithSSUsage, function(x) mean(x$ss3usage))
medSS5Usage <- sapply(organizedExonsWithSSUsage, function(x) median(x$ss5usage))
medSS3Usage <- sapply(organizedExonsWithSSUsage, function(x) median(x$ss3usage))
minSS5Usage <- sapply(organizedExonsWithSSUsage, function(x) min(x$ss5usage))
minSS3Usage <- sapply(organizedExonsWithSSUsage, function(x) min(x$ss3usage))
usageDf <- data.frame(ss5usage=hgmdMerge$ss5usage, ss3usage=hgmdMerge$ss3usage, meanSS5usage=avgSS5Usage, meanSS3usage=avgSS3Usage, medSS5usage=medSS5Usage, medSS3usage=medSS3Usage)

nRow <- length(indices)
groupChars <- c(rep("HGMD SS5 Usage", nRow),
                rep("Mean Ctrl SS5 Usage", nRow),
                rep("Median Ctrl SS5 Usage", nRow),
                rep("HGMD SS3 Usage", nRow),
                rep("Mean Ctrl SS3 Usage", nRow),
                rep("Median Ctrl SS3 Usage", nRow))
boxDf <- data.frame(grp=factor(groupChars,
                               levels=unique(groupChars)),
                    usage=c(hgmdMerge$ss5usage,
                            avgSS5Usage,
                            medSS5Usage,
                            hgmdMerge$ss3usage,
                            avgSS3Usage,
                            medSS3Usage))
dfSummary <- boxDf %>% group_by(grp) %>%
  summarise(meanUsage=mean(usage),
            seUsage=sd(usage)/sqrt(n()))
ggplot(dfSummary, aes(grp, meanUsage)) +
  geom_col() +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.8, 1)) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        text=element_text(size=18))
ggsave(file.path("..", "plots", "usagePlot.pdf"))
