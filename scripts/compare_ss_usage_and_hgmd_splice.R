#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  source(file.path("..", "lib", "granges_util.R"))
  library(sqldf)
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
})

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)
## Add complete splice site regions to each exon
exsWithSSUsageAndSites <- AddSpliceSites(exsWithSSUsage)

## Get HGMD variants with associated SS usage values
hgmdGr <- Df2Gr(read.delim(file.path("..", "data", "hgmd_splice.txt"), colClasses=c("acc_num"="character")))
hgmdMerge <- mergeByOverlaps(exsWithSSUsageAndSites, hgmdGr)
hgmdMerge <- hgmdMerge[!duplicated(hgmdMerge$hgmdGr), ]

## Get 'control' variants (HGMD variants within the critical GT/AG nucleotides)
hgmdCritGr <- Df2Gr(read.delim(file.path("..", "data", "hgmd_splice_critical.txt"), colClasses=c("acc_num"="character")))
hgmdCritMerge <- mergeByOverlaps(exsWithSSUsageAndSites, hgmdCritGr)
hgmdCritMerge <- hgmdCritMerge[!duplicated(hgmdCritMerge$hgmdCritGr), ]

## There are more control variants than experimental variants. Randomly select from control variants
nRow <- nrow(hgmdMerge)
randomIdxes <- sample(nrow(hgmdCritMerge), nRow)
hgmdCritMerge <- hgmdCritMerge[randomIdxes, ]

## Calculate statistics
avgSS5Usage <- mean(hgmdMerge$ss5usage)
avgSS3Usage <- mean(hgmdMerge$ss3usage)
avgCritSS5Usage <- mean(hgmdCritMerge$ss5usage)
avgCritSS3Usage <- mean(hgmdCritMerge$ss3usage)

## Plot
groupChars <- c(rep("Non-GT 5'SS SNVs", nRow),
                rep("GT 5'SS SNVs", nRow),
                rep("Non-AG 3'SS SNVs", nRow),
                rep("AG 3'SS SNVs", nRow))
boxDf <- data.frame(grp=factor(groupChars,
                               levels=unique(groupChars)),
                    usage=c(hgmdMerge$ss5usage,
                            hgmdCritMerge$ss5usage,
                            hgmdMerge$ss3usage,
                            hgmdCritMerge$ss3usage))
dfSummary <- boxDf %>% group_by(grp) %>%
  summarise(meanUsage=mean(usage),
            seUsage=sd(usage)/sqrt(n()))
ggplot(dfSummary, aes(grp, meanUsage)) +
  geom_col() +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.8, 1)) +
  labs(x="", y="Mean Usage (%)") +
  geom_signif(comparisons=list(c("Non-GT 5'SS SNVs", "GT 5'SS SNVs"),
                               c("Non-AG 3'SS SNVs", "AG 3'SS SNVs")),
              annotation="***") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        text=element_text(size=18))
ggsave(file.path("..", "plots", "usagePlot.pdf"))
