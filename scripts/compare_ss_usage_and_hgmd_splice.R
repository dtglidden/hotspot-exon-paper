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

## Incorporate RBP binding usage data
rbpUsage <- read.delim(file.path("..", "data", "ss_usage_by_binding.txt"), col.names=c("grp", "meanUsage", "seUsage"))
rbpUsage$meanUsage <- rbpUsage$meanUsage / 100
rbpUsage$seUsage <- rbpUsage$seUsage / 100
rbpUsage$grp <- c("3'SS Binding",
                  "5'SS Binding",
                  "3'SS No Binding",
                  "5'SS No Binding")

SE <- function(x) sqrt(var(x)/length(x))

## Incorporate amiloride usage data
amiloUsageLst <- readRDS(file.path("..", "data", "amilo_usage.rds"))
#amiloUsage <- data.frame(meanSS5Usage=,
#                         meanSS3Usage=,
#                         meanUsSS5Usage=,
#                         meanDsSS3Usage=mean(amiloUsageLst$amilo_ds_ss3_usage$usage),
#                         meanSS5NonUsage=,
#                         meanSS3NonUsage=,
#                         meanUsSS5NonUsage=,
#                         meanDsSS3NonUsage=mean(amiloUsageLst$amilo_ds_ss3_non_usage$usage),
#                         seSS5Usage=SE(amiloUsageLst$amilo_ss5_usage$usage),
#                         seSS3Usage=SE(amiloUsageLst$amilo_ss3_usage$usage),
#                         seUsSS5Usage=SE(amiloUsageLst$amilo_us_ss5_usage$usage),
#                         seDsSS3Usage=SE(amiloUsageLst$amilo_ds_ss3_usage$usage),
#                         seSS5NonUsage=SE(amiloUsageLst$amilo_ss5_non_usage$usage),
#                         seSS3NonUsage=SE(amiloUsageLst$amilo_ss3_non_usage$usage),
#                         seUsSS5NonUsage=SE(amiloUsageLst$amilo_us_ss5_non_usage$usage),
#                         seDsSS3NonUsage=SE(amiloUsageLst$amilo_ds_ss3_non_usage$usage))
amiloUsage <- data.frame(
  grp=c("ss5Amilo",
        "ss5NonAmilo",
        "ss3Amilo",
        "ss3NonAmilo"),
  meanUsage=c(mean(amiloUsageLst$amilo_ss5_usage$usage),
              mean(amiloUsageLst$amilo_ss5_non_usage$usage),
              mean(amiloUsageLst$amilo_ss3_usage$usage),
              mean(amiloUsageLst$amilo_ss3_non_usage$usage)),
  seUsage=c(SE(amiloUsageLst$amilo_ss5_usage$usage),
              SE(amiloUsageLst$amilo_ss5_non_usage$usage),
              SE(amiloUsageLst$amilo_ss3_usage$usage),
              SE(amiloUsageLst$amilo_ss3_non_usage$usage)))

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
## Reorder the dataframe so all the 5'SS elements come before the 3'SS ones
dfSummary <- rbind(dfSummary[1:2, ],
                   rbpUsage[c(2, 4), ],
                   amiloUsage[1:2, ],
                   dfSummary[3:4, ],
                   rbpUsage[c(1, 3), ],
                   amiloUsage[3:4, ])
dfSummary$grp <- factor(dfSummary$grp, levels=as.character(dfSummary$grp))
halfNRow <- nrow(dfSummary) / 2
dfSummary$ss <- factor(c(rep("5'SS", halfNRow), rep("3'SS", halfNRow)), levels=c("5'SS", "3'SS"))
ggplot(dfSummary, aes(grp, meanUsage, fill=ss)) +
  geom_col() +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.8, 1)) +
  geom_signif(comparisons=list(c("Non-GT 5'SS SNVs", "GT 5'SS SNVs"),
                               c("Non-AG 3'SS SNVs", "AG 3'SS SNVs")),
              annotation="***") +
  labs(x="", y="Splice Site Usage (%)", fill="Splice Site") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=18))
ggsave(file.path("..", "plots", "usagePlot.pdf"))
