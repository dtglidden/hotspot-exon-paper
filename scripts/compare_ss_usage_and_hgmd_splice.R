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
rbpUsage$grp <- c("Resistant 3'SS", #3'ss
                  "Resistant 5'SS", #5'ss
                  "Sensitive 3'SS", #3'ss
                  "Sensitive 5'SS") #5'ss

SE <- function(x) sqrt(var(x)/length(x))

## Plot
groupChars <- c(rep("Sensitive 5'SS", nRow),
                rep("Resistant 5'SS", nRow),
                rep("Sensitive 3'SS", nRow),
                rep("Resistant 3'SS", nRow))
boxDf <- data.frame(grp=factor(groupChars,
                               levels=unique(groupChars)),
                    usage=c(hgmdMerge$ss5usage,
                            hgmdCritMerge$ss5usage,
                            hgmdMerge$ss3usage,
                            hgmdCritMerge$ss3usage))
dfSummary <- boxDf %>% group_by(grp) %>%
  summarise(meanUsage=mean(usage),
            seUsage=sd(usage)/sqrt(n()))
dfSummary$grp <- factor(dfSummary$grp, levels=as.character(dfSummary$grp))
halfNRow <- nrow(dfSummary) / 2
dfSummary$ss <- factor(c(rep("5'SS", halfNRow), rep("3'SS", halfNRow)), levels=c("5'SS", "3'SS"))
ggplot(dfSummary, aes(grp, meanUsage, fill=ss)) +
  geom_col(color="black") +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.7, 1)) +
  geom_signif(comparisons=list(c("Sensitive 5'SS", "Resistant 5'SS"),
                               c("Sensitive 3'SS", "Resistant 3'SS")),
              annotation="***",
              textsize=7,
              margin_top=0.2) +
  labs(x="", y="Splice Site Usage (%)", fill="Splice Site") +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text=element_text(size=24),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
ggsave(file.path("..", "plots", "usagePlotHGMD.pdf"))

## RBP plot
## Reorder the dataframe so all the 5'SS elements come before the 3'SS ones
dfSummary <- rbind(rbpUsage[c(4, 2), ],
                   rbpUsage[c(3, 1), ])
dfSummary$grp <- factor(dfSummary$grp, levels=as.character(unique(dfSummary$grp)))
halfNRow <- nrow(dfSummary) / 2
dfSummary$ss <- factor(c(rep("5'SS", halfNRow), rep("3'SS", halfNRow)), levels=c("5'SS", "3'SS"))
ggplot(dfSummary, aes(grp, meanUsage, fill=ss)) +
  geom_col(color="black") +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.7, 1)) +
  geom_signif(comparisons=list(c("Resistant 5'SS", "Sensitive 5'SS"),
                               c("Resistant 3'SS", "Sensitive 3'SS")),
              annotation="***",
              textsize=7,
              margin_top=0.05) +
  labs(x="", y="Splice Site Usage (%)", fill="Splice Site") +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text=element_text(size=24),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
ggsave(file.path("..", "plots", "usagePlotRBP.pdf"))

## Incorporate amiloride usage data
amiloUsageLst <- readRDS(file.path("..", "data", "amilo_usage.rds"))
amiloUsage <- data.frame(
  grp=c("Sensitive 5'SS",
        "Resistant 5'SS",
        "Sensitive 3'SS",
        "Resistant 3'SS"),
  meanUsage=c(mean(amiloUsageLst$amilo_ss5_usage$usage),
              mean(amiloUsageLst$amilo_ss5_non_usage$usage),
              mean(amiloUsageLst$amilo_ss3_usage$usage),
              mean(amiloUsageLst$amilo_ss3_non_usage$usage)),
  seUsage=c(SE(amiloUsageLst$amilo_ss5_usage$usage),
              SE(amiloUsageLst$amilo_ss5_non_usage$usage),
              SE(amiloUsageLst$amilo_ss3_usage$usage),
              SE(amiloUsageLst$amilo_ss3_non_usage$usage)))
dfSummary <- amiloUsage
dfSummary$grp <- factor(dfSummary$grp, levels=as.character(unique(dfSummary$grp)))
halfNRow <- nrow(dfSummary) / 2
dfSummary$ss <- factor(c(rep("5'SS", halfNRow), rep("3'SS", halfNRow)), levels=c("5'SS", "3'SS"))
ggplot(dfSummary, aes(grp, meanUsage, fill=ss)) +
  geom_col(color="black") +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(0.7, 1)) +
  geom_signif(comparisons=list(c("Resistant 5'SS", "Sensitive 5'SS"),
                               c("Resistant 3'SS", "Sensitive 3'SS")),
              annotation="***",
              textsize=7,
              margin_top=0.2) +
  labs(x="", y="Splice Site Usage (%)", fill="Splice Site") +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        text=element_text(size=24),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
ggsave(file.path("..", "plots", "usagePlotAmilo.pdf"))
