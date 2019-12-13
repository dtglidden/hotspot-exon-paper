#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
  source(file.path("..", "lib", "granges_util.R"), chdir=T)
  library(getopt)
  library(sqldf)
  library(ggplot2)
  library(ggsignif)
  library(dplyr)
})

opt <- getopt(matrix(c("test", "t", 0, "logical"),
                     byrow=T, ncol=4))
if (is.null(opt$test)) opt$test <- F

set.seed(1)

exsWithSSUsage <- if (opt$test) {
  readRDS(file.path("..", "example_data", "exons_with_ss_usage.rds"))
} else {
  QueryExonsWithSSUsage(as.GRanges=T)
}

## Add complete splice site regions to each exon
exsWithSSUsageAndSites <- AddSpliceSites(exsWithSSUsage)

## Get HGMD variants with associated SS usage values
hgmdGr <- if (opt$test) {
  Df2Gr(read.delim(file.path("..", "example_data", "splice_alleles_test.txt")))
} else {
  Df2Gr(read.delim(file.path("..", "data", "hgmd_splice.txt"), colClasses=c("acc_num"="character")))
}
hgmdMerge <- mergeByOverlaps(exsWithSSUsageAndSites, hgmdGr)
hgmdMerge <- hgmdMerge[!duplicated(hgmdMerge$hgmdGr), ]

## Get 'control' variants (HGMD variants within the critical GT/AG nucleotides)
hgmdCritGr <- if (opt$test) {
  Df2Gr(read.delim(file.path("..", "example_data", "splice_critical_alleles_test.txt")))
} else {
  Df2Gr(read.delim(file.path("..", "data", "hgmd_splice_critical.txt"), colClasses=c("acc_num"="character")))
}
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

## t-tests between HGMD variants and RBP knock-downs
## 'Func_num' > 0 => Sensitive; 'Func_num == 0 => Resistant
rbpUsage <- if (opt$test) {
  read.delim(file.path("..", "example_data", "HepG2_exons_by_func_rbp_num_usage.txt"))
} else {
  read.delim(file.path("..", "data", "HepG2_exons_by_func_rbp_num_usage.txt"))
}
tSens5SS <- t.test(hgmdMerge$ss5usage, rbpUsage$Fivep_usage[rbpUsage$Func_num > 0])
tRes5SS <- t.test(hgmdCritMerge$ss5usage, rbpUsage$Fivep_usage[rbpUsage$Func_num == 0])
tSens3SS <- t.test(hgmdMerge$ss3usage, rbpUsage$Threep_usage[rbpUsage$Func_num > 0])
tRes3SS <- t.test(hgmdCritMerge$ss3usage, rbpUsage$Threep_usage[rbpUsage$Func_num == 0])

## Incorporate RBP binding usage data
rbpUsage <- if (opt$test) {
  read.delim(file.path("..", "example_data", "ss_usage_by_binding.txt"), col.names=c("grp", "meanUsage", "seUsage"))
} else {
  read.delim(file.path("..", "data", "ss_usage_by_binding.txt"), col.names=c("grp", "meanUsage", "seUsage"))
}
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
yMin <- if (opt$test) 0 else 0.7
ggplot(dfSummary, aes(grp, meanUsage, fill=ss)) +
  geom_col(color="black") +
  geom_errorbar(aes(ymin=meanUsage-seUsage, ymax=meanUsage+seUsage), width=0.2) +
  coord_cartesian(ylim=c(yMin, 1)) +
  geom_signif(comparisons=list(c("Sensitive 5'SS", "Resistant 5'SS"),
                               c("Sensitive 3'SS", "Resistant 3'SS")),
              annotation="***",
              textsize=7,
              margin_top=0.2) +
  labs(x="", y="Splice Site Usage", fill="Splice Site") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text=element_text(size=24))
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
  labs(x="", y="Splice Site Usage", fill="Splice Site") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(color="white",
                                  margin=margin(t=0, r=20, b=0, l=0)),
        text=element_text(size=24))
ggsave(file.path("..", "plots", "usagePlotRBP.pdf"))

## Incorporate amiloride usage data
amiloUsageLst <- if (opt$test) {
  readRDS(file.path("..", "example_data", "amilo_usage.rds"))
} else {
  readRDS(file.path("..", "data", "amilo_usage.rds"))
}
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
  labs(x="", y="Splice Site Usage", fill="Splice Site") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(color="white",
                                  margin=margin(t=0, r=20, b=0, l=0)),
        text=element_text(size=24))
ggsave(file.path("..", "plots", "usagePlotAmilo.pdf"))
