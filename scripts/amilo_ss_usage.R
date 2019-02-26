#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
})

usage <- readRDS(file.path("..", "data", "ss_usage.rds"))

df <- read.delim(file.path("..", "data", "amilo_se_rmats.txt"))
df <- filter(df,
             PValue < 0.05,
             abs(IncLevelDifference) > 0.15)

sq <- Seqinfo(genome="hg19")

gr <- GRanges(seqnames=df$chr,
              ranges=IRanges(start=df$exonStart_0base + 1,
                             end=df$exonEnd),
              strand=df$strand,
              seqinfo=sq,
              gene=df$GeneID,
              deltaPsi=-df$IncLevelDifference)

usGr <- GRanges(seqnames=df$chr,
                ranges=IRanges(start=df$upstreamES,
                               end=df$upstreamEE),
                strand=df$strand,
                seqinfo=sq,
                gene=df$GeneID)

dsGr <- GRanges(seqnames=df$chr,
                ranges=IRanges(start=df$downstreamES + 1,
                               end=df$downstreamEE),
                strand=df$strand,
                seqinfo=sq,
                gene=df$GeneID)

## Get the usage of exons targeted by amiloride
amilo_ss5_usage <- mergeByOverlaps(gr, usage$ss5, type="end")
amilo_ss3_usage <- mergeByOverlaps(gr, usage$ss3, type="start")
amilo_us_ss5_usage <- mergeByOverlaps(usGr, usage$ss5, type="end")
amilo_ds_ss3_usage <- mergeByOverlaps(dsGr, usage$ss3, type="start")

## Get the usage of exons not targeted by amiloride
amilo_ss5_non_usage <- usage$ss5[!usage$ss5 %in% amilo_ss5_usage$usage.ss5]
amilo_ss3_non_usage <- usage$ss3[!usage$ss3 %in% amilo_ss3_usage$usage.ss3]
amilo_us_ss5_non_usage <- usage$ss5[!usage$ss5 %in% amilo_us_ss5_usage$usage.ss5]
amilo_ds_ss3_non_usage <- usage$ss3[!usage$ss3 %in% amilo_ds_ss3_usage$usage.ss3]

saveRDS(list(amilo_ss5_usage=amilo_ss5_usage,
             amilo_ss3_usage=amilo_ss3_usage,
             amilo_us_ss5_usage=amilo_us_ss5_usage,
             amilo_ds_ss3_usage=amilo_ds_ss3_usage,
             amilo_ss5_non_usage=amilo_ss5_non_usage,
             amilo_ss3_non_usage=amilo_ss3_non_usage,
             amilo_us_ss5_non_usage=amilo_us_ss5_non_usage,
             amilo_ds_ss3_non_usage=amilo_ds_ss3_non_usage),
        file.path("..", "data", "amilo_usage.rds"))
