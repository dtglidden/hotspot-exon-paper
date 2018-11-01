#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  library(sqldf)
})

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

## Get HGMD variants with associated SS usage values
hgmdGr <- Df2Gr(read.delim(file.path("..", "data", "hgmd_splice.txt"), colClasses=c("acc_num"="character")))

hgmdMerge <- mergeByOverlaps(exsWithSSUsage, hgmdGr)
hgmdMerge <- hgmdMerge[!duplicated(hgmdMerge$hgmdGr), ]

## Get control variants (from ExAC) with associated SS usage values
## 5'SS variants
ctrl5Gr <- Df2Gr(read.delim(file.path("..", "data", "ExAC_ss5vars"), header=F,
                            col.names=c("chrom",
                                        "start",
                                        "end",
                                        "snp",
                                        "ref",
                                        "alt",
                                        "chrom2",
                                        "ssStart",
                                        "ssEnd",
                                        "refGene",
                                        "score",
                                        "strand"),
                            colClasses=c("snp"="character")))
ctrl5Gr <- ctrl5Gr[nchar(as.character(ctrl5Gr$ref)) == 1 &
                     nchar(as.character(ctrl5Gr$alt)) == 1]
start(ctrl5Gr) <- start(ctrl5Gr) + 1

ctrl5Merge <- mergeByOverlaps(exsWithSSUsage, ctrl5Gr)
ctrl5Merge <- ctrl5Merge[!duplicated(ctrl5Merge$ctrl5Gr), ]
dsIdx <- sample(1:nrow(ctrl5Merge), nrow(hgmdMerge))
ctrl5MergeDS <- ctrl5Merge[dsIdx, ]

## 3'SS variants
ctrl3Gr <- Df2Gr(read.delim(file.path("..", "data", "ExAC_ss3vars"), header=F,
                            col.names=c("chrom",
                                        "start",
                                        "end",
                                        "snp",
                                        "ref",
                                        "alt",
                                        "chrom2",
                                        "ssStart",
                                        "ssEnd",
                                        "refGene",
                                        "score",
                                        "strand"),
                            colClasses=c("snp"="character")))
ctrl3Gr <- ctrl3Gr[nchar(as.character(ctrl3Gr$ref)) == 1 &
                     nchar(as.character(ctrl3Gr$alt)) == 1]
start(ctrl3Gr) <- start(ctrl3Gr) + 1

ctrl3Merge <- mergeByOverlaps(exsWithSSUsage, ctrl3Gr)
ctrl3Merge <- ctrl3Merge[!duplicated(ctrl3Merge$ctrl3Gr), ]
dsIdx <- sample(1:nrow(ctrl3Merge), nrow(hgmdMerge))
ctrl3MergeDS <- ctrl3Merge[dsIdx, ]

tTest5 <- t.test(hgmdMerge$ss5usage, ctrl5MergeDS$ss5usage)
tTest5Maxent <- t.test(hgmdMerge$ss5score, ctrl5MergeDS$ss5score)
tTest3 <- t.test(hgmdMerge$ss3usage, ctrl3MergeDS$ss3usage)
tTest3Maxent <- t.test(hgmdMerge$ss3score, ctrl3MergeDS$ss3score)
