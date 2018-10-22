#!/usr/bin/env Rscript
## Calculate mutation features for MaPSy alleles

suppressPackageStartupMessages({
source(file.path("..", "lib", "feature_gen.R"))
library(randomForest)
library(ROCR)
library(ggplot2)
})

df <- read.csv(file.path("..", "data", "pMaster_020615_1.csv"))
df <- df[grepl("^HGMD", df$classifier), ]
nRows <- nrow(df)
splitIdx <- floor(0.7 * nRows)
trainIndices <- 1:splitIdx
testIndices <- (splitIdx+1):nRows
gr <- GRanges(paste(sub("-", ":", df$ex_hg19ann), df$strand, sep=":"))
allelic_skew <- AllelicSkew(df$hek_ws, df$hek_wu, df$hek_ms, df$hek_mu)
affects_splicing <- as.factor(allelic_skew >= log2(1.5)) # p-value?
response <- affects_splicing[trainIndices]
gr$mutation_base_change <- as.factor(paste0(df$allele_w, df$allele_m))
gr$w5score <- df$w5score
gr$w3score <- df$w3score
gr$m5score <- df$m5score
gr$m3score <- df$m3score
gr$mwdif_5score <- df$mwdif_5score
gr$mwdif_3score <- df$mwdif_3score

train <- as.data.frame(mcols(gr))[trainIndices, ]
test <- as.data.frame(mcols(gr))[testIndices, ]

model <- randomForest(train, response, importance=T)
predNums <- predict(model, test, type="prob")

realityLabels <- affects_splicing[testIndices]
rocrPred <- prediction(predNums[, 2], realityLabels)
auc <- performance(rocrPred, "auc")
myRoc <- performance(rocrPred, "tpr", "fpr")
ggDf <- data.frame(tpr=myRoc@x.values[[1]], fpr=myRoc@y.values[[1]], auc=auc@y.values[[1]])
pdf("mapsy_roc.pdf")
ggplot(ggDf, aes(tpr, fpr, color=auc)) +
  geom_line(size=2, alpha=0.7) +
  labs(x="False Positive Rate (1-Specificity)",
    y="True Positive Rate (Sensitivity)") +
  theme(text=element_text(size=20))
dev.off()
