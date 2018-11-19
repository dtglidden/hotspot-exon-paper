#!/usr/bin/env Rscript
## Calculate mutation features for MaPSy alleles

suppressPackageStartupMessages({
source(file.path("..", "lib", "feature_gen.R"))
library(sqldf)
library(randomForest)
library(ROCR)
library(ggplot2)
library(RColorBrewer)
})

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

df <- read.csv(file.path("..", "data", "pMaster_020615_1.csv"))
df <- df[grepl("^HGMD", df$classifier), ]

gr <- GRanges(paste(sub("-", ":", df$ex_hg19ann), df$strand, sep=":"))
start(gr) <- start(gr) + 1
gr$hek_ws <- df$hek_ws
gr$hek_wu <- df$hek_wu
gr$hek_ms <- df$hek_ms
gr$hek_mu <- df$hek_mu
gr$mutation_base_change <- as.factor(paste0(df$allele_w, df$allele_m))
gr$w5score <- df$w5score
gr$w3score <- df$w3score
gr$m5score <- df$m5score
gr$m3score <- df$m3score
gr$mwdif_5score <- df$mwdif_5score
gr$mwdif_3score <- df$mwdif_3score
gr <- AddExonId(gr)

## Have to merge dataframes via sqldf because the regular merge function affects the ordering.
## This way, the exon IDs are in the same order as those in 'gr'
grMcols <- as.data.frame(mcols(gr))
usageMcols <- as.data.frame(mcols(exsWithSSUsage))
mcols(gr) <- sqldf(paste("SELECT grMcols.*, ss5usage, ss3usage, chasin_ese_density, chasin_ess_density",
                         "FROM grMcols LEFT JOIN usageMcols",
                         "ON grMcols.exon_id=usageMcols.exon_id;"))
gr <- gr[complete.cases(mcols(gr))]

nRows <- length(gr)
splitIdx <- floor(0.7 * nRows)
trainIndices <- 1:splitIdx
testIndices <- (splitIdx+1):nRows

allelic_skew <- AllelicSkew(gr$hek_ws, gr$hek_wu, gr$hek_ms, gr$hek_mu)
affects_splicing <- as.factor(allelic_skew >= log2(1.5)) # p-value?
gr$affects_splicing <- affects_splicing
response <- affects_splicing[trainIndices]

## Variations of the data, so we can compare how some features affect AUC
mlDf1 <- as.data.frame(mcols(gr))[, c(
  "hek_ws",
  "hek_wu",
  "mutation_base_change",
  "w5score",
  "w3score",
  "m5score",
  "m3score",
  "mwdif_5score",
  "mwdif_3score",
  "ss5usage",
  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

mlDfNoMuts <- as.data.frame(mcols(gr))[, c(
  "hek_ws",
  "hek_wu",
#  "mutation_base_change",
  "w5score",
  "w3score",
#  "m5score",
#  "m3score",
#  "mwdif_5score",
#  "mwdif_3score",
  "ss5usage",
  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

mlDfNoMotifs <- as.data.frame(mcols(gr))[, c(
  "hek_ws",
  "hek_wu",
  "mutation_base_change",
#  "w5score",
#  "w3score",
  "m5score",
  "m3score",
  "mwdif_5score",
  "mwdif_3score",
  "ss5usage",
  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

mlDfNoExons <- as.data.frame(mcols(gr))[, c(
#  "hek_ws",
#  "hek_wu",
  "mutation_base_change",
  "w5score",
  "w3score",
  "m5score",
  "m3score",
  "mwdif_5score",
  "mwdif_3score",
  "ss5usage",
  "ss3usage"
#  "chasin_ese_density",
#  "chasin_ess_density"
)]

mlDfNoTranscripts <- as.data.frame(mcols(gr))[, c(
  "hek_ws",
  "hek_wu",
  "mutation_base_change",
  "w5score",
  "w3score",
  "m5score",
  "m3score",
  "mwdif_5score",
  "mwdif_3score",
#  "ss5usage",
#  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

featureSets <- list(mlDf1, mlDfNoMuts, mlDfNoMotifs, mlDfNoExons, mlDfNoTranscripts)

runRF <- function(df) {
  train <- df[trainIndices, ]
  test <- df[testIndices, ]

  model <- randomForest(train, response, importance=T)
  predNums <- predict(model, test, type="prob")

  realityLabels <- affects_splicing[testIndices]
  rocrPred <- prediction(predNums[, 2], realityLabels)
  auc <- performance(rocrPred, "auc")
  myRoc <- performance(rocrPred, "tpr", "fpr")
  return(data.frame(tpr=myRoc@x.values[[1]], fpr=myRoc@y.values[[1]], auc=auc@y.values[[1]]))
}

rocs <- lapply(featureSets, runRF)

ggplot() +
  geom_line(data=rocs[[1]], aes(tpr, fpr),
            size=2, alpha=0.7) +
  labs(x="False Positive Rate (1-Specificity)",
       y="True Positive Rate (Sensitivity)",
       color="Feature Sets") +
  theme(text=element_text(size=14)) +
  geom_abline(slope=1, linetype="dotted") +
  coord_fixed()
ggsave(file.path("..", "plots", "mapsy_roc_single.pdf"))

## Barplot of AUCs
featureSet <- c(
  "All Features",
  "No Mutation Features",
  "No Motif Features",
  "No Exon Features",
  "No Transcript Features"
)
featureSet <- factor(featureSet, levels=featureSet)
barDf <- data.frame(
  FeatureSet=featureSet,
  AUC=c(
    rocs[[1]]$auc[[1]],
    rocs[[2]]$auc[[1]],
    rocs[[3]]$auc[[1]],
    rocs[[4]]$auc[[1]],
    rocs[[5]]$auc[[1]]
  ))
cols <- c("#585858", brewer.pal(4, "Set2"))
ggplot(barDf, aes(FeatureSet, AUC)) +
  geom_bar(stat="identity", fill=cols) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        text=element_text(size=18)) +
  coord_cartesian(ylim=c(0.6, 0.865))
ggsave(file.path("..", "plots", "mapsy_auc_feature_levels.pdf"))
