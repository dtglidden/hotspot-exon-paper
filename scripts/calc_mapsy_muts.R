#!/usr/bin/env Rscript
## Calculate mutation features for MaPSy alleles

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  library(sqldf)
  library(randomForest)
  library(gbm)
  library(ROCR)
  library(ggplot2)
  library(RColorBrewer)
  library(parallel)
})

nCores <- detectCores()
set.seed(1)

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

gr <- readRDS(file.path("..", "data", "mapsy_features_gr.rds"))
## Temp fix: change booleans to integers in order to train ML models
gr$isWtPaired <- as.integer(gr$isWtPaired)
gr$isMutPaired <- as.integer(gr$isMutPaired)
## Have to merge dataframes via sqldf because the regular merge function affects the ordering.
## This way, the exon IDs are in the same order as those in 'gr'
grMcols <- as.data.frame(mcols(gr))
exomeMcols <- as.data.frame(mcols(exsWithSSUsage))
exomMcolNames <- paste(c("ss5seq", "ss5score", "ss3seq", "ss3score",
                         "ss5usage", "ss3usage",
                         "chasin_ese_density", "chasin_ess_density",
                         "ei"), collapse=", ")
mcols(gr) <- sqldf(sprintf(paste("SELECT grMcols.*, %s",
                                 "FROM grMcols LEFT JOIN exomeMcols",
                                 "ON grMcols.exon_id=exomeMcols.exon_id;"), exomMcolNames))
gr <- gr[complete.cases(mcols(gr))]

nRows <- length(gr)
splitIdx <- floor(0.7 * nRows)
trainIndices <- 1:splitIdx
testIndices <- (splitIdx+1):nRows

allelic_skew <- AllelicSkew(gr$hek_ws, gr$hek_wu, gr$hek_ms, gr$hek_mu)
affects_splicing <- as.factor(allelic_skew >= log2(1.5)) # p-value?
gr$affects_splicing <- affects_splicing
response <- affects_splicing[trainIndices]
responseGBM <- as.integer(response) - 1

featureDf <- read.csv(file.path("..", "data", "feature_table.csv"), stringsAsFactors=F)

## Variations of the data, so we can compare how some features affect AUC
mlDfAll <- as.data.frame(mcols(gr))[, featureDf$feature]
mlDfNoMuts <- as.data.frame(mcols(gr))[, featureDf$feature[featureDf$level != "mutation"]]
mlDfNoMotifs <- as.data.frame(mcols(gr))[, featureDf$feature[featureDf$level != "motif"]]
mlDfNoExons <- as.data.frame(mcols(gr))[, featureDf$feature[featureDf$level != "exon"]]
mlDfNoTranscripts <- as.data.frame(mcols(gr))[, featureDf$feature[featureDf$level != "transcript"]]

featureSets <- list(mlDfAll, mlDfNoMuts, mlDfNoMotifs, mlDfNoExons, mlDfNoTranscripts)

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

runGBM <- function(df) {
  train <- df[trainIndices, ]
  test <- df[testIndices, ]
  gbmDf <- cbind(responseGBM, train)

  #model <- gbm(responseGBM ~ ., data=gbmDf, distribution="bernoulli", cv.folds=5, n.cores=nCores)
  ## Most function arguments were taken from the example on the man page for 'gbm'
  model <- gbm(responseGBM ~ .,             # formula
               data=gbmDf,                  # dataset
               #var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
                                            # +1: monotone increase,
                                            #  0: no monotone restrictions
               distribution="bernoulli",    # see the help for other choices
               n.trees=1000,                # number of trees
               shrinkage=0.05,              # shrinkage or learning rate,
                                            # 0.001 to 0.1 usually work
               interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
               bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
               #train.fraction = 0.5,        # fraction of data for training,
                                            # first train.fraction*N used for training
               n.minobsinnode = 10,         # minimum total weight needed in each node
               cv.folds = 3,                # do 3-fold cross-validation
               keep.data=F,                 # don't keep a copy of the dataset with the object
               verbose=F,                   # don't print out progress
               n.cores=nCores)
  best.iter <- gbm.perf(model, method="cv", plot.it=F)
  ## 'response' type scales values to the 0-1 interval
  ## (i.e. the probability of affecting splicing or not)
  predNums <- predict(model, test, n.trees=best.iter, type="response")

  realityLabels <- affects_splicing[testIndices]
  rocrPred <- prediction(predNums, realityLabels)
  auc <- performance(rocrPred, "auc")
  myRoc <- performance(rocrPred, "tpr", "fpr")
  return(data.frame(tpr=myRoc@x.values[[1]], fpr=myRoc@y.values[[1]], auc=auc@y.values[[1]]))
}

rocs <- lapply(featureSets, runGBM)

ggplot() +
  geom_line(data=rocs[[1]], aes(tpr, fpr),
            size=2, alpha=0.7) +
  labs(x="False Positive Rate (1-Specificity)",
       y="True Positive Rate (Sensitivity)",
       color="Feature Sets") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        text=element_text(size=24)) +
  geom_abline(slope=1, linetype="dotted") +
  coord_fixed()
ggsave(file.path("..", "plots", "mapsy_roc_extra_features_gbm.pdf"))

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
  geom_bar(stat="identity", fill=cols, color="black") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        text=element_text(size=24)) +
  coord_cartesian(ylim=c(0.6, 0.9))
ggsave(file.path("..", "plots", "mapsy_auc_feature_levels_gbm_new_features.pdf"))

## Random Forest
rocs <- lapply(featureSets, runRF)

ggplot() +
  geom_line(data=rocs[[1]], aes(tpr, fpr),
            size=2, alpha=0.7) +
  labs(x="False Positive Rate (1-Specificity)",
       y="True Positive Rate (Sensitivity)",
       color="Feature Sets") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        text=element_text(size=24)) +
  geom_abline(slope=1, linetype="dotted") +
  coord_fixed()
ggsave(file.path("..", "plots", "mapsy_roc_extra_features_rf.pdf"))

## Barplot of AUCs

barDf <- data.frame(
  FeatureSet=featureSet,
  AUC=c(
    rocs[[1]]$auc[[1]],
    rocs[[2]]$auc[[1]],
    rocs[[3]]$auc[[1]],
    rocs[[4]]$auc[[1]],
    rocs[[5]]$auc[[1]]
  ))
ggplot(barDf, aes(FeatureSet, AUC)) +
  geom_bar(stat="identity", fill=cols, color="black") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        text=element_text(size=24)) +
  coord_cartesian(ylim=c(0.6, 0.9))
ggsave(file.path("..", "plots", "mapsy_auc_feature_levels_rf_new_features.pdf"))
