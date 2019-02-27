#!/usr/bin/env Rscript
## Train GBM model on MaPSy alleles and store it

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"))
  library(sqldf)
  library(gbm)
  library(parallel)
  library(ROCR)
})

nCores <- detectCores()

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

#allelic_skew <- AllelicSkew(gr$hek_ws, gr$hek_wu, gr$hek_ms, gr$hek_mu)
#affects_splicing <- as.factor(allelic_skew >= log2(1.5)) # p-value?
#gr$affects_splicing <- affects_splicing
#responseGBM <- as.integer(affects_splicing) - 1
splitIdx <- floor(0.7 * nRows)
trainIndices <- 1:splitIdx
testIndices <- (splitIdx+1):nRows

allelic_skew <- AllelicSkew(gr$hek_ws, gr$hek_wu, gr$hek_ms, gr$hek_mu)
affects_splicing <- as.factor(allelic_skew >= log2(1.5)) # p-value?
gr$affects_splicing <- affects_splicing
response <- affects_splicing[trainIndices]
responseGBM <- as.integer(response) - 1

mlDf <- as.data.frame(mcols(gr))[, c(
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
  "ss3usage",
  "chasin_ese_density",
  "chasin_ess_density"
)]

gbmDf <- cbind(responseGBM, mlDf[trainIndices, ])
test <- mlDf[testIndices, ]

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
predNums <- predict(model, test, n.trees=best.iter, type="response")

realityLabels <- affects_splicing[testIndices]
rocrPred <- prediction(predNums, realityLabels)
auc <- performance(rocrPred, "auc")
myRoc <- performance(rocrPred, "tpr", "fpr")

saveRDS(list(model=model, n.trees=best.iter), file.path("..", "data", "ml_model_no_hek.rds"))
