#!/usr/bin/env Rscript
## Train GBM model on MaPSy alleles and store it

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
  library(sqldf)
  library(dplyr)
  library(gbm)
  library(parallel)
  library(ROCR)
})

nCores <- detectCores()

set.seed(1337)

exsWithSSUsage <- QueryExonsWithSSUsage(as.GRanges=T)

gr <- readRDS(file.path("..", "data", "mapsy_features_gr.rds"))
## Temp fix: change booleans to integers in order to train ML models
gr$isWtPaired <- as.integer(gr$isWtPaired)
gr$isMutPaired <- as.integer(gr$isMutPaired)
## Have to merge dataframes via sqldf because the regular merge function affects the ordering.
## This way, the exon IDs are in the same order as those in 'gr'
grMcols <- as.data.frame(mcols(gr))
exomeMcols <- as.data.frame(mcols(exsWithSSUsage))
exomeMcolNames <- paste(c("ss5seq", "ss5score", "ss3seq", "ss3score",
                         "ss5usage", "ss3usage",
                         "chasin_ese_density", "chasin_ess_density",
                         "ei"), collapse=", ")
mcols(gr) <- sqldf(sprintf(paste("SELECT grMcols.*, %s",
                                 "FROM grMcols LEFT JOIN exomeMcols",
                                 "ON grMcols.exon_id=exomeMcols.exon_id;"), exomeMcolNames))
gr <- gr[complete.cases(mcols(gr))]

allelic_skew <- AllelicSkew(gr$hek_ws, gr$hek_wu, gr$hek_ms, gr$hek_mu)
affects_splicing <- as.factor(allelic_skew >= log2(1.5))
gr$affects_splicing <- affects_splicing
response <- affects_splicing
responseGBM <- as.integer(response) - 1

## Using MaPSy count severely limits the number of exons we can use for predictions,
## so we will exclude it to make predictions for all possible exonic variants.
## Splice site usage also restricts the number of exons to a lesser degree, but
## we will keep it in this model
featureDf <- read.csv(file.path("..", "data", "feature_table.csv"), stringsAsFactors=F) %>%
  filter(!feature %in% c("hek_wu", "hek_ws"))

mlDf <- cbind(responseGBM, as.data.frame(mcols(gr))[, featureDf$feature])

## Most function arguments were taken from the example on the man page for 'gbm'
model <- gbm(responseGBM ~ .,             # formula
              data=mlDf,                  # dataset
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

saveRDS(list(model=model, n.trees=best.iter), file.path("..", "data", "ml_model_no_hek.rds"))
