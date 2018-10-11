#!/usr/bin/env Rscript
## Writes genomic exon features that are saved as RDS files into an SQL table
library(getopt)
library(RSQLite)
library(GenomicRanges, quietly=T)

opt <- getopt(matrix(c("outfile", "o", 1, "character",
                       "overwrite", "w", 0, "logical"),
                     byrow=T, ncol=4))
if (is.null(opt$outfile)) opt$outfile <- file.path("..", "data", "features.db")
if (is.null(opt$overwrite)) opt$overwrite <- F

## Load the saved features
maxent <- readRDS(file.path("..", "data", "wt_maxent.rds"))
ss_usage <- readRDS(file.path("..", "data", "ss_usage.rds"))
chasin <- as.data.frame(mcols(readRDS(file.path("..", "data", "chasin_cis_element_densities.rds"))))

## Join the features into a single data frame to be written into the SQL database
df <- merge(maxent$ss5, maxent$ss3, by="exon_id", all=T)
colnames(df) <- c("exon_id", "ss5seq", "ss5score", "ss3seq", "ss3score")
df2 <- as.data.frame(merge(mcols(ss_usage$ss5), mcols(ss_usage$ss3), by="exon_id", all=T))
colnames(df2) <- c("exon_id", "ss5usage", "ss3usage")
df3 <- merge(df2, chasin, by="exon_id", all=T)

DF <- merge(df, df3, by="exon_id", all=T)

## Write the table to the db
sqlTable <- "exons"
db <- dbConnect(SQLite(), opt$outfile)

if (opt$overwrite) {
  if (dbExistsTable(db, sqlTable)) dbRemoveTable(db, sqlTable)

  s <- sprintf("CREATE TABLE %s(%s, PRIMARY KEY(%s))", sqlTable,
                  paste(names(DF), collapse=", "),
                  names(DF)[1])
  dbSendStatement(db, s)
}

dbAppendTable(db, sqlTable, DF)
dbDisconnect(db)
