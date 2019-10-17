#!/usr/bin/env Rscript
## Writes genomic exon features that are saved as RDS files into an SQL table
suppressPackageStartupMessages({
  library(getopt)
  library(RSQLite)
  library(sqldf)
  library(GenomicRanges)
})

opt <- getopt(matrix(c("outfile", "o", 1, "character",
                       "overwrite", "w", 0, "logical"),
                     byrow=T, ncol=4))
if (is.null(opt$outfile)) opt$outfile <- file.path("..", "data", "features.db")
if (is.null(opt$overwrite)) opt$overwrite <- F

## Load the saved features
maxent <- readRDS(file.path("..", "data", "wt_maxent.rds"))
ss_usage <- readRDS(file.path("..", "data", "ss_usage.rds"))
chasin <- as.data.frame(mcols(readRDS(file.path("..", "data", "chasin_cis_element_densities.rds"))))
avg_ei_scores <- as.data.frame(mcols(readRDS(file.path("..", "data", "exonic_avg_ei.rds"))))

## Join the features into a single data frame to be written into the SQL database
maxent_df <- merge(maxent$ss5, maxent$ss3, by="exon_id", all=T)
colnames(maxent_df) <- c("exon_id", "ss5seq", "ss5score", "ss3seq", "ss3score")
ss_usage_df <- as.data.frame(merge(mcols(ss_usage$ss5), mcols(ss_usage$ss3), by="exon_id", all=T))
colnames(ss_usage_df) <- c("exon_id", "ss5usage", "ss3usage")
#df3 <- merge(ss_usage_df, chasin, by="exon_id", all=T)

#DF <- merge(maxent_df, df3, by="exon_id", all=T)

allCols <- paste(c("maxent_df.exon_id",
                   "ss5seq", "ss5score", "ss3seq", "ss3score",
                   "ss5usage", "ss3usage",
                   "chasin_ese_density", "chasin_ess_density",
                   "ei"), collapse=", ")
## sqldf cannot do OUTER JOINs, but the maxent_df should have all exon_ids in it, so LEFT JOIN should give the same results
query <- sprintf(paste("SELECT %s FROM maxent_df",
                       "LEFT JOIN ss_usage_df",
                       "ON maxent_df.exon_id=ss_usage_df.exon_id",
                       "LEFT JOIN chasin",
                       "ON maxent_df.exon_id=chasin.exon_id",
                       "LEFT JOIN avg_ei_scores",
                       "ON maxent_df.exon_id=avg_ei_scores.exon_id;"), allCols)

DF <- sqldf(query)


## Write the table to the db
sqlTable <- "exons"
db <- dbConnect(SQLite(), opt$outfile)

if (opt$overwrite) {
  if (dbExistsTable(db, sqlTable)) dbRemoveTable(db, sqlTable)

  s <- sprintf("CREATE TABLE %s(%s, PRIMARY KEY(%s));", sqlTable,
                  paste(names(DF), collapse=", "),
                  names(DF)[1])
  dbSendStatement(db, s)
}

dbAppendTable(db, sqlTable, DF)
dbDisconnect(db)
