## This file has functions to query an SQL database of exonic features across the genome

library(RSQLite)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

dbSource <- file.path("..", "data", "features.db")

## Helper function to convert a vector into a string separated by commas
mkString <- function(vec) { return(paste(vec, collapse=", ")) }

## Adds new metadata columns in bulk to a GRanges object
AddMcols <- function(gr, data, by="exon_id") {
  newMcols <- merge(mcols(gr), data, by=by, all.x=T)
  mcols(gr) <- newMcols
  gr <- gr[complete.cases(mcols(gr))]
  return(gr)
}

## Attempts to add an exon_id column to a GRanges object
AddExonId <- function(gr) {
  exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
  ovs <- mergeByOverlaps(gr, exs, type="equal")
  gr2 <- ovs$gr
  gr2$exon_id <- ovs$exon_id
  return(gr2)
}


## Gets 5'SS Maxent scores by exon_id
## exon_ids: Numeric vector of exon ids
QuerySS5Scores <- function(exon_ids) {
  con <- dbConnect(SQLite(), dbSource)
  query <- sprintf("SELECT exon_id, ss5score FROM exons WHERE exon_id IN (%s);",
                   mkString(exon_ids))
  result <- dbGetQuery(con, query)
  dbDisconnect(con)
  return(result)
}

## Gets 3'SS Maxent scores by exon_id
## exon_ids: Numeric vector of exon ids
QuerySS3Scores <- function(exon_ids) {
  con <- dbConnect(SQLite(), dbSource)
  query <- sprintf("SELECT exon_id, ss3score FROM exons WHERE exon_id IN (%s);",
                   mkString(exon_ids))
  result <- dbGetQuery(con, query)
  dbDisconnect(con)
  return(result)
}

## Get all exon info for those that have both 5' and 3'SS usage
## as.GRanges: Logical determining if the result should be returned as a GRanges object
QueryExonsWithSSUsage <- function(as.GRanges=F) {
  con <- dbConnect(SQLite(), dbSource)
  query <- "SELECT * FROM exons WHERE ss5usage NOT NULL AND ss3usage NOT NULL;"
  result <- dbGetQuery(con, query)
  dbDisconnect(con)
  if (as.GRanges) {
    exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
    filteredExs <- exs[exs$exon_id %in% result$exon_id]
    mcols(filteredExs) <- merge(mcols(filteredExs), result)
    result <- filteredExs
  }

  return(result)
}
