## This file has functions to query an SQL database of exonic features across the genome

library(RSQLite)

dbSource <- file.path("..", "data", "features.db")

## Helper function to convert a vector into a string separated by commas
mkString <- function(vec) { return(paste(vec, collapse=", ")) }

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
