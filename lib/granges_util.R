suppressPackageStartupMessages({
  library(GenomicRanges)
})

## Grow/shrink both sides of a GRanges object by specified amounts (can be negative)
## If the range length will become < 0, return the original GRange and add an mcol column
## to mark which ranges were unaltered
## gr: A GRanges object
## left: Numeric to describe how much to extend the left side further (strand-aware)
## right: Numeric to describe how much to extend the right side further (strand-aware)
AlterBothSides <- function(gr, left, right) {
  ## Doesn't work for negatives right now
  ## What about 0's?
  lgr <- flank(gr, left)
  if (left < 0) lgr <- flank(lgr, 1, start=F)
  newWidth <- width(gr) + left + right
  return(resize(lgr, newWidth))
}

AddSpliceSites <- function(gr) AlterBothSides(gr, 20, 6)
