#!/usr/bin/env Rscript
## Loads MaPSy master file, calculates additional features, and saves it out as a GRanges object

## Additional features:
## 1) is the wt or mutant paired based on RNA-Fold data
## 2) Avg delta EI score for all hexamers that include the mutation position
## 3) Avg delta Shendure A5SS and A3SS scores (exonic only)

## Use ViennaRNA RNA-Fold predictions to determine if a mutation position is pair/unpaired in
## both the mutant and the wild-type

suppressPackageStartupMessages({
  source(file.path("..", "lib", "feature_gen.R"), chdir=T)
})

## Load the data
df <- read.csv(file.path("..", "data", "pMaster_020615_1.csv"), stringsAsFactors=F)
df <- df[grepl("^HGMD", df$classifier), ]
mutLocs <- FindMutPos(df$viv_seq_w, df$viv_seq_m)


## RNA-Fold pairings
isWtPaired <- substr(df$viv_vienna_w, mutLocs, mutLocs) != "."
isMutPaired <- substr(df$viv_vienna_m, mutLocs, mutLocs) != "."


## Avg delta EI scores
wtHexamers <- lapply(seq_along(df$viv_seq_w), function(i) {
  loc <- mutLocs[[i]]
  up <- loc - 5
  down <- loc + 5
  substring(df[i, "viv_seq_w"], up:loc, loc:down)
})

mutHexamers <- lapply(seq_along(df$viv_seq_m), function(i) {
  loc <- mutLocs[[i]]
  up <- loc - 5
  down <- loc + 5
  substring(df[i, "viv_seq_m"], up:loc, loc:down)
})

hexamerEIScores <- read.table(file.path("..", "data", "chasin_ESS_ESE_need.txt"),
                              row.names=1, col.names=c("seq", "score", "type"), sep="\t", na.strings="none")
hexamerEIScores[is.na(hexamerEIScores)] <- 0

avgDeltaEIScores <- sapply(seq_along(wtHexamers), function(i) {
  wt <- hexamerEIScores[wtHexamers[[i]], "score"]
  mut <- hexamerEIScores[mutHexamers[[i]], "score"]
  mean(mut-wt)
})


## Shendure scores
## There are exonic and intronic scores
## The exonic A3SS is "R2" and the exonic A5SS is "R1"
shendureA3SS <- read.delim(file.path("..", "data", "mean_effects_sdpos_A3SS_R2.txt"), row.names="motif", stringsAsFactors=F)
shendureA5SS <- read.delim(file.path("..", "data", "mean_effects_sdpos_A5SS_R1.txt"), row.names="motif", stringsAsFactors=F)

avgDeltaA3SSScores <- sapply(seq_along(wtHexamers), function(i) {
  wt <- shendureA3SS[wtHexamers[[i]], "A3SS_R2_mean_effects"]
  mut <- shendureA3SS[mutHexamers[[i]], "A3SS_R2_mean_effects"]
  mean(mut-wt)
})
avgDeltaA5SSScores <- sapply(seq_along(wtHexamers), function(i) {
  wt <- shendureA5SS[wtHexamers[[i]], "A5SS_R1_mean_effects"]
  mut <- shendureA5SS[mutHexamers[[i]], "A5SS_R1_mean_effects"]
  mean(mut-wt)
})


## GRanges object of alternative exons in MaPSy
gr <- GRanges(seqnames=df$chrom,
              ranges=IRanges(start=df$exst+1, end=df$exen),
              strand=df$strand,
              seqinfo=Seqinfo(genome="hg19"))
gr$Idw <- df$Idw
gr$Idm <- df$Idm
gr$isWtPaired <- isWtPaired
gr$isMutPaired <- isMutPaired
gr$avgDeltaEIScores <- avgDeltaEIScores
gr$avgDeltaA3SSScores <- avgDeltaA3SSScores
gr$avgDeltaA5SSScores <- avgDeltaA5SSScores
gr$geneLength <- df$genelength
gr$nExons <- df$gene_Nex
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

saveRDS(gr, file.path("..", "data", "mapsy_features_gr.rds"))
