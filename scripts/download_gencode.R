#!/usr/bin/env Rscript

# Downloads the Genome reference from Gencode in order to calculate genomic features

hg38 <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
hg19 <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz"

download.file(hg38, file.path("..", "data", "hg38.gtf.gz"))
download.file(hg19, file.path("..", "data", "hg19.gtf.gz"))
