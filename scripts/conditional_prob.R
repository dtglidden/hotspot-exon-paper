#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(RColorBrewer)
  library(tidyverse)
})

set.seed(1)

df <- read.delim(file.path("..", "data", "input_ML_complete_table.txt"))

## Only consider exons that have multiple MaPSy mutations in them
df <- df[duplicated(df$ref) | duplicated(df$ref, fromLast=T), ]

## Remove species with more than 1 mutation in the individual construct
df <- df[!grepl("*", df$mut_pos_construct, fixed=T), ]

df2 <- df %>%
  group_by(chrom_exon, start_exon, end_exon) %>%
  summarize(std_vivo_ratio=sd(vivo_ratio)) %>%
  filter(!is.na(std_vivo_ratio))

ggdensity(df2, "std_vivo_ratio", fill="gray") +
  geom_vline(aes(xintercept=sd(df$vivo_ratio)),
             linetype="dashed",
             color="red") +
  labs(x="Std(Allelic Ratio)", y="Density") +
  theme(text=element_text(size=24))
ggsave(file.path("..", "plots", "conditional_prob_std_by_exon.pdf"))
