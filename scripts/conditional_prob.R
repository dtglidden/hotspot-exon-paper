#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(getopt)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(RColorBrewer)
  library(tidyverse)
})

opt <- getopt(matrix(c("test", "t", 0, "logical"),
                     byrow=T, ncol=4))
if (is.null(opt$test)) opt$test <- F

set.seed(1)

df <- if (opt$test) {
  read.delim(file.path("..", "example_data", "conditional_prob_test.txt"))
} else {
  read.delim(file.path("..", "data", "input_ML_complete_table.txt"))
}

## Only consider exons that have multiple MaPSy mutations in them
df <- df[duplicated(df$ref) | duplicated(df$ref, fromLast=T), ]

## Remove species with more than 1 mutation in the individual construct
df <- df[!grepl("*", df$mut_pos_construct, fixed=T), ]

df2 <- df %>%
  mutate(rand_vivo_ratio=sample(vivo_ratio)) %>%
  group_by(chrom_exon, start_exon, end_exon) %>%
  summarize(std_vivo_ratio=sd(vivo_ratio),
            std_rand_vivo_ratio=sd(rand_vivo_ratio)) %>%
  drop_na

nRow <- nrow(df2)
df3 <- data.frame(x=c(df2$std_vivo_ratio, df2$std_rand_vivo_ratio),
                  group=factor(rep(c("Unshuffled", "Shuffled"), rep(nRow, 2))))

ggdensity(df3, "x", color="group", fill="group", alpha=0.5, add="mean") +
  labs(x="Std(Allelic Ratio)", y="Density", fill="") +
  xlim(c(0, 8)) +
  theme(text=element_text(size=28)) +
  scale_color_discrete(guide=F)
ggsave(file.path("..", "plots", "conditional_prob_std_by_exon_rand.pdf"))

## Wilcoxon Signed Rank Test: p-value < 2.2e-16
wilcox.test(x ~ group, df3, paired=T, conf.int=T)


## Conditional probability analysis
## Sample 2 mutations from each exon and determine the probability that the second mutation affect splicing
## given the effect of the first mutation
df3 <- df %>%
  group_by(ref) %>%
  summarize(samples=list(sample(vivo_ratio))) %>%
  mutate(first_affects_splicing=sapply(samples, function(s) abs(s[[1]]) >= log2(1.5)),
         second_affects_splicing=sapply(samples, function(s) abs(s[[2]]) >= log2(1.5))) %>%
  group_by(first_affects_splicing) %>%
  summarize(prob_second_affects_splicing=sum(second_affects_splicing)/length(second_affects_splicing)) %>%
  mutate(first_affects_splicing=fct_recode(as.factor(first_affects_splicing),
                                           "Affects splicing"="TRUE", "Does not affect splicing"="FALSE"))

ggbarplot(df3, "first_affects_splicing", "prob_second_affects_splicing",
          fill="first_affects_splicing",
          order=c("Affects splicing", "Does not affect splicing"),
          xlab="First mutation",
          ylab="Likelihood second\nmutation affects splicing") +
  theme(text=element_text(size=20),
        axis.text.y=element_text(size=24),
        axis.title=element_text(size=28),
        legend.position="none")
ggsave(file.path("..", "plots", "conditional_prob_first_affects_second.pdf"))
