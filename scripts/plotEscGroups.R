#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(getopt)
  library(ggplot2)
})

opt <- getopt(matrix(c("test", "t", 0, "logical"),
                     byrow=T, ncol=4))
if (is.null(opt$test)) opt$test <- F

df <- if (opt$test) {
  read.delim(file.path("..", "example_data", "esco_test.txt"))
} else {
  read.delim(file.path("..", "data", "see-saw.txt"))
}

ggplot(df, aes(v_exon_eff, skip.incl, color=group)) +
  geom_point(size=2.5) +
  labs(x="Splicing Efficiency of flanking exon",
       y="Skipping:Inclusion ratio of test exon",
       color="") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        text=element_text(size=18)) +
  stat_smooth(method="lm", se=F)
ggsave(file.path("..", "plots", "see-saw_with_groups.pdf"), width=6, height=5, useDingbats=F)
