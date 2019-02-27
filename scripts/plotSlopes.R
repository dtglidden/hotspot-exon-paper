#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

Slope <- function(data) {
  return(coefficients(lm(Y ~ X, data=data))["X"])
}

df <- read.delim(file.path("..", "data", "see-saw.txt"))
dfSummary <- df %>% group_by(f_exon) %>%
  summarize(slope=Slope(data.frame(X=v_exon_eff, Y=skip.incl)),
            f_exon_eff=mean(f_exon_eff)) %>%
  arrange(f_exon_eff)

ggplot(dfSummary, aes(f_exon_eff, slope)) +
  geom_point(size=2.5) +
  labs(x="Splicing Efficiency of test exon",
       y="Slope of linear regression model") +
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(margin=margin(t=0, r=5, b=0, l=0)),
        text=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
ggsave(file.path("..", "plots", "see-saw_slopes.pdf"), width=3, height=3, useDingbats=F)
