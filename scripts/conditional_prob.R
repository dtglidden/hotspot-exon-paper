#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

df <- read.delim(file.path("..", "data", "input_ML_complete_table.txt"))

## Only consider exons that have multiple MaPSy mutations in them
df <- df[duplicated(df$ref) | duplicated(df$ref, fromLast=T), ]

## Remove species with more than 1 mutation in the individual construct
df <- df[!grepl("*", df$mut_pos_construct, fixed=T), ]

## Creates quartiles for mw ratio
quarts <- summary(df$vivo_ratio)

between <- function(e, x, y) {
  return(e >= x & e < y)
}

## Assigns a quartile for all vivo_ratios
df$bin <- sapply(df$vivo_ratio, function(e) {
  if (between(e, quarts["Min."], quarts["1st Qu."]))
    return(1)
  else if (between(e, quarts["1st Qu."], quarts["Median"]))
    return(2)
  else if (between(e, quarts["Median"], quarts["3rd Qu."]))
    return(3)
  else if (between(e, quarts["3rd Qu."], quarts["Max."]))
    return(4)
  else ## This value is the max value
    return(4)
})

## Group the dataframe by the WT species (creates a list of dataframes)
dfByWt <- lapply(unique(df$ref), function(x) df[df$ref == x, c("alt", "vivo_ratio", "bin")])
names(dfByWt) <- unique(df$ref)

## Randomly sample 2 items from each group
samps <- lapply(1:length(dfByWt), function(i) sample(dfByWt[[i]]$bin, 2))

breaks <- c(0.5, 1.5, 2.5, 3.5, 4.5)
ones <- hist(sapply(Filter(function(x) x[[1]]==1, samps), function(x) x[[2]]),
             breaks=breaks,
             plot=F)
ones <- ones$counts / sum(ones$counts)
twos <- hist(sapply(Filter(function(x) x[[1]]==2, samps), function(x) x[[2]]),
             breaks=breaks,
             plot=F)
twos <- twos$counts / sum(twos$counts)
threes <- hist(sapply(Filter(function(x) x[[1]]==3, samps), function(x) x[[2]]),
               breaks=breaks,
               plot=F)
threes <- threes$counts / sum(threes$counts)
fours <- hist(sapply(Filter(function(x) x[[1]]==4, samps), function(x) x[[2]]),
              breaks=breaks,
              plot=F)
fours <- fours$counts / sum(fours$counts)

p1 <- ggplot(data.frame(x=1:4, ones), aes(x, ones)) +
  geom_col(width=1, color="black") +
  ylab("Probability") +
  xlab("Quartile of SNV 2") +
  ylim(NA, 0.5) +
  theme(text=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p2 <- ggplot(data.frame(x=1:4, twos), aes(x, twos)) +
  geom_col(width=1, color="black") +
  ylab("Probability") +
  xlab("Quartile of SNV 2") +
  ylim(NA, 0.5) +
  theme(text=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p3 <- ggplot(data.frame(x=1:4, threes), aes(x, threes)) +
  geom_col(width=1, color="black") +
  ylab("Probability") +
  xlab("Quartile of SNV 2") +
  ylim(NA, 0.5) +
  theme(text=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p4 <- ggplot(data.frame(x=1:4, fours), aes(x, fours)) +
  geom_col(width=1, color="black") +
  ylab("Probability") +
  xlab("Quartile of SNV 2") +
  ylim(NA, 0.5) +
  theme(text=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))

pdf(file.path("..", "plots", "conditional_prob_ggplot_hist.pdf"), height=2)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), size="last"))
dev.off()
