#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
})

set.seed(1)

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
samps <- lapply(1:length(dfByWt), function(i) dfByWt[[i]][sample(nrow(dfByWt[[i]]), 2), ])

firsts <- sapply(samps, function(x) x[1, "vivo_ratio"])

ones   <- sapply(Filter(function(x) x[1, "bin"]==1, samps), function(x) x[2, "vivo_ratio"])
twos   <- sapply(Filter(function(x) x[1, "bin"]==2, samps), function(x) x[2, "vivo_ratio"])
threes <- sapply(Filter(function(x) x[1, "bin"]==3, samps), function(x) x[2, "vivo_ratio"])
fours  <- sapply(Filter(function(x) x[1, "bin"]==4, samps), function(x) x[2, "vivo_ratio"])

allSnv1 <- sapply(samps, function(x) x[1, "vivo_ratio"])
allSnv1Sum <- summary(allSnv1)
allSnv <- c(allSnv1, sapply(samps, function(x) x[2, "vivo_ratio"]))
alls <- c(ones, twos, threes, fours, allSnv1)
xAxis <- seq(floor(min(alls)), ceiling(max(alls)), by=0.5)

colorPal <- brewer.pal(4, "Set3")
colors <- sapply(xAxis, function(x) {
  if (x < allSnv1Sum["1st Qu."]) colorPal[[1]]
  else if (x < allSnv1Sum["Median"]) colorPal[[2]]
  else if (x < allSnv1Sum["3rd Qu."]) colorPal[[3]]
  else colorPal[[4]]
})

myXlim <- c(-5, 5)
myYlim <- c(0, 400)

SE <- function(x) sd(x) / sqrt(length(x))
barNames <- 1:4
barHeights <- c(mean(ones), mean(twos), mean(threes), mean(fours))
barSE <- sapply(list(ones, twos, threes, fours), SE)
barDf <- data.frame(quartile=barNames,
                    ratio=barHeights,
                    se=barSE)
ggplot(barDf, aes(quartile, ratio)) +
  geom_col(color="black") +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                width=0.2) +
  labs(x="Quartile of SNV 1",
       y="Average M/W Ratio of SNV 2") +
  theme_classic() +
  theme(text=element_text(size=24))
ggsave(file.path("..", "plots", "conditional_prob_bar.pdf"))
p1 <- ggplot(data.frame(ones), aes(ones)) +
  geom_histogram(data=data.frame(allSnv1), mapping=aes(allSnv1), breaks=xAxis, color=NA, fill="gray35") +
  geom_histogram(breaks=xAxis, color=NA, fill="#8DD3C7") +
  ylab("Frequency") +
  xlab("M/W Ratio") +
  xlim(myXlim) +
  ylim(myYlim) +
  theme(text=element_text(size=8, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p2 <- ggplot(data.frame(twos), aes(twos)) +
  geom_histogram(data=data.frame(allSnv1), mapping=aes(allSnv1), breaks=xAxis, color=NA, fill="gray35") +
  geom_histogram(breaks=xAxis, color=NA, fill="#FFFFB3") +
  ylab("Frequency") +
  xlab("M/W Ratio") +
  xlim(myXlim) +
  ylim(myYlim) +
  theme(text=element_text(size=8, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p3 <- ggplot(data.frame(threes), aes(threes)) +
  geom_histogram(data=data.frame(allSnv1), mapping=aes(allSnv1), breaks=xAxis, color=NA, fill="gray35") +
  geom_histogram(breaks=xAxis, color=NA, fill="#BEBADA") +
  ylab("Frequency") +
  xlab("M/W Ratio") +
  xlim(myXlim) +
  ylim(myYlim) +
  theme(text=element_text(size=8, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))
p4 <- ggplot(data.frame(fours), aes(fours)) +
  geom_histogram(data=data.frame(allSnv1), mapping=aes(allSnv1), breaks=xAxis, color=NA, fill="gray35") +
  geom_histogram(breaks=xAxis, color=NA, fill="#FB8072") +
  ylab("Frequency") +
  xlab("M/W Ratio") +
  xlim(myXlim) +
  ylim(myYlim) +
  theme(text=element_text(size=8, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA))

pdf(file.path("..", "plots", "conditional_prob_ggplot_hist.pdf"), height=2)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), size="last"))
dev.off()



## Plot the overall distribution of all SNVs
ggplot(data.frame(allSnv), aes(allSnv)) +
  geom_histogram(breaks=xAxis, color=NA) +
  ylab("Frequency") +
  xlab("M/W Ratio") +
#  ylim(myYlim) +
  theme_classic() +
  theme(text=element_text(size=24))
ggsave(file.path("..", "plots", "conditional_prob_all_snv.pdf"))
