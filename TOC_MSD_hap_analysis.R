# TOC_MSD: 1A_block_20384 与 2D_block_177 单倍型 P7 ANOVA + 育种历程占比
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

## ---- 读基因型（0/0, 1/1）----
g <- read.delim("TOC_MSD.genotype", header=TRUE, comment.char="", check.names=FALSE)
gt <- as.matrix(g[, -(1:9)])
storage.mode(gt) <- "character"
samp_ids <- colnames(gt)
hap_all <- sub("_number\\.txt-.*", "", g$ID)
hap_id  <- sub(".*_number\\.txt-", "", g$ID)
blocks  <- unique(hap_all)
stopifnot(length(blocks) == 2)

hap_list <- list()
for (bl in blocks) {
  h <- rep(NA_character_, length(samp_ids))
  for (k in which(hap_all == bl)) h[gt[k, ] == "1/1"] <- hap_id[k]
  nh <- sum(hap_all == bl)
  hap_list[[bl]] <- data.frame(ID=samp_ids,
                               hap=factor(h, levels=paste0("Hap", 1:nh)))
}
cat("单倍型分布:\n")
for (bl in blocks) { cat(bl, ":\n"); print(table(hap_list[[bl]]$hap, useNA="ifany")) }

## ---- 表型 P7 ----
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe$V1 <- as.character(phe$V1)
phe <- phe[, c(1, 8)]
names(phe) <- c("ID", "P7")

## ---- 分组 ----
grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID","group"), comment.char="#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace","Introduced","Pre-1960","1961-1980","1981-2000","2001-2020")
grp$group <- factor(grp$group, levels=grp_order)

pal <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69")

for (bl in blocks) {
  blname <- sub("^X", "", bl)
  d1 <- merge(hap_list[[bl]], phe, by="ID")
  d1 <- d1[!is.na(d1$P7) & !is.na(d1$hap), ]
  cat("\n===", blname, "P7 有效材料:", nrow(d1), "===\n")

  ## ---- 任务1: ANOVA + 小提琴图 ----
  fit <- aov(P7 ~ hap, data=d1)
  print(summary(fit))
  tuk <- agricolae::HSD.test(fit, "hap", console=FALSE)
  ltr <- tolower(tuk$groups$groups); names(ltr) <- rownames(tuk$groups)
  cat("Tukey HSD:\n"); print(tuk$groups)

  ms1 <- d1 %>% group_by(hap) %>%
    summarise(n=n(), mean=mean(P7), sd=sd(P7), median=median(P7),
              sem=sd(P7)/sqrt(n()), .groups="drop") %>% as.data.frame()
  ms1$letter <- ltr[as.character(ms1$hap)]
  ms1$mean_sd <- sprintf("%.3f±%.3f", ms1$mean, ms1$sd)
  ms1$anova_P <- signif(summary(fit)[[1]][1,5], 3)
  write.table(ms1[, c("hap","n","mean","sd","median","sem","letter","mean_sd","anova_P")],
              paste0("TOC_MSD_", blname, "_P7_anova_means.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  nh <- nlevels(d1$hap)
  p1 <- ggplot(d1, aes(x=hap, y=P7, fill=hap)) +
    geom_violin(trim=FALSE, scale="width", linewidth=0.3, alpha=0.8,
                colour="grey25") +
    geom_boxplot(width=0.12, outlier.shape=NA, fill="white",
                 alpha=0.9, linewidth=0.3) +
    geom_text(data=ms1, aes(x=hap, y=max(d1$P7)*1.06, label=letter),
              size=5, fontface="bold", inherit.aes=FALSE) +
    scale_fill_manual(values=pal[1:nh]) +
    labs(x=paste0("Haplotype of ", blname), y="P7 phenotype",
         title=paste0("P7 phenotype by haplotype of ", blname,
                      " (ANOVA + Tukey HSD)")) +
    annotate("text", x=0.6, y=max(d1$P7)*1.10, hjust=0, size=3.2,
             label=paste0("ANOVA P = ", signif(summary(fit)[[1]][1,5], 3))) +
    coord_cartesian(ylim=c(min(d1$P7), max(d1$P7)*1.15)) +
    theme_bw() + theme(legend.position="none", panel.grid.minor=element_blank())
  ggsave(paste0("TOC_MSD_", blname, "_P7_violin.pdf"), p1, width=7.5, height=6)

  ## ---- 任务2: 分组占比堆积图 ----
  d2 <- merge(hap_list[[bl]], grp, by="ID") %>%
    filter(!is.na(group), !is.na(hap))
  prop2 <- d2 %>% group_by(group, hap) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(group) %>% mutate(total=sum(n), prop=n/total) %>% as.data.frame()
  prop2$prop_pct <- sprintf("%.1f%%", prop2$prop*100)
  write.table(prop2[, c("group","hap","n","total","prop","prop_pct")],
              paste0("TOC_MSD_", blname, "_group_prop.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  w2 <- prop2 %>% select(group, hap, prop_pct) %>%
    pivot_wider(names_from=hap, values_from=prop_pct, values_fill="0.0%") %>%
    as.data.frame()
  write.table(w2, paste0("TOC_MSD_", blname, "_group_prop_wide.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  p2 <- ggplot(prop2, aes(x=group, y=prop, fill=hap)) +
    geom_col(position="stack", colour="black", linewidth=0.2, width=0.75) +
    geom_text(aes(label=ifelse(prop>=0.05, paste0(round(prop*100), "%"), "")),
              position=position_stack(vjust=0.5), size=3) +
    scale_fill_manual(values=pal[1:nh]) +
    scale_x_discrete(limits=grp_order) +
    scale_y_continuous(labels=scales::percent_format()) +
    labs(x="Breeding era", y="Proportion", fill="Haplotype",
         title=paste0("Haplotype composition of ", blname, " by breeding era")) +
    theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1),
                       panel.grid.minor=element_blank())
  ggsave(paste0("TOC_MSD_", blname, "_group_stack.pdf"), p2, width=8, height=6)
  print(w2)
}
cat("DONE\n")