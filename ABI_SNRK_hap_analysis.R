# ABI_SNRK: 两个block单倍型 ANOVA + 育种历程占比
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

tr7 <- paste0("P", 1:7)

## ---- 读基因型（0/0 与 1/1）----
g <- read.delim("ABI_SNRK.genotype", header=TRUE, comment.char="", check.names=FALSE)
gt <- as.matrix(g[, -(1:9)])
storage.mode(gt) <- "character"
hap_all <- sub("_number\\.txt-.*", "", g$ID)
hap_id  <- sub(".*_number\\.txt-", "", g$ID)
blocks  <- unique(hap_all)

## 样本携带的单倍型：1/1
geno_mat <- matrix(NA_character_, nrow=ncol(gt), ncol=0)
samp_ids <- colnames(gt)
hap_list <- list()
for (bl in blocks) {
  idx <- which(hap_all == bl)
  h <- rep(NA_character_, length(samp_ids))
  for (k in idx) h[gt[k, ] == "1/1"] <- hap_id[k]
  hap_list[[bl]] <- data.frame(ID=samp_ids, hap=h)
  hap_list[[bl]]$block <- bl
}
names(hap_list) <- blocks

n_hap_of <- sapply(blocks, function(bl) sum(hap_all==bl))
for (bl in blocks) {
  lv <- paste0("Hap", 1:n_hap_of[bl])
  hap_list[[bl]]$hap <- factor(hap_list[[bl]]$hap, levels=lv)
}

## ---- 表型 ----
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe$V1 <- as.character(phe$V1)
names(phe) <- c("ID", tr7)

## ---- 分组 ----
grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID","group"), comment.char="#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace","Introduced","Pre-1960","1961-1980","1981-2000","2001-2020")
grp$group <- factor(grp$group, levels=grp_order)

## ---- 全性状 ANOVA 筛查 ----
res <- NULL
for (bl in blocks) {
  d0 <- merge(hap_list[[bl]], phe, by="ID")
  for (t in tr7) {
    dd <- d0[!is.na(d0[[t]]) & !is.na(d0$hap), ]
    if (nlevels(droplevels(dd$hap)) < 2) next
    f <- aov(dd[[t]] ~ hap, data=dd)
    s <- summary(f)[[1]]
    p <- s[1,5]
    res <- rbind(res, data.frame(block=bl, trait=t, n=nrow(dd),
                                 F=s[1,4], P=p))
  }
}
print(res %>% mutate(P=signif(P,3)) %>% as.data.frame(), right=FALSE)
write.table(res %>% mutate(across(where(is.numeric), ~signif(.x,5))),
            "ABI_SNRK_hap_anova_all_traits.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: ABI_SNRK_hap_anova_all_traits.txt\n")

## ================= 绘图与保存 =================
pal <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69")
sig_traits <- c("P1","P6","P7")

for (bl in blocks) {
  blname <- sub("^X", "", bl)
  d0 <- merge(hap_list[[bl]], phe, by="ID")
  nh <- nlevels(d0$hap)

  ## ---- 1) 显著性状的 ANOVA + 小提琴图 ----
  ms_all <- NULL
  for (t in sig_traits) {
    dd <- d0[!is.na(d0[[t]]) & !is.na(d0$hap), ]
    fit <- aov(dd[[t]] ~ hap, data=dd)
    tuk <- agricolae::HSD.test(fit, "hap", console=FALSE)
    ltr <- tolower(tuk$groups$groups); names(ltr) <- rownames(tuk$groups)
    ms <- dd %>% group_by(hap) %>%
      summarise(n=n(), mean=mean(.data[[t]]), sd=sd(.data[[t]]),
                median=median(.data[[t]]), sem=sd(.data[[t]])/sqrt(n()),
                .groups="drop") %>% as.data.frame()
    ms$letter <- ltr[as.character(ms$hap)]
    ms$mean_sd <- sprintf("%.3f±%.3f", ms$mean, ms$sd)
    ms$trait <- t; ms$anova_P <- signif(summary(fit)[[1]][1,5], 3)
    ms_all <- rbind(ms_all, ms)

    p <- ggplot(dd, aes(x=hap, y=.data[[t]], fill=hap)) +
      geom_violin(trim=FALSE, scale="width", linewidth=0.3, alpha=0.8,
                  colour="grey25") +
      geom_boxplot(width=0.12, outlier.shape=NA, fill="white",
                   alpha=0.9, linewidth=0.3) +
      geom_text(data=ms, aes(x=hap, y=max(dd[[t]])*1.06, label=letter),
                size=4.6, fontface="bold", inherit.aes=FALSE) +
      scale_fill_manual(values=pal[1:nh]) +
      scale_x_discrete(drop=FALSE) +
      labs(x=paste0("Haplotype of ", blname), y=paste(t, "phenotype"),
           title=paste0(t, " phenotype by haplotype of ", blname,
                        " (ANOVA + Tukey HSD)")) +
      annotate("text", x=0.6, y=max(dd[[t]])*1.10, hjust=0, size=3.2,
               label=paste0("ANOVA P = ", signif(summary(fit)[[1]][1,5], 3))) +
      coord_cartesian(ylim=c(min(dd[[t]]), max(dd[[t]])*1.15)) +
      theme_bw() + theme(legend.position="none", panel.grid.minor=element_blank())
    ggsave(paste0("ABI_SNRK_", blname, "_", t, "_violin.pdf"), p,
           width=7.5, height=6)
  }
  write.table(ms_all, paste0("ABI_SNRK_", blname, "_hap_pheno_means.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  ## ---- 2) 分组单倍型占比（与性状无关，每 block 一张） ----
  d2 <- merge(hap_list[[bl]], grp, by="ID") %>%
    filter(!is.na(group), !is.na(hap))
  prop2 <- d2 %>% group_by(group, hap) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(group) %>% mutate(total=sum(n), prop=n/total) %>% as.data.frame()
  prop2$prop_pct <- sprintf("%.1f%%", prop2$prop*100)
  write.table(prop2[, c("group","hap","n","total","prop","prop_pct")],
              paste0("ABI_SNRK_", blname, "_group_prop.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  w2 <- prop2 %>% select(group, hap, prop_pct) %>%
    pivot_wider(names_from=hap, values_from=prop_pct, values_fill="0.0%") %>%
    as.data.frame()
  write.table(w2, paste0("ABI_SNRK_", blname, "_group_prop_wide.txt"),
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
  ggsave(paste0("ABI_SNRK_", blname, "_group_stack.pdf"), p2, width=8, height=6)
}
cat("DONE all plots/tables\n")
