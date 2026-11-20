# 有利上位效应对数目(edge_count)与表型的关系：分箱箱线图 + 光滑曲线
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

tr7 <- paste0("P", 1:7)

## 表型
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe$V1 <- as.character(phe$V1)
names(phe) <- c("ID", tr7)

## 每材料每性状的有利上位网络 edge_count
dat <- list()
for (i in 1:7) {
  st <- read.delim(paste0("gwas_epistasis_netstats_P", i, "_txt"), row.names=1)
  dat[[i]] <- data.frame(ID=rownames(st), trait=tr7[i], edge=st$edge_count)
}
dat <- do.call(rbind, dat) %>% merge(phe, by="ID") %>%
  pivot_longer(all_of(tr7), names_to="pheno", values_to="value") %>%
  filter(trait == pheno) %>%          # 性状与其对应的有利网络配对
  select(ID, trait, edge, value) %>%
  filter(!is.na(value)) %>%
  mutate(trait=factor(trait, levels=tr7))

cat("pairs per trait:\n"); print(table(dat$trait))

## 按每性状 edge_count 分位数分 5 组（各组成员数相近），记录区间范围
dat <- dat %>% group_by(trait) %>%
  mutate(bin=cut(edge,
                 breaks=unique(quantile(edge, probs=seq(0,1,0.2))),
                 include.lowest=TRUE)) %>%
  ungroup()

## 分组区间标签与组中点
lv <- levels(dat$bin)
lab <- gsub("\\(", "[", lv)
bin_lab <- setNames(lab, lab)
dat <- dat %>%
  mutate(bin=factor(bin, levels=lv, labels=lab),
         lo=as.numeric(sub("^\\[([^,]+),.*", "\\1", lab[as.numeric(factor(bin))])) )
## 组中点（用于数值x轴，便于叠加光滑曲线）
edges <- dat %>% distinct(trait, bin)
mid <- dat %>% group_by(trait, bin) %>%
  summarise(bin_min=min(edge), bin_max=max(edge), .groups="drop") %>%
  mutate(bin_mid=(bin_min+bin_max)/2)
dat <- dat %>% left_join(mid %>% select(trait, bin, bin_mid),
                         by=c("trait","bin"))

## 分组统计表
summ <- dat %>% group_by(trait, bin) %>%
  summarise(n=n(), bin_min=min(edge), bin_max=max(edge),
            mean_val=mean(value), sd_val=sd(value), median_val=median(value),
            .groups="drop") %>%
  as.data.frame() %>%
  mutate(mean_sd=sprintf("%.3f±%.3f", mean_val, sd_val),
         across(where(is.numeric), ~round(.x, 3)))
write.table(summ, "favorable_edge_pheno_group_stats.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: favorable_edge_pheno_group_stats.txt\n")

## 散点图 + 线性回归（仅保留 edge<1000，横坐标 0-1000；每性状一页）
pdf("favorable_edge_vs_pheno_scatter_lm_edge_lt1000.pdf", width=9, height=7)
for (t in tr7) {
  d <- dat %>% filter(trait==t, edge<1000)
  ## 线性回归
  fit <- lm(value ~ edge, data=d)
  s   <- summary(fit)
  r   <- signif(sqrt(s$r.squared) * sign(s$coefficients[2,1]), 3)
  pv  <- signif(s$coefficients[2,4], 3)
  cat(sprintf("%s: n=%d, r=%.3f, P=%.3g\n", t, nrow(d), r, pv))
  p <- ggplot(d, aes(x=edge, y=value)) +
    geom_point(shape=21, fill="#4DAF4A", colour="grey30",
               alpha=0.55, size=1.8, stroke=0.3) +
    geom_smooth(method="lm", se=TRUE, colour="#E41A1C", fill="#E41A1C",
                alpha=0.12, linewidth=0.7) +
    annotate("text", x=Inf, y=Inf, hjust=1.08, vjust=1.6, size=4.2,
             label=sprintf("r = %.3f,  P = %.2e", r, pv),
             fontface="italic") +
    scale_x_continuous(limits=c(0, 1000)) +
    coord_cartesian(xlim=c(0, 1000)) +
    labs(x="Number of favorable epistasis pairs (<1000)",
         y="Phenotype value",
         title=paste0("Favorable epistasis pair count vs phenotype (", t, ")"),
         subtitle="points = materials; red line = linear fit with 95% CI") +
    theme_bw()
  print(p)
}
dev.off()
cat("saved: favorable_edge_vs_pheno_scatter_lm_edge_lt1000.pdf (7 pages)\n")
