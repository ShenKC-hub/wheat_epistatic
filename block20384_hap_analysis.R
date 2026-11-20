# X1A_block_20384 五种单倍型：P7性状ANOVA+小提琴图；分组单倍型占比堆积图
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

block_id <- "X1A_block_20384"

## ---- 读基因型：5行单倍型 × 681样本 ----
g <- read.delim("gwas_all_block_genotype.txt", header=TRUE, comment.char="#",
                check.names=FALSE)
g <- g[grepl(block_id, g$ID), , drop=FALSE]
hap_name <- sub(paste0("^", block_id, "_number\\.txt-"), "", g$ID)  # Hap1..Hap5
gt <- as.matrix(g[, -(1:9)])
storage.mode(gt) <- "integer"
stopifnot(nrow(gt) == 5)

## 每个样本携带的单倍型（值为2）
hap <- rep(NA_character_, ncol(gt))
for (k in 1:5) hap[gt[k, ] == 2] <- hap_name[k]
dat <- data.frame(ID=colnames(gt), hap=factor(hap, levels=paste0("Hap", 1:5)))
cat("样品单倍型分布:\n"); print(table(dat$hap, useNA="ifany"))

## ---- 表型 P7 ----
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe$V1 <- as.character(phe$V1)
phe <- phe[, c(1, 8)]
names(phe) <- c("ID", "P7")

d1 <- merge(dat, phe, by="ID", all=FALSE)
d1 <- d1[!is.na(d1$P7) & !is.na(d1$hap), ]

## ================= 任务1：ANOVA + 小提琴图 =================
fit <- aov(P7 ~ hap, data=d1)
print(summary(fit))
tuk <- agricolae::HSD.test(fit, "hap", console=FALSE)
lev <- rownames(tuk$groups)
ltr <- tolower(tuk$groups$groups)
names(ltr) <- lev
cat("\nTukey HSD 分组:\n"); print(tuk$groups)

ms1 <- d1 %>% group_by(hap) %>%
  summarise(n=n(), mean=mean(P7), sd=sd(P7), median=median(P7),
            sem=sd(P7)/sqrt(n()), .groups="drop") %>%
  as.data.frame()
ms1$letter <- ltr[as.character(ms1$hap)]
ms1$mean_sd <- sprintf("%.3f±%.3f", ms1$mean, ms1$sd)
write.table(ms1[, c("hap","n","mean","sd","median","sem","letter","mean_sd")],
            "block20384_hap_P7_anova_means.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: block20384_hap_P7_anova_means.txt\n")

ylab <- expression("P7 phenotype ("*10^6~"µm"^3*")")  # 若单位不同可修改
p1 <- ggplot(d1, aes(x=hap, y=P7, fill=hap)) +
  geom_violin(trim=FALSE, scale="width", linewidth=0.3, alpha=0.8,
              colour="grey25") +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white",
               alpha=0.9, linewidth=0.3) +
  geom_text(data=ms1, aes(x=hap, y=max(d1$P7)*1.06, label=letter),
            size=5, fontface="bold", inherit.aes=FALSE) +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3")) +
  labs(x="Haplotype of 1A_block_20384", y=ylab,
       title="P7 phenotype by haplotype (ANOVA + Tukey HSD)") +
  annotate("text", x=0.6, y=max(d1$P7)*1.10, hjust=0, size=3.2,
           label=paste0("ANOVA P = ", signif(summary(fit)[[1]][1,5], 3))) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.minor=element_blank())
ymax <- max(d1$P7)*1.15
p1 <- p1 + coord_cartesian(ylim=c(min(d1$P7), ymax))
ggsave("block20384_hap_P7_violin.pdf", p1, width=7, height=6)
cat("saved: block20384_hap_P7_violin.pdf\n")

## ================= 任务2：分组单倍型占比堆积图 =================
grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID","group"), comment.char="#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace","Introduced","Pre-1960","1961-1980","1981-2000","2001-2020")
grp$group <- factor(grp$group, levels=grp_order)

d2 <- merge(dat, grp, by="ID") %>% filter(!is.na(group), !is.na(hap))

prop2 <- d2 %>% group_by(group, hap) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(group) %>%
  mutate(total=sum(n), prop=n/total) %>%
  as.data.frame()
prop2$prop_pct <- sprintf("%.1f%%", prop2$prop*100)

out2 <- prop2 %>%
  select(group, hap, n, total, prop, prop_pct) %>%
  mutate(across(c(prop), ~round(.x, 4)))
write.table(out2, "block20384_hap_by_group_prop.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: block20384_hap_by_group_prop.txt\n")

## 宽表形式（行=分组，列=5种单倍型占比）
w2 <- prop2 %>% select(group, hap, prop_pct) %>%
  pivot_wider(names_from=hap, values_from=prop_pct, values_fill="0.0%") %>%
  as.data.frame()
write.table(w2, "block20384_hap_by_group_prop_wide.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

p2 <- ggplot(prop2, aes(x=group, y=prop, fill=hap)) +
  geom_col(position="stack", colour="black", linewidth=0.2, width=0.75) +
  geom_text(aes(label=ifelse(prop>=0.05, paste0(round(prop*100), "%"), "")),
            position=position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3")) +
  scale_x_discrete(limits=grp_order) +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(x="Breeding era", y="Proportion", fill="Haplotype",
       title="Haplotype composition of 1A_block_20384 by breeding era") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        panel.grid.minor=element_blank())
ggsave("block20384_hap_by_group_stack.pdf", p2, width=8, height=6)
cat("saved: block20384_hap_by_group_stack.pdf\n")

print(w2)