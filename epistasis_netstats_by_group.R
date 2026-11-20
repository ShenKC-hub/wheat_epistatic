# 育种历程分组统计：各性状(P1-P7)有利上位网络的 node_count / edge_count
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID", "group"), comment.char="#")
grp$ID <- as.character(grp$ID)

## 分组按育种历程排序
grp_order <- c("Landrace", "Introduced", "Pre-1960",
               "1961-1980", "1981-2000", "2001-2020")

## 读取 7 个性状的网络统计
dat <- list()
for (i in 1:7) {
  st <- read.delim(paste0("gwas_epistasis_netstats_P", i, "_txt"), row.names=1)
  stopifnot(colnames(st)[1] == "node_count", colnames(st)[2] == "edge_count")
  d <- data.frame(ID=rownames(st), trait=paste0("P", i),
                  node_count=st$node_count, edge_count=st$edge_count)
  dat[[i]] <- d
}
dat <- do.call(rbind, dat)
dat <- merge(dat, grp, by="ID")
dat$trait  <- factor(dat$trait, levels=paste0("P", 1:7))
dat$group  <- factor(dat$group, levels=grp_order)

## ---- 汇总表：均值 ± 标准差 ----
summ <- dat %>%
  group_by(trait, group) %>%
  summarise(n=n(),
            node_mean=mean(node_count), node_sd=sd(node_count),
            edge_mean=mean(edge_count), edge_sd=sd(edge_count),
            .groups="drop") %>%
  as.data.frame()
summ$node_ms <- ifelse(is.na(summ$node_sd), summ$node_mean,
                       sprintf("%.2f±%.2f", summ$node_mean, summ$node_sd))
summ$edge_ms <- ifelse(is.na(summ$edge_sd), summ$edge_mean,
                       sprintf("%.2f±%.2f", summ$edge_mean, summ$edge_sd))

out_wide <- summ %>%
  transmute(trait, group, n,
            node_mean, node_sd, node_mean_sd=node_ms,
            edge_mean, edge_sd, edge_mean_sd=edge_ms) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
write.table(out_wide, "epistasis_netstats_by_group.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: epistasis_netstats_by_group.txt\n")
print(out_wide, right=FALSE)

## ---- 柱状图：均值 ± 标准差 ----
pd <- position_dodge(0.8)
p1 <- ggplot(summ, aes(x=group, y=node_mean, fill=trait)) +
  geom_bar(stat="identity", position=pd, width=0.8, colour="black", linewidth=0.2) +
  geom_errorbar(aes(ymin=node_mean-node_sd, ymax=node_mean+node_sd),
                position=pd, width=0.25, linewidth=0.3) +
  labs(title="Favorable epistasis network: node_count by breeding era",
       x="Breeding era", y="Mean node_count (±SD)", fill="Trait") +
  theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1))

p2 <- ggplot(summ, aes(x=group, y=edge_mean, fill=trait)) +
  geom_bar(stat="identity", position=pd, width=0.8, colour="black", linewidth=0.2) +
  geom_errorbar(aes(ymin=edge_mean-edge_sd, ymax=edge_mean+edge_sd),
                position=pd, width=0.25, linewidth=0.3) +
  labs(title="Favorable epistasis network: edge_count by breeding era",
       x="Breeding era", y="Mean edge_count (±SD)", fill="Trait") +
  theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1))

ggsave("epistasis_netstats_node_by_group.pdf", p1, width=10, height=6)
ggsave("epistasis_netstats_edge_by_group.pdf", p2, width=10, height=6)
cat("saved: epistasis_netstats_node_by_group.pdf / epistasis_netstats_edge_by_group.pdf\n")
