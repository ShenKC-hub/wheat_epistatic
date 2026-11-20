# 三部分统计：分性状绘图 + 农艺性状分组 + 不利上位网络
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(ggplot2)
library(dplyr)
library(tidyr)

grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID", "group"), comment.char="#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace", "Introduced", "Pre-1960",
               "1961-1980", "1981-2000", "2001-2020")
grp$group <- factor(grp$group, levels=grp_order)
tr7 <- paste0("P", 1:7)

## 通用：均值±SD汇总
summarise_ms <- function(d, val) {
  nm <- val
  d %>% group_by(trait, group) %>%
    summarise(n=n(),
              mean_val=mean(.data[[val]], na.rm=TRUE),
              sd_val=sd(.data[[val]], na.rm=TRUE), .groups="drop") %>%
    rename(!!nm := mean_val) %>%
    rename(!!paste0(nm, "_sd") := sd_val) %>%
    as.data.frame()
}

## 通用：分面柱状图（每性状一个面，y 轴自由）
plot_facet <- function(summ, mcol, scol, ylab, file) {
  summ$.sd <- summ[[scol]]
  p <- ggplot(summ, aes(x=group, y=.data[[mcol]], fill=trait)) +
    geom_bar(stat="identity", colour="black", linewidth=0.2) +
    geom_errorbar(aes(ymin=.data[[mcol]]-.sd, ymax=.data[[mcol]]+.sd),
                  width=0.25, linewidth=0.3) +
    facet_wrap(~trait, scales="free_y") +
    labs(x="Breeding era", y=ylab, fill="Trait") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          strip.text=element_text(face="bold"),
          legend.position="none")
  ggsave(file, p, width=12, height=9)
}

## ================= 1. 有利上位效应：分面绘图 =================
dat <- list()
for (i in 1:7) {
  st <- read.delim(paste0("gwas_epistasis_netstats_P", i, "_txt"), row.names=1)
  dat[[i]] <- data.frame(ID=rownames(st), trait=tr7[i],
                         node_count=st$node_count, edge_count=st$edge_count)
}
fav <- do.call(rbind, dat) %>% merge(grp, by="ID") %>%
  mutate(trait=factor(trait, levels=tr7))

fav_n <- summarise_ms(fav, "node_count"); names(fav_n)[4:5] <- c("node_mean","node_sd")
fav_e <- summarise_ms(fav, "edge_count"); names(fav_e)[4:5] <- c("edge_mean","edge_sd")

plot_facet(fav_n, "node_mean", "node_sd",
           "Mean node_count (±SD)", "epistasis_netstats_node_by_group_facet.pdf")
plot_facet(fav_e, "edge_mean", "edge_sd",
           "Mean edge_count (±SD)", "epistasis_netstats_edge_by_group_facet.pdf")
cat("1 done: favorable facet plots\n")

## ================= 2. 7个农艺性状(phe2.txt)按育种历程统计 =================
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe$V1 <- as.character(phe$V1)
names(phe) <- c("ID", tr7)

phe_l <- phe %>%
  pivot_longer(all_of(tr7), names_to="trait", values_to="value") %>%
  merge(grp, by="ID") %>%
  mutate(trait=factor(trait, levels=tr7))

phe_s <- phe_l %>% group_by(trait, group) %>%
  summarise(n=n(), mean_val=mean(value, na.rm=TRUE),
            sd_val=sd(value, na.rm=TRUE), .groups="drop") %>%
  as.data.frame() %>%
  mutate(mean_sd=sprintf("%.3f±%.3f", mean_val, sd_val),
         across(where(is.numeric), ~round(.x, 3)))
write.table(phe_s, "phenotype_by_group.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("2 done: phenotype_by_group.txt\n")

plot_facet(phe_s, "mean_val", "sd_val",
           "Mean trait value (±SD)", "phenotype_by_group_facet.pdf")
cat("  facet plot saved\n")

## ================= 3. 不利上位效应(inferior) =================
inf <- list()
for (i in 1:7) {
  st <- read.delim(paste0("inferior/gwas_epistasis_netstats_inferior_P", i, "_txt"),
                   row.names=1)
  inf[[i]] <- data.frame(ID=rownames(st), trait=tr7[i],
                         node_count=st$node_count, edge_count=st$edge_count)
}
inf <- do.call(rbind, inf) %>% merge(grp, by="ID") %>%
  mutate(trait=factor(trait, levels=tr7))

inf_n <- summarise_ms(inf, "node_count"); names(inf_n)[4:5] <- c("node_mean","node_sd")
inf_e <- summarise_ms(inf, "edge_count"); names(inf_e)[4:5] <- c("edge_mean","edge_sd")

inf_n$node_mean_sd <- sprintf("%.2f±%.2f", inf_n$node_mean, inf_n$node_sd)
inf_e$edge_mean_sd <- sprintf("%.2f±%.2f", inf_e$edge_mean, inf_e$edge_sd)
out_inf <- merge(inf_n, inf_e, by=c("trait","group","n")) %>%
  mutate(across(all_of(c("node_mean","node_sd","edge_mean","edge_sd")), ~round(.x, 3))) %>%
  select(trait, group, n, node_mean, node_sd, node_mean_sd,
         edge_mean, edge_sd, edge_mean_sd)
write.table(out_inf, "epistasis_netstats_inferior_by_group.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("3 done: epistasis_netstats_inferior_by_group.txt\n")

plot_facet(inf_n, "node_mean", "node_sd",
           "Mean node_count (±SD)", "epistasis_netstats_inferior_node_by_group_facet.pdf")
plot_facet(inf_e, "edge_mean", "edge_sd",
           "Mean edge_count (±SD)", "epistasis_netstats_inferior_edge_by_group_facet.pdf")
cat("  inferior facet plots saved\n")

## 概览
cat("\n--- inferior summary (first rows) ---\n")
print(head(out_inf, 6), right=FALSE)
