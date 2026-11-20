# qinben_progeny_pheno_correlation.R
# 任务1: 每群体每性状, 子代上位网络指标(n_edges_carried, n_modules_covered,
#        mean_completeness_effective) 与表型的相关性 (P1-P7 分组)
# 任务2: 亲本互补程度(3指标) 与 上位效应遗传力 h2_epi / 子代表型增加潜力的
#        相关性 (每性状 6 组), 点图 + 相关统计
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
pops$p0 <- as.character(pops$p0); pops$p2 <- as.character(pops$p2)
tr7 <- paste0("P", 1:7)

## ---- 表型: 双亲群体子代 ----
sp <- read.delim("shuangqin_phe2.txt", header=TRUE, check.names=FALSE)
names(sp)[1:2] <- c("group","ID")

## ---- ID 规范化: 去掉数字前的X, _或. 转 - ----
norm_id <- function(x) gsub("[_.]", "-", sub("^X(?=[0-9])", "", x, perl=TRUE))

## ---- 亲本表型查询: 先 shuangqin_phe2, 后 phe2.txt ----
phe_nat <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe_nat$V1 <- as.character(phe_nat$V1)
names(phe_nat) <- c("ID", tr7)
get_parent <- function(pid, t) {
  v1 <- sp[[t]][sp$ID == pid]
  if (length(v1) == 1 && !is.na(v1)) return(v1)
  v2 <- phe_nat[[t]][phe_nat$ID == pid]
  if (length(v2) == 1 && !is.na(v2)) return(v2)
  NA_real_
}

## ================= 任务1 =================
cor_long <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]
  for (i in 1:7) {
    st <- tryCatch(read.delim(file.path("qinben_progeny_module", pop,
                                        paste0(pop, "_progeny_stats_P", i, "_txt"))),
                   error=function(e) NULL)
    if (is.null(st)) next
    st$ID <- norm_id(st$progeny)
    d <- merge(st, sp[sp$group==pop, c("ID", tr7[i])], by.x="ID", by.y="ID")
    if (nrow(d) < 10) next
    for (m in c("n_edges_carried","n_modules_covered","mean_completeness_effective")) {
      ct <- cor.test(d[[m]], d[[tr7[i]]], method="spearman", exact=FALSE)
      cp <- cor.test(d[[m]], d[[tr7[i]]], method="pearson")
      cor_long <- rbind(cor_long, data.frame(
        pop=pop, trait=tr7[i], metric=m, n=nrow(d),
        rho_spear=round(ct$estimate,3), p_spear=signif(ct$p.value,3),
        r_pearson=round(cp$estimate,3), p_pearson=signif(cp$p.value,3)))
    }
  }
}
cor_long <- cor_long %>% mutate(sig_spear=ifelse(p_spear<0.001,"***",
                        ifelse(p_spear<0.01,"**",ifelse(p_spear<0.05,"*","ns"))))
write.table(cor_long, "qinben_progeny_pheno_correlation_txt",
            sep="\t", quote=FALSE, row.names=FALSE)

summ1 <- cor_long %>% group_by(trait, metric) %>%
  summarise(n_pops=n(), n_sig=sum(p_spear<0.05),
            median_rho=round(median(rho_spear),3),
            mean_rho=round(mean(rho_spear),3),
            min_rho=round(min(rho_spear),3), max_rho=round(max(rho_spear),3),
            .groups="drop") %>% as.data.frame()
write.table(summ1, "qinben_progeny_pheno_correlation_summary_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
print(summ1, right=FALSE)

p1 <- ggplot(cor_long, aes(x=trait, y=rho_spear, fill=metric)) +
  geom_hline(yintercept=0, colour="grey40") +
  geom_boxplot(position=position_dodge(0.75), width=0.7, outlier.size=0.8) +
  scale_fill_manual(values=c("#8DD3C7","#FB8072","#80B1D3")) +
  labs(title="Correlation of progeny epistasis-network metrics with phenotype",
       subtitle="Each point = one population (Spearman rho within population)",
       x="Trait", y="Spearman rho", fill="Metric") +
  theme_bw()
ggsave("qinben_progeny_pheno_correlation_boxplot.pdf", p1, width=9, height=6)

## ================= 任务2 =================
h2 <- read.csv("biparental_23pops_h2_epi_npairs_csv")

comp_all <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]
  for (i in 1:7) {
    ov <- tryCatch(read.delim(file.path("qinben_progeny_module", pop,
                      paste0(pop, "_overlap_modules_P", i, "_txt"))),
                   error=function(e) NULL)
    if (is.null(ov)) next
    n_both <- sum(ov$edges_p0>0 & ov$edges_p2>0)
    n_p0only <- sum(ov$edges_p0>0 & ov$edges_p2==0)
    n_p2only <- sum(ov$edges_p0==0 & ov$edges_p2>0)
    present <- sum(ov$edges_p0>0 | ov$edges_p2>0)
    metric_module_comp <- ifelse(present>0, (n_p0only+n_p2only)/present, NA)
    metric_edge_comp <- mean(1-ov$overlap_coeff, na.rm=TRUE)
    metric_comp_union <- mean(ov$comp_union)
    ## 遗传力与子代增益
    hh <- h2$h2_epi[h2$pop==pop & h2$trait==tr7[i]]
    hh <- ifelse(length(hh)==1, hh, NA)
    spv <- sp[[tr7[i]]][sp$group==pop]
    prog_mean <- mean(spv, na.rm=TRUE)
    a <- get_parent(pops$p0[r], tr7[i]); b <- get_parent(pops$p2[r], tr7[i])
    gain <- ifelse(!is.na(a) & !is.na(b), prog_mean - (a+b)/2, NA)
    comp_all <- rbind(comp_all, data.frame(
      pop=pop, trait=tr7[i],
      n_mod_present=present, n_mod_both=n_both,
      n_mod_p0only=n_p0only, n_mod_p2only=n_p2only,
      metric_module_comp=round(metric_module_comp,4),
      metric_edge_comp=round(metric_edge_comp,4),
      metric_comp_union=round(metric_comp_union,4),
      h2_epi=hh, prog_mean=round(prog_mean,3),
      parent_p0=round(a,3), parent_p2=round(b,3), gain=round(gain,3)))
  }
}
write.table(comp_all, "parent_complementarity_h2_gain_txt",
            sep="\t", quote=FALSE, row.names=FALSE)

## 相关性: 每性状 3指标 × 2目标 = 6 组
cors2 <- NULL
for (t in tr7) {
  d <- comp_all %>% filter(trait==t) %>%
    select(pop, metric_module_comp, metric_edge_comp, metric_comp_union,
           h2_epi, gain) %>% as.data.frame()
  for (m in c("metric_module_comp","metric_edge_comp","metric_comp_union")) {
    for (y in c("h2_epi","gain")) {
      dd <- d[!is.na(d[[m]]) & !is.na(d[[y]]), ]
      if (nrow(dd) < 5) {
        cors2 <- rbind(cors2, data.frame(trait=t, x=m, y=y, n=nrow(dd),
          rho=NA, p=NA, r=NA, p_pearson=NA)); next
      }
      ct <- cor.test(dd[[m]], dd[[y]], method="spearman", exact=FALSE)
      cp <- cor.test(dd[[m]], dd[[y]], method="pearson")
      cors2 <- rbind(cors2, data.frame(trait=t, x=m, y=y, n=nrow(dd),
        rho=round(ct$estimate,3), p=signif(ct$p.value,3),
        r=round(cp$estimate,3), p_pearson=signif(cp$p.value,3)))
    }
  }
}
write.table(cors2, "parent_complementarity_correlations_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
print(cors2, right=FALSE)

## 点图: 每性状一页, 2行(h2/gain) × 3列(指标)
pdf("parent_complementarity_scatter_7traits.pdf", width=12, height=7)
for (t in tr7) {
  d <- comp_all %>% filter(trait==t) %>%
    select(pop, metric_module_comp, metric_edge_comp, metric_comp_union,
           h2_epi, gain) %>%
    pivot_longer(-c(pop, h2_epi, gain),
                 names_to="metric", values_to="value") %>% as.data.frame()
  p <- ggplot(d, aes(x=value, y=h2_epi)) +
    geom_point(colour="#E41A1C", size=2.4) +
    geom_smooth(method="lm", se=TRUE, colour="grey30", linewidth=0.5) +
    facet_wrap(~metric, scales="free_x", nrow=1) +
    labs(title=paste0(t, ": parent complementarity vs h2_epi (red) / progeny gain"),
         x="Complementarity metric", y="h2_epi / progeny gain") +
    theme_bw()
  d2 <- d; d2$h2_epi <- d2$gain
  p2 <- p + geom_point(data=d2, aes(x=value, y=h2_epi), colour="#377EB8",
                       size=2.4) +
    geom_smooth(data=d2, aes(x=value, y=h2_epi), method="lm", se=TRUE,
                colour="grey30", linewidth=0.5, linetype=2)
  print(p2)
}
dev.off()
cat("saved: parent_complementarity_scatter_7traits.pdf\nDONE\n")