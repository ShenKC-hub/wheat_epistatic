# qinben_progeny_pheno_multivariate.R
# 三指标多元回归 + LMG 表型解释率分解 (22群体 × P1-P7, 群体内分析)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
pops$p0 <- as.character(pops$p0); pops$p2 <- as.character(pops$p2)
tr7 <- paste0("P", 1:7)
METRICS <- c("n_edges_carried","n_modules_covered","mean_completeness_effective")

sp <- read.delim("shuangqin_phe2.txt", header=TRUE, check.names=FALSE)
names(sp)[1:2] <- c("group","ID")
norm_id <- function(x) gsub("[_.]", "-", sub("^X(?=[0-9])", "", x, perl=TRUE))

## LMG: 3指标所有进入顺序的偏R²平均 (手工实现)
lmg3 <- function(d, t) {
  R2full <- summary(lm(d[[t]] ~ ., data=d[, METRICS]))$r.squared
  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  inc <- integer(3)
  for (k in seq_len(6)) {
    R2prev <- 0
    ord <- perms[[k]]
    for (j in seq_len(3)) {
      cols <- METRICS[ord[1:j]]
      R2 <- summary(lm(d[[t]] ~ ., data=d[, cols, drop=FALSE]))$r.squared
      inc[ord[j]] <- inc[ord[j]] + (R2 - R2prev)
      R2prev <- R2
    }
  }
  inc <- inc / 6
  names(inc) <- METRICS
  list(R2full=R2full, lmg=inc)
}

res <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]
  for (i in 1:7) {
    f <- file.path("qinben_progeny_module", pop,
                   paste0(pop, "_progeny_stats_P", i, "_txt"))
    st <- tryCatch(read.delim(f), error=function(e) NULL)
    if (is.null(st)) next
    st$ID <- norm_id(st$progeny)
    d <- merge(st, sp[sp$group==pop, c("ID", tr7[i])], by="ID")
    d <- d[, c(METRICS, tr7[i])]
    d <- d[complete.cases(d), ]
    if (nrow(d) < 30) next
    ## 标准化
    dz <- as.data.frame(scale(d))
    names(dz) <- c(METRICS, "y")
    fit <- lm(y ~ ., data=dz)
    s <- summary(fit)
    R2 <- s$r.squared
    model_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                  lower.tail=FALSE)
    cf <- s$coefficients
    ## VIF
    vif <- sapply(METRICS, function(m) {
      fm <- lm(dz[[m]] ~ ., data=dz[, setdiff(METRICS, m), drop=FALSE])
      1/(1-summary(fm)$r.squared) })
    ## LMG 分解
    L <- suppressWarnings(lmg3(dz, "y"))
    ## 偏 R² (控制其余两个指标)
    pr2 <- sapply(METRICS, function(m) {
      f1 <- lm(y ~ ., data=dz[, c(setdiff(METRICS, m), "y")])
      f2 <- lm(y ~ ., data=dz[, c(METRICS, "y")])
      an <- anova(f1, f2)
      an$"Sum of Sq"[2] / (an$"Sum of Sq"[2] + an$"RSS"[2]) })
    for (m in METRICS) {
      res <- rbind(res, data.frame(
        pop=pop, trait=tr7[i], n=nrow(dz),
        model_R2=round(R2,4), model_p=signif(model_p,3),
        metric=m,
        beta=round(cf[m,"Estimate"],4), beta_p=signif(cf[m,"Pr(>|t|)"],3),
        vif=round(vif[m],2),
        lmg=round(L$lmg[m],4),
        lmg_pct=round(L$lmg[m]/max(R2,1e-9)*100,2),
        partial_R2=round(pr2[m],4)))
    }
  }
}
write.table(res, "qinben_progeny_pheno_multivariate_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: qinben_progeny_pheno_multivariate_txt\n")

## ---- 按性状汇总 ----
summ <- res %>% group_by(trait) %>%
  summarise(n_models=length(unique(pop)),
            n_sig_model=sum(model_p<0.05, na.rm=TRUE),
            median_R2=round(median(model_R2),4),
            mean_R2=round(mean(model_R2),4),
            lmg_edges=round(mean(lmg[metric=="n_edges_carried"]),4),
            lmg_modules=round(mean(lmg[metric=="n_modules_covered"]),4),
            lmg_comp=round(mean(lmg[metric=="mean_completeness_effective"]),4),
            median_vif=round(median(vif),2),
            .groups="drop") %>% as.data.frame()
write.table(summ, "qinben_progeny_pheno_multivariate_summary_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
print(summ, right=FALSE)

## ---- 图 ----
lv <- paste0("P", 1:7)
res$trait <- factor(res$trait, levels=lv)
p1 <- ggplot(unique(res[, c("pop","trait","model_R2","model_p")]),
             aes(x=trait, y=model_R2)) +
  geom_hline(yintercept=0, colour="grey40") +
  geom_boxplot(fill="#8DD3C7", outlier.size=0.9, width=0.6) +
  geom_jitter(width=0.12, size=0.5, alpha=0.4, colour="grey25") +
  labs(title="Joint explanatory power of 3 epistasis-network metrics (multiple regression)",
       subtitle="Each point = one population; y = model R2 of phenotype ~ 3 metrics",
       x="Trait", y="Model R2") + theme_bw()
ggsave("qinben_progeny_pheno_multivariate_R2_boxplot.pdf", p1, width=8, height=6)

lmg_bar <- summ %>%
  select(trait, lmg_edges, lmg_modules, lmg_comp) %>%
  pivot_longer(-trait, names_to="metric", values_to="lmg_mean")
p2 <- ggplot(lmg_bar, aes(x=trait, y=lmg_mean, fill=metric)) +
  geom_col(position="stack", colour="black", linewidth=0.2, width=0.7) +
  scale_fill_manual(values=c("#8DD3C7","#FB8072","#80B1D3"),
                    labels=c("n_edges_carried","n_modules_covered",
                             "mean_completeness_effective")) +
  labs(title="Average LMG share of phenotype variance explained",
       x="Trait", y="Mean LMG (variance explained)", fill="Metric") +
  theme_bw()
ggsave("qinben_progeny_pheno_multivariate_LMG_bar.pdf", p2, width=8, height=6)
cat("saved: 2 pdfs\nDONE\n")