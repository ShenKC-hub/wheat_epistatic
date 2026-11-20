# qinben_progeny_pheno_8term_model.R
# 在 m123_i3 模型(3个上位网络指标 + 3个两两交互 + 1个三元交互)基础上,
# 加入加性效应指标 n_favorable (qinben_progeny_favorable_hap/{POP}/{POP}_progeny_favorable_P{t}_txt)
# 共 8 项: n_favorable, edges, modules, completeness,
#          edges×modules, edges×completeness, modules×completeness, edges×modules×completeness
# 用精确 LMG(全部 2^8-1=255 个子集 R², 加权公式)分解各项对表型的方差解释率
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
pops$p0 <- as.character(pops$p0); pops$p2 <- as.character(pops$p2)
tr7 <- paste0("P", 1:7)
EMETRICS <- c("n_edges_carried","n_modules_covered","mean_completeness_effective")
TERMS8 <- c("n_favorable", EMETRICS,
            "int_em_mo","int_em_cpe","int_mo_cpe","int_all")

sp <- read.delim("shuangqin_phe2.txt", header=TRUE, check.names=FALSE)
names(sp)[1:2] <- c("group","ID")
norm_id <- function(x) gsub("[_.]", "-", sub("^X(?=[0-9])", "", x, perl=TRUE))

## ---- 精确 LMG: 8 项, 枚举全部 255 个非空子集的 R², 加权公式
## LMG_j = sum_{S 不含 j} |S|!(p-|S|-1)!/p! * (R2(S∪{j}) - R2(S)),  R2(空集)=0
lmg8 <- function(dz, y) {
  p <- length(TERMS8); TSS <- sum((y - mean(y))^2)
  bits <- 2^(0:(p-1))
  R2 <- numeric(2^p); R2[1] <- 0   # mask+1 索引; mask=0 -> 空
  for (m in 1:(2^p-1)) {
    cols <- TERMS8[which(bitwAnd(m, bits) > 0)]
    X <- cbind(1, as.matrix(dz[, cols, drop=FALSE]))
    f <- lm.fit(X, y)
    R2[m+1] <- 1 - sum(f$residuals^2)/TSS
  }
  fact <- factorial(0:p)   # fact[i] = (i-1)!
  lmg <- numeric(p); names(lmg) <- TERMS8
  for (j in seq_len(p)) {
    bitj <- bits[j]; tot <- 0
    for (m in 0:(2^p-1)) {
      if (bitwAnd(m, bitj) != 0) next
      k <- sum(bitwAnd(m, bits) > 0)
      w <- fact[k+1] * fact[p-k] / fact[p+1]
      tot <- tot + w * (R2[m + bitj + 1] - R2[m+1])
    }
    lmg[j] <- tot
  }
  list(R2full = R2[2^p], lmg = lmg)
}

term_rows <- NULL; model_rows <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]
  pdir1 <- file.path("qinben_progeny_module", pop)
  pdir2 <- file.path("qinben_progeny_favorable_hap", pop)
  for (i in 1:7) {
    f1 <- file.path(pdir1, paste0(pop, "_progeny_stats_P", i, "_txt"))
    f2 <- file.path(pdir2, paste0(pop, "_progeny_favorable_P", i, "_txt"))
    if (!file.exists(f1) || !file.exists(f2)) next
    s1 <- read.delim(f1); s2 <- read.delim(f2)
    dm <- merge(s1[, c("progeny", EMETRICS)],
                s2[, c("progeny","n_favorable")], by="progeny")
    dm$ID <- norm_id(dm$progeny)
    d <- merge(dm, sp[sp$group==pop, c("ID", tr7[i])], by="ID")
    d <- d[, c("n_favorable", EMETRICS, tr7[i])]
    d <- d[complete.cases(d), ]
    if (nrow(d) < 30) next
    dz <- as.data.frame(scale(d[, c("n_favorable", EMETRICS)]))
    names(dz) <- c("n_favorable", EMETRICS)
    dz$y <- d[[tr7[i]]]
    dz$int_em_mo  <- dz$n_edges_carried * dz$n_modules_covered
    dz$int_em_cpe <- dz$n_edges_carried * dz$mean_completeness_effective
    dz$int_mo_cpe <- dz$n_modules_covered * dz$mean_completeness_effective
    dz$int_all    <- dz$n_edges_carried * dz$n_modules_covered * dz$mean_completeness_effective

    ## 完整 8 项模型
    f8 <- lm(as.formula(paste("y ~", paste(TERMS8, collapse="+"))), data=dz)
    s8 <- summary(f8); R2_8 <- s8$r.squared
    p8 <- pf(s8$fstatistic[1], s8$fstatistic[2], s8$fstatistic[3], lower.tail=FALSE)
    adj8 <- 1 - (1-R2_8)*(nrow(dz)-1)/(nrow(dz)-length(TERMS8)-1)
    cf8 <- s8$coefficients
    ## 7 项对照模型 (m123_i3, 不含 n_favorable)
    t7 <- setdiff(TERMS8, "n_favorable")
    f7 <- lm(as.formula(paste("y ~", paste(t7, collapse="+"))), data=dz)
    R2_7 <- summary(f7)$r.squared
    an <- anova(f7, f8)
    delta_p <- an$"Pr(>F)"[2]

    ## LMG 分解
    L <- lmg8(dz[, TERMS8], dz$y)
    for (tt in TERMS8) {
      term_rows <- rbind(term_rows, data.frame(
        pop=pop, trait=tr7[i], n=nrow(dz),
        model_R2=round(R2_8,4), model_adjR2=round(adj8,4), model_p=signif(p8,3),
        term=tt,
        beta=ifelse(tt %in% rownames(cf8), round(cf8[tt,"Estimate"],4), NA),
        beta_p=ifelse(tt %in% rownames(cf8), signif(cf8[tt,"Pr(>|t|)"],3), NA),
        lmg=round(L$lmg[tt],4),
        lmg_pct=round(L$lmg[tt]/max(R2_8,1e-9)*100,2)))
    }
    model_rows <- rbind(model_rows, data.frame(
      pop=pop, trait=tr7[i], n=nrow(dz),
      R2_7term=round(R2_7,4), R2_8term=round(R2_8,4),
      delta_R2=round(R2_8-R2_7,4), delta_p=signif(delta_p,3),
      model_p=signif(p8,3)))
  }
  cat(sprintf("%s done\n", pop))
}
term_rows <- term_rows[order(term_rows$pop, term_rows$trait, term_rows$term), ]
write.table(term_rows, "qinben_progeny_pheno_8term_terms_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(model_rows, "qinben_progeny_pheno_8term_models_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: qinben_progeny_pheno_8term_terms_txt / _models_txt\n")

## ---- 按性状汇总 ----
summ <- model_rows %>% group_by(trait) %>%
  summarise(n_models=n(), n_sig_model=sum(model_p<0.05),
            n_sig_delta=sum(delta_p<0.05),
            median_R2_8=round(median(R2_8term),4),
            mean_R2_8=round(mean(R2_8term),4),
            mean_R2_7=round(mean(R2_7term),4),
            mean_delta=round(mean(delta_R2),4),
            .groups="drop") %>% as.data.frame()
lmg_mean <- term_rows %>% group_by(trait, term) %>%
  summarise(mean_lmg=round(mean(lmg),4), .groups="drop") %>%
  pivot_wider(names_from=term, values_from=mean_lmg,
              names_prefix="lmg_") %>% as.data.frame()
summ <- merge(summ, lmg_mean, by="trait")
write.table(summ, "qinben_progeny_pheno_8term_summary_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
print(summ, right=FALSE)

## ---- 图: LMG 堆叠条形图 + R² 箱线图 ----
lmg_bar <- term_rows %>% group_by(trait, term) %>%
  summarise(mean_lmg=mean(lmg), .groups="drop") %>% as.data.frame()
p1 <- ggplot(lmg_bar, aes(x=trait, y=mean_lmg, fill=term)) +
  geom_col(position="stack", colour="black", linewidth=0.2, width=0.7) +
  labs(title="LMG decomposition of 8-term model (3 epistasis metrics + interactions + n_favorable)",
       x="Trait", y="Mean LMG (phenotype variance explained)", fill="Term") +
  theme_bw() + theme(axis.text.x=element_text(angle=0))
ggsave("qinben_progeny_pheno_8term_LMG_bar.pdf", p1, width=10, height=6)

p2 <- ggplot(model_rows, aes(x=trait, y=R2_8term)) +
  geom_boxplot(fill="#8DD3C7", outlier.size=0.9, width=0.6) +
  geom_jitter(width=0.12, size=0.5, alpha=0.4, colour="grey25") +
  labs(title="8-term model R2 (with n_favorable added)",
       subtitle="7-term baseline (m123_i3) shown in red",
       x="Trait", y="Model R2") +
  geom_hline(data=. %>% group_by(trait) %>% summarise(m=mean(R2_7term)),
             aes(yintercept=m), colour="#E41A1C", linetype=2) +
  theme_bw()
ggsave("qinben_progeny_pheno_8term_R2_boxplot.pdf", p2, width=8, height=6)
cat("saved: 2 pdfs\nDONE\n")