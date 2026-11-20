# parent_complementarity_h2_gain_v2.R
# 任务2更新版:
#   (a) 23GZM-295 的田间编号映射: SZS->23SZS-219, ZJ->24ZJ-249, JZ->23JZ-96,
#       SZM->23SZM-1, SZY->23SZY-1, EZX->24EZX-1, GYGZ->23GYGZ-19, GZH->23GZH-29
#   (b) 子代网络指标(22群体×7性状): 1/2/3指标+交互效应的多元回归,
#       LMG 分解表型解释率
#   (c) 亲本互补程度(3指标) 与 h2_epi / gain 的相关性 (Pearson + Spearman)
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

phe_nat <- read.delim("phe2.txt", header=FALSE, comment.char="#")
phe_nat$V1 <- as.character(phe_nat$V1)
names(phe_nat) <- c("ID", tr7)

## ---- 23GZM-295 田间编号映射 (按群体) ----
GZM_MAP <- c(SZS="23SZS-219", ZJ="24ZJ-249", JZ="23JZ-96", SZM="23SZM-1",
             SZY="23SZY-1", EZX="24EZX-1", GYGZ="23GYGZ-19", GZH="23GZH-29")
## ---- 其他亲本补充映射 (phe2.txt中为NA者) ----
EXTRA_MAP <- c(EZX="24EZX-181")   # 24GYL-181 -> 24EZX-181 (EZX系内田间编号)

get_parent <- function(pid, pop, t) {
  if (pid == "23GZM-295" && pop %in% names(GZM_MAP)) pid <- GZM_MAP[[pop]]
  if (pid %in% names(EXTRA_MAP)) pid <- EXTRA_MAP[[pid]]
  v1 <- sp[[t]][sp$ID == pid]
  if (length(v1) == 1 && !is.na(v1)) return(v1)
  v2 <- phe_nat[[t]][phe_nat$ID == pid]
  if (length(v2) == 1 && !is.na(v2)) return(v2)
  NA_real_
}

## ================= (b) 子代指标多元回归: 1/2/3指标 + 交互 =================
lmg <- function(dz, terms) {
  ## terms: 字符向量(模型右侧项名, 含交互)
  if (length(terms) == 0) return(setNames(numeric(0), NULL))
  n <- length(terms)
  if (n <= 6) {
    ## 精确: 全排列
    ord_list <- gtools_perm(terms)
  } else {
    ## 随机排列近似 (n>6, 如含交互项的7项模型): 500次随机排序
    set.seed(42)
    ord_list <- lapply(1:500, function(x) sample(terms))
  }
  inc <- setNames(rep(0, n), terms)
  for (ord in ord_list) {
    R2prev <- 0
    for (j in seq_along(ord)) {
      cols <- ord[1:j]
      f <- as.formula(paste("y ~", paste(cols, collapse="+")))
      R2 <- summary(lm(f, data=dz))$r.squared
      inc[ord[j]] <- inc[ord[j]] + (R2 - R2prev)
      R2prev <- R2
    }
  }
  inc <- inc / length(ord_list)
  inc
}
gtools_perm <- function(v) {
  ## 生成 v 的全排列(列表)
  n <- length(v)
  if (n <= 1) return(list(v))
  res <- list()
  for (k in seq_len(n)) {
    rest <- v[-k]
    if (length(rest) == 0) { res[[length(res)+1]] <- v[k]; next }
    sub <- Recall(rest)
    for (p in sub) res[[length(res)+1]] <- c(v[k], p)
  }
  res
}

mod_res <- NULL
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
    dz <- as.data.frame(scale(d)); names(dz) <- c(METRICS, "y")
    ## 交互项
    dz$int_em_mo  <- dz$n_edges_carried * dz$n_modules_covered
    dz$int_em_cpe <- dz$n_edges_carried * dz$mean_completeness_effective
    dz$int_mo_cpe <- dz$n_modules_covered * dz$mean_completeness_effective
    dz$int_all    <- dz$n_edges_carried * dz$n_modules_covered *
                     dz$mean_completeness_effective

    models <- list(
      m1 = METRICS[1],
      m2 = METRICS[2],
      m3 = METRICS[3],
      m12 = c(METRICS[1], METRICS[2]),
      m13 = c(METRICS[1], METRICS[3]),
      m23 = c(METRICS[2], METRICS[3]),
      m123 = METRICS,
      m123_i2 = c(METRICS, "int_em_mo","int_em_cpe","int_mo_cpe"),
      m123_i3 = c(METRICS, "int_em_mo","int_em_cpe","int_mo_cpe","int_all"))
    for (mn in names(models)) {
      terms <- models[[mn]]
      fit <- lm(as.formula(paste("y ~", paste(terms, collapse="+"))), data=dz)
      s <- summary(fit)
      R2 <- s$r.squared
      mp <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail=FALSE)
      npar <- length(terms)+1
      adjR2 <- 1-(1-R2)*(nrow(dz)-1)/(nrow(dz)-npar)
      L <- lmg(dz, terms)
      lsum <- sum(L)
      for (tt in terms) {
        mod_res <- rbind(mod_res, data.frame(
          pop=pop, trait=tr7[i], n=nrow(dz), model=mn,
          model_R2=round(R2,4), model_adjR2=round(adjR2,4), model_p=signif(mp,3),
          term=tt,
          lmg=round(L[[tt]],4),
          lmg_pct=round(L[[tt]]/R2*100,2),
          share_of_all=round(L[[tt]]/lsum*100,2)))
      }
    }
  }
}
write.table(mod_res, "qinben_progeny_pheno_nested_models_lmg_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: qinben_progeny_pheno_nested_models_lmg_txt\n")

## 模型水平汇总
ms <- mod_res %>% distinct(pop, trait, model, model_R2, model_adjR2, model_p)
summ_mod <- ms %>% group_by(trait, model) %>%
  summarise(n_models=n(),
            n_sig=sum(model_p<0.05),
            median_R2=round(median(model_R2),4),
            mean_R2=round(mean(model_R2),4),
            mean_adjR2=round(mean(model_adjR2),4),
            .groups="drop") %>% as.data.frame()
write.table(summ_mod, "qinben_progeny_pheno_nested_models_summary_txt",
            sep="\t", quote=FALSE, row.names=FALSE)

## 交互项贡献(仅含交互的模型 vs m123 差值)
comp <- ms %>% select(pop, trait, model, model_R2) %>%
  pivot_wider(names_from=model, values_from=model_R2) %>% as.data.frame()
comp$delta_vs_m123 <- round(comp$m123_i2 - comp$m123, 4)
comp$delta_all3 <- round(comp$m123_i3 - comp$m123, 4)
write.table(comp, "qinben_progeny_pheno_interaction_delta_txt",
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n=== 模型比较(按性状×模型平均R2) ===\n")
print(summ_mod %>% pivot_wider(names_from=model,
      values_from=c(n_sig,mean_R2)) %>% as.data.frame(), right=FALSE)

## ================= (c) 亲本互补程度 vs h2_epi / gain =================
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
    hh <- h2$h2_epi[h2$pop==pop & h2$trait==tr7[i]]
    hh <- ifelse(length(hh)==1, hh, NA)
    spv <- sp[[tr7[i]]][sp$group==pop]
    prog_mean <- mean(spv, na.rm=TRUE)
    a <- get_parent(pops$p0[r], pop, tr7[i])
    b <- get_parent(pops$p2[r], pop, tr7[i])
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
write.table(comp_all, "parent_complementarity_h2_gain_v2_txt",
            sep="\t", quote=FALSE, row.names=FALSE)

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
          r=NA, p_pearson=NA, rho=NA, p_spear=NA)); next
      }
      cp <- cor.test(dd[[m]], dd[[y]], method="pearson")
      ct <- cor.test(dd[[m]], dd[[y]], method="spearman", exact=FALSE)
      cors2 <- rbind(cors2, data.frame(trait=t, x=m, y=y, n=nrow(dd),
        r=round(cp$estimate,3), p_pearson=signif(cp$p.value,3),
        rho=round(ct$estimate,3), p_spear=signif(ct$p.value,3)))
    }
  }
}
write.table(cors2, "parent_complementarity_correlations_v2_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\n=== 亲本互补 vs h2/gain 相关 ===\n")
print(cors2, right=FALSE)

## 图: 互补指标 vs h2_epi / gain
pdf("parent_complementarity_scatter_v2.pdf", width=12, height=7)
for (t in tr7) {
  d <- comp_all %>% filter(trait==t) %>%
    select(pop, metric_module_comp, metric_edge_comp, metric_comp_union,
           h2_epi, gain) %>%
    pivot_longer(-c(pop, h2_epi, gain), names_to="metric", values_to="value") %>%
    as.data.frame()
  p <- ggplot(d, aes(x=value, y=h2_epi)) +
    geom_point(colour="#E41A1C", size=2.6) +
    geom_smooth(method="lm", se=TRUE, colour="grey30", linewidth=0.5) +
    geom_point(data=d, aes(x=value, y=gain), colour="#377EB8", size=2.6) +
    geom_smooth(data=d, aes(x=value, y=gain), method="lm", se=TRUE,
                colour="grey30", linewidth=0.5, linetype=2) +
    facet_wrap(~metric, scales="free_x", nrow=1) +
    labs(title=paste0(t, ": complementarity vs h2_epi (red, solid) / gain (blue, dashed)"),
         x="Complementarity metric", y="h2_epi / progeny gain") +
    theme_bw()
  print(p)
}
dev.off()
cat("saved: parent_complementarity_scatter_v2.pdf\nDONE\n")