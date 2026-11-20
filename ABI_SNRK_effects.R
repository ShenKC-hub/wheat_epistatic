# ABI_SNRK_effects.R
# 在25个双亲群体中计算 ABI(3A_block_11916) 与 SNRK(5A_block_4075) 两个block的
# 加性效应与上位效应(交互项)。基因型与表型获取规则与 double_pop.R 一致。
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "ABI_SNRK"
dir.create(outdir, showWarnings = FALSE)

pops <- c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
          "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK")

## 两个目标block的坐标区间
targets <- data.frame(chr = c("3A", "5A"),
                      start = c(371606574, 78152983),
                      end   = c(371635959, 78164825),
                      block = c("3A_block_11916(ABI)", "5A_block_4075(SNRK)"),
                      stringsAsFactors = FALSE)

res <- NULL
for (p in pops) {
  df <- read.csv(paste0(p, "_qtl_blue.csv"), header = FALSE)
  df$V3 <- df$V3 * 1000000
  geno <- df[!is.na(df$V3), ]
  df_odd  <- geno[seq(1, nrow(geno), by = 2), ]
  df_even <- geno[(1:nrow(geno)) %% 2 == 0, ]
  df_odd <- cbind(df_odd[, 1:3], end = df_even$V3, df_odd[, 4:ncol(df_odd)])

  ## 每个目标区间取第一个重叠QTL行的基因型（与double_pop.R选点规则一致）
  gvec <- c(); okv <- c()
  for (k in 1:nrow(targets)) {
    d2 <- df_odd[df_odd$V2 == targets$chr[k] & df_odd$V3 <= targets$start[k] &
                 df_odd$end >= targets$end[k], ]
    d3 <- df_odd[df_odd$V2 == targets$chr[k] & df_odd$V3 >  targets$start[k] &
                 df_odd$V3 <  targets$end[k], ]
    d4 <- rbind(d2, d3)
    if (nrow(d4) == 0) { gvec[[k]] <- NULL; okv[k] <- FALSE; next }
    q <- d4[1, ]
    gvec[[k]] <- as.numeric(q[, 5:ncol(q)])
    okv[k] <- TRUE
  }
  if (!all(okv)) {
    for (i in 1:7) res <- rbind(res, data.frame(pop = p, trait = paste0("P", i),
      n = NA, note = "interval QTL not found"))
    next
  }
  g1 <- gvec[[1]]; g2 <- gvec[[2]]   # 3A(ABI), 5A(SNRK)
  poly1 <- length(unique(g1)) > 1; poly2 <- length(unique(g2)) > 1

  for (i in 1:7) {
    trait_raw <- t(df[i, 4:ncol(df)])
    y <- as.numeric(trait_raw)
    keep <- !is.na(y) & !is.na(g1) & !is.na(g2)
    y <- y[keep]; a <- g1[keep]; b <- g2[keep]
    n <- length(y)
    note <- paste0(ifelse(poly1, "", "3A_mono;"), ifelse(poly2, "", "5A_mono"))
    if (n < 8 || (!poly1 && !poly2)) {
      res <- rbind(res, data.frame(pop = p, trait = paste0("P", i), n = n, note = note,
        beta_3A = NA, se_3A = NA, p_3A = NA, alpha_3A = NA,
        beta_5A = NA, se_5A = NA, p_5A = NA, alpha_5A = NA,
        beta_epi = NA, se_epi = NA, p_epi = NA, R2 = NA))
      next
    }
    X <- cbind(1, a, b, a * b)
    fit <- lm.fit(X, y)
    np <- 4; dfr <- n - np
    cf <- fit$coefficients
    if (any(is.na(cf))) {
      res <- rbind(res, data.frame(pop = p, trait = paste0("P", i), n = n, note = "rank_deficient",
        beta_3A = NA, se_3A = NA, p_3A = NA, alpha_3A = NA,
        beta_5A = NA, se_5A = NA, p_5A = NA, alpha_5A = NA,
        beta_epi = NA, se_epi = NA, p_epi = NA, R2 = NA))
      next
    }
    se <- sqrt(diag(chol2inv(qr.R(fit$qr))) * sum(fit$residuals^2) / dfr)
    tv <- cf / se
    pv <- 2 * pt(-abs(tv), dfr)
    r2 <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
    res <- rbind(res, data.frame(pop = p, trait = paste0("P", i), n = n, note = note,
      beta_3A = cf[2], se_3A = se[2], p_3A = pv[2], alpha_3A = 2 * cf[2],
      beta_5A = cf[3], se_5A = se[3], p_5A = pv[3], alpha_5A = 2 * cf[3],
      beta_epi = cf[4], se_epi = se[4], p_epi = pv[4], R2 = round(r2, 4)))
  }
  cat(p, "done\n")
}
res$trait <- factor(res$trait, levels = paste0("P", 1:7))
res <- res[order(res$pop, res$trait), ]
write.table(res, file.path(outdir, "ABI_SNRK_effects_25pops_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 汇总 ----
cat("\n===== 汇总 =====\n")
cat("群体×性状组合:", nrow(res), "\n")
sig_epi <- res[!is.na(res$p_epi) & res$p_epi < 0.05, ]
cat("上位效应显著(P<0.05):", nrow(sig_epi), "个组合\n")
if (nrow(sig_epi) > 0) print(sig_epi[, c("pop","trait","n","beta_epi","se_epi","p_epi","R2")], row.names = FALSE, digits = 3)
sig_a1 <- res[!is.na(res$p_3A) & res$p_3A < 0.05, ]
cat("\nABI(3A)加性效应显著(P<0.05):", nrow(sig_a1), "个组合\n")
if (nrow(sig_a1) > 0) print(sig_a1[, c("pop","trait","beta_3A","p_3A")], row.names = FALSE, digits = 3)
sig_a2 <- res[!is.na(res$p_5A) & res$p_5A < 0.05, ]
cat("\nSNRK(5A)加性效应显著(P<0.05):", nrow(sig_a2), "个组合\n")
if (nrow(sig_a2) > 0) print(sig_a2[, c("pop","trait","beta_5A","p_5A")], row.names = FALSE, digits = 3)
cat("DONE\n")
