# double_pop_effect_h2AVS.R
# 双亲群体 AVS：7个性状的显著SNP对效应量（无协变量）+ sommer遗传力分解
# 数据提取与 double_pop.R 完全一致；忽略群体结构（不加PCA协变量）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages(library(sommer))
pop <- "AVS"
t0 <- Sys.time()
cat("===", pop, "start:", format(t0), "===\n")

## ---------- 数据提取（与 double_pop.R 完全一致） ----------
df <- read.csv(paste0(pop, "_qtl_blue.csv"), header=FALSE)
df$V3 <- df$V3 * 1000000

all_block <- read.delim("all_block.bed", header=FALSE, comment.char="#")

geno <- df[!is.na(df$V3), ]
df_odd  <- geno[seq(1, nrow(geno), by=2), ]
df_even <- geno[(1:nrow(geno)) %% 2 == 0, ]
df_odd <- cbind(df_odd[, 1:3], end=df_even$V3, df_odd[, 4:ncol(df_odd)])  # 等价于 mutate(end, .after=3)

data <- NULL
for (k in 1:nrow(all_block)) {
  d2 <- df_odd[df_odd$V2 == all_block$V1[k] & df_odd$V3 <= all_block$V2[k] & df_odd$end >= all_block$V3[k], ]
  d3 <- df_odd[df_odd$V2 == all_block$V1[k] & df_odd$V3 >  all_block$V2[k] & df_odd$V3 <  all_block$V3[k], ]
  d4 <- rbind(d2, d3)
  data <- rbind(data, d4[1, ])
}
stopifnot(nrow(data) == nrow(all_block))
G_all <- t(data[, 5:ncol(data)])     # 与原脚本每轮 t(data[,5:ncol(data)]) 等价
colnames(G_all) <- all_block$V4
cat("samples:", ncol(df) - 3, "  blocks:", nrow(all_block), "\n")

## ---------- 逐性状处理 ----------
for (i in 1:7) {
  ti <- Sys.time()
  trait_raw <- t(df[i, 4:ncol(df)])
  snpA <- cbind(trait_raw, G_all)
  snpA <- na.omit(snpA)
  trait <- as.numeric(snpA[, 1])
  G <- snpA[, -1, drop=FALSE]
  n <- length(trait)
  cat("P", i, ": n =", n, "\n")

  pairs <- read.delim(paste0("Epistasis_", pop, "_P", i, "_txt"))
  pairs$snp1 <- as.character(pairs$snp1)
  pairs$snp2 <- as.character(pairs$snp2)
  ok <- pairs$snp1 %in% colnames(G) & pairs$snp2 %in% colnames(G)
  if (any(!ok)) cat("WARNING: dropping", sum(!ok), "pairs with SNP id not found\n")
  pairs <- pairs[ok, ]
  cat("P", i, ": pairs =", nrow(pairs), "\n")

  ## ---- 第一步：逐对效应量（无协变量: y ~ snp1 + snp2 + snp1:snp2）----
  ## 四组合均值直接用原始表型（无群体结构可校正）
  COLS <- c("snp1","snp2",
            "n_BiBj","n_BiAj","n_AiBj","n_AiAj","m_BiBj","m_BiAj","m_AiBj","m_AiAj",
            "aaij","D","rel_D",
            "beta_epi","se_epi","p_epi","partialR2_epi",
            "beta_snp1","se_snp1","p_snp1","alpha_snp1",
            "beta_snp2","se_snp2","p_snp2","alpha_snp2",
            "pvalue")
  out <- matrix(NA_real_, nrow=nrow(pairs), ncol=length(COLS)-3)
  Xcov <- matrix(1, n, 1)              # 仅截距
  nc <- 1
  np  <- nc + 3
  meanY <- mean(trait)
  CODES <- c(6, 2, 4, 0)               # a+2b: 6=BiBj, 2=BiAj, 4=AiBj, 0=AiAj
  idx1 <- match(pairs$snp1, colnames(G))
  idx2 <- match(pairs$snp2, colnames(G))

  for (k in seq_len(nrow(pairs))) {
    a <- G[, idx1[k]]
    b <- G[, idx2[k]]

    ns <- integer(4); ms <- rep(NA_real_, 4)
    code <- a + 2*b
    for (g in 1:4) {
      sel <- code == CODES[g]
      ns[g] <- sum(sel)
      if (ns[g] > 0) ms[g] <- mean(trait[sel])
    }

    X   <- cbind(Xcov, a, b, a*b)
    fit <- lm.fit(X, trait)
    cf  <- fit$coefficients[(nc+1):(nc+3)]
    if (any(is.na(cf))) next
    dfr  <- n - np
    s2   <- sum(fit$residuals^2) / dfr
    se   <- sqrt(s2 * diag(chol2inv(qr.R(fit$qr)))[(nc+1):(nc+3)])
    tv   <- cf / se
    pv   <- 2 * pt(-abs(tv), dfr)
    pR2  <- tv[3]^2 / (tv[3]^2 + dfr)

    aaij <- (ms[1] - ms[2] - ms[3] + ms[4]) / 4
    D    <- 4 * aaij
    relD <- ifelse(meanY != 0, D / meanY * 100, NA)

    out[k, ] <- c(ns, ms, aaij, D, relD,
                  cf[3], se[3], pv[3], pR2,
                  cf[1], se[1], pv[1], 2*cf[1],
                  cf[2], se[2], pv[2], 2*cf[2])
  }
  res <- data.frame(snp1=pairs$snp1, snp2=pairs$snp2, out, pvalue=pairs$pvalue)
  colnames(res) <- COLS
  write.table(res, paste0("Epistasis_effect_", pop, "_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("P", i, "effect done, elapsed:",
      round(as.numeric(difftime(Sys.time(), ti, units="mins")), 1), "min\n")

  ## ---- 第二步：遗传力分解（sommer, 无协变量: y ~ 1）----
  ti2 <- Sys.time()
  h2 <- data.frame(trait=paste0(pop, "_P", i),
                   n_pairs=nrow(pairs), n_snps=NA_integer_,
                   Va=NA_real_, V_epi=NA_real_, Ve=NA_real_, Vp=NA_real_,
                   h2_add=NA_real_, h2_epi=NA_real_, h2_total=NA_real_,
                   note="no significant pairs")
  if (nrow(pairs) > 0) {
    snps <- unique(c(pairs$snp1, pairs$snp2))
    Gs   <- G[, snps, drop=FALSE]
    sds  <- apply(Gs, 2, sd)
    if (any(sds == 0)) {
      cat("WARNING: dropping", sum(sds == 0), "zero-variance SNPs for GRM\n")
      Gs <- Gs[, sds > 0, drop=FALSE]
    }
    h2$n_snps <- ncol(Gs)
    if (ncol(Gs) > 1) {
      h2r <- tryCatch({
        Gsc  <- scale(Gs)
        K.add <- tcrossprod(Gsc) / (ncol(Gsc) - 1)
        K.epi <- K.add * K.add
        rownames(K.add) <- colnames(K.add) <- as.character(seq_len(nrow(K.add)))
        rownames(K.epi) <- colnames(K.epi) <- rownames(K.add)
        dt <- data.frame(y = trait)
        dt$user  <- factor(seq_len(n), levels = as.character(seq_len(n)))
        mod <- mmer(y ~ 1,
                    random = ~ vsr(user, Gu=K.add) + vsr(user, Gu=K.epi),
                    data = dt)   # 残差由 sommer 默认 rcov=~units 估计
        cat("P", i, "convergence:", mod$convergence, "\n")
        print(mod$sigma_scaled)
        sig <- unlist(mod$sigma)
        stopifnot(length(sig) == 3)
        Va <- sig[1]; Vepi <- sig[2]; Ve <- sig[3]
        Vp <- Va + Vepi + Ve
        data.frame(Va=Va, V_epi=Vepi, Ve=Ve, Vp=Vp,
                   h2_add=Va/Vp, h2_epi=Vepi/Vp, h2_total=(Va+Vepi)/Vp,
                   note="OK")
      }, error = function(e) { cat("P", i, "mmer ERROR:", conditionMessage(e), "\n")
                              data.frame(Va=NA, V_epi=NA, Ve=NA, Vp=NA,
                                         h2_add=NA, h2_epi=NA, h2_total=NA,
                                         note=paste("mmer error:", conditionMessage(e))) })
      cols <- c("Va","V_epi","Ve","Vp","h2_add","h2_epi","h2_total","note")
      h2[, cols] <- h2r[, cols]
    }
  }
  write.table(h2, paste0("Epistasis_h2_", pop, "_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("P", i, "h2 done, elapsed:",
      round(as.numeric(difftime(Sys.time(), ti2, units="mins")), 1), "min\n")
}
cat("=== DONE", pop, "total:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min ===\n")
