# GWAS_epistasis_effect_h2_P1.R
# 性状 P1：(1) 显著SNP对的加性效应、上位效应、四基因型组合均值
#              (2) sommer 混合模型遗传力分解（加性 + 上位）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(sommer)
i <- 1
t0 <- Sys.time()
cat("=== P", i, "start:", format(t0), "===\n")

## ---------- 数据读取（与 GWAS_epistasis_P*.R 完全一致） ----------
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])

phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
pca <- read.delim("pca.txt", header=FALSE, comment.char="#", sep="\t")
phe <- merge(phe, pca, by="V1")

colnames(genotype) <- genotype[1, ]
genotype <- genotype[-1, ]
phe <- merge(phe, genotype, by.x="V1", by.y="ID")

snpA <- phe[, c(i+1, 9:ncol(phe))]
snpA <- na.omit(snpA)
trait     <- as.numeric(snpA[, 1])
covariate <- as.matrix(sapply(snpA[, 2:4], as.numeric))   # 3 PCs 群体结构矩阵
G <- as.matrix(sapply(snpA[, -c(1:4)], as.numeric))       # 0/2 编码
cat("samples:", nrow(G), "  SNPs:", ncol(G), "\n")

pairs <- read.delim(paste0("gwas_Epistasis_P", i, "_txt"))
pairs$snp1 <- as.character(pairs$snp1)
pairs$snp2 <- as.character(pairs$snp2)
ok <- pairs$snp1 %in% colnames(G) & pairs$snp2 %in% colnames(G)
if (any(!ok)) cat("WARNING: dropping", sum(!ok), "pairs with SNP id not found\n")
pairs <- pairs[ok, ]
cat("pairs:", nrow(pairs), "\n")

## ---------- 第一步：逐对效应量 ----------
## 四基因型组合均值使用协变量校正后的表型: y_adj = 残差(y~PCs) + 总体均值
Xcov  <- cbind(1, covariate)
fitc  <- lm.fit(Xcov, trait)
y_adj <- fitc$residuals + mean(trait)
meanY <- mean(trait)

nc  <- ncol(Xcov)          # 4: 截距 + 3 PC
n   <- length(trait)
np  <- nc + 3              # 全模型参数个数
CODES <- c(6, 2, 4, 0)     # a+2b: 6=BiBj, 2=BiAj, 4=AiBj, 0=AiAj

COLS <- c("snp1","snp2",
          "n_BiBj","n_BiAj","n_AiBj","n_AiAj","m_BiBj","m_BiAj","m_AiBj","m_AiAj",
          "aaij","D","rel_D",
          "beta_epi","se_epi","p_epi","partialR2_epi",
          "beta_snp1","se_snp1","p_snp1","alpha_snp1",
          "beta_snp2","se_snp2","p_snp2","alpha_snp2",
          "pvalue")
out <- matrix(NA_real_, nrow=nrow(pairs), ncol=length(COLS)-3)

t1 <- Sys.time()
for (k in seq_len(nrow(pairs))) {
  a <- G[, pairs$snp1[k]]
  b <- G[, pairs$snp2[k]]

  ## 四组合均值（协变量校正表型）与样本数: 顺序 BiBj, BiAj, AiBj, AiAj
  ns <- integer(4); ms <- rep(NA_real_, 4)
  code <- a + 2*b
  for (g in 1:4) {
    sel <- code == CODES[g]
    ns[g] <- sum(sel)
    if (ns[g] > 0) ms[g] <- mean(y_adj[sel])
  }

  ## 回归: y ~ cov + snp1 + snp2 + snp1:snp2
  X   <- cbind(Xcov, a, b, a*b)
  fit <- lm.fit(X, trait)
  cf  <- fit$coefficients[(nc+1):(nc+3)]
  if (any(is.na(cf))) next
  dfr  <- n - np
  s2   <- sum(fit$residuals^2) / dfr
  se   <- sqrt(s2 * diag(chol2inv(qr.R(fit$qr)))[(nc+1):(nc+3)])
  tv   <- cf / se
  pv   <- 2 * pt(-abs(tv), dfr)
  pR2  <- tv[3]^2 / (tv[3]^2 + dfr)   # 交互项偏R²

  aaij <- (ms[1] - ms[2] - ms[3] + ms[4]) / 4   # 经典上位效应
  D    <- 4 * aaij                               # BiBj组合超加性偏移量
  relD <- ifelse(meanY != 0, D / meanY * 100, NA)

  out[k, ] <- c(ns, ms, aaij, D, relD,
                cf[3], se[3], pv[3], pR2,
                cf[1], se[1], pv[1], 2*cf[1],
                cf[2], se[2], pv[2], 2*cf[2])
}
res <- data.frame(snp1=pairs$snp1, snp2=pairs$snp2, out, pvalue=pairs$pvalue)
colnames(res) <- COLS

write.table(res, paste0("gwas_Epistasis_effect_P", i, "_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("effect done:", nrow(res), "pairs, elapsed:",
    round(as.numeric(difftime(Sys.time(), t1, units="mins")), 1), "min\n")

## ---------- 第二步：遗传力分解（sommer, GRM 方差组分） ----------
t2 <- Sys.time()
snps <- unique(c(pairs$snp1, pairs$snp2))
Gs   <- scale(G[, snps, drop=FALSE])

K.add <- tcrossprod(Gs) / (ncol(Gs) - 1)   # 加性 GRM（显著位点）
K.epi <- K.add * K.add                     # 上位 GRM（Hadamard 积）
rownames(K.add) <- colnames(K.add) <- as.character(seq_len(nrow(K.add)))
rownames(K.epi) <- colnames(K.epi) <- rownames(K.add)

dt <- data.frame(y = trait,
                 PC1 = covariate[, 1], PC2 = covariate[, 2], PC3 = covariate[, 3])
dt$user  <- factor(seq_len(n), levels = as.character(seq_len(n)))
dt$units <- factor(seq_len(n))

mod <- mmer(y ~ PC1 + PC2 + PC3,
            random = ~ vsr(user, Gu=K.add) + vsr(user, Gu=K.epi),
            data = dt)   # 残差由 sommer 默认 rcov=~units 估计
cat("convergence:", mod$convergence, "\n")
print(mod$sigma_scaled)
sig <- unlist(mod$sigma)
stopifnot(length(sig) == 3)
Va   <- sig[1]    # K.add  -> 加性方差
Vepi <- sig[2]    # K.epi  -> 上位方差
Ve   <- sig[3]    # 残差
Vp   <- Va + Vepi + Ve

h2 <- data.frame(trait = paste0("P", i),
                 n_pairs = nrow(pairs), n_snps = length(snps),
                 Va = Va, V_epi = Vepi, Ve = Ve, Vp = Vp,
                 h2_add = Va / Vp, h2_epi = Vepi / Vp,
                 h2_total = (Va + Vepi) / Vp)
write.table(h2, paste0("gwas_epistasis_h2_P", i, "_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("h2 done, elapsed:", round(as.numeric(difftime(Sys.time(), t2, units="mins")), 1), "min\n")
cat("=== DONE P", i, "total:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min ===\n")
