# GWAS_epistasis_network_P1.R
# 性状 P1：上位效应网络分析
#   (1) 样本×互作 presence-absence 矩阵（等价于 HB-HB 携带矩阵）
#   (2) 有利/不利/中等上位分类（BiBj携带者 vs 非携带者 t检验 + D方向）
#   (3) 有利上位网络：全局 + 个体水平 igraph 拓扑指标与模块（fast greedy, >=5节点）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(Matrix)
library(igraph)
i <- 1
t0 <- Sys.time()
cat("=== Network P", i, "start:", format(t0), "===\n")

## ---------- 数据读取（与 GWAS_epistasis_effect_h2_P*.R 一致，附加样本ID） ----------
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
pca <- read.delim("pca.txt", header=FALSE, comment.char="#", sep="\t")
phe <- merge(phe, pca, by="V1")
colnames(genotype) <- genotype[1, ]
genotype <- genotype[-1, ]
phe <- merge(phe, genotype, by.x="V1", by.y="ID")

snpA <- na.omit(phe[, c(1, i+1, 9:ncol(phe))])
ids    <- as.character(snpA$V1)
trait  <- as.numeric(snpA[, 2])
G <- as.matrix(sapply(snpA[, -(1:5)], as.numeric))
n <- length(trait)
cat("samples:", n, " SNPs:", ncol(G), "\n")

pairs <- read.delim(paste0("gwas_Epistasis_P", i, "_txt"))
pairs$snp1 <- as.character(pairs$snp1)
pairs$snp2 <- as.character(pairs$snp2)
eff <- read.delim(paste0("gwas_Epistasis_effect_P", i, "_txt"))
eff$snp1 <- as.character(eff$snp1)
eff$snp2 <- as.character(eff$snp2)
pairs$key <- paste(pairs$snp1, pairs$snp2)
eff$key   <- paste(eff$snp1, eff$snp2)
eff <- eff[match(pairs$key, eff$key), ]
stopifnot(identical(eff$key, pairs$key))
np <- nrow(pairs)
cat("pairs:", np, "\n")

i1 <- match(pairs$snp1, colnames(G))
i2 <- match(pairs$snp2, colnames(G))

## ---------- 第一步：presence-absence 矩阵 + 边分类 ----------
PA <- Matrix(0L, nrow=n, ncol=np, sparse=TRUE)   # 1 = 样本携带该互作(BiBj)
cls   <- character(np)
testP <- rep(NA_real_, np)
mean_carrier <- rep(NA_real_, np)
mean_non     <- rep(NA_real_, np)

for (k in seq_len(np)) {
  mask <- G[, i1[k]] == 2 & G[, i2[k]] == 2
  rc <- which(mask)
  PA[rc, k] <- 1L
  r0 <- which(!mask)
  nc1 <- length(rc); nc0 <- length(r0)
  if (nc1 > 1 && nc0 > 1) {
    m1 <- mean(trait[rc]); v1 <- var(trait[rc])
    m0 <- mean(trait[r0]); v0 <- var(trait[r0])
    mean_carrier[k] <- m1; mean_non[k] <- m0
    den <- v1/nc1 + v0/nc0
    if (den > 0) {
      tst <- (m1 - m0) / sqrt(den)
      dfw <- den^2 / ((v1/nc1)^2/(nc1-1) + (v0/nc0)^2/(nc0-1))
      testP[k] <- 2 * pt(-abs(tst), dfw)
    }
  }
  Dk <- eff$D[k]
  cls[k] <- ifelse(!is.na(testP[k]) && testP[k] < 0.05,
                   ifelse(Dk > 0, "favorable", "inferior"),
                   "moderate")
}
edge_tab <- data.frame(snp1=pairs$snp1, snp2=pairs$snp2,
                       pvalue=pairs$pvalue, class=cls, testP=testP,
                       D=eff$D, aaij=eff$aaij, beta_epi=eff$beta_epi,
                       partialR2_epi=eff$partialR2_epi,
                       n_BiBj=eff$n_BiBj, m_BiBj=eff$m_BiBj,
                       mean_non=mean_non, mean_carrier=mean_carrier)
write.table(edge_tab, paste0("gwas_epistasis_edge_class_P", i, "_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("classes:", paste(names(table(cls)), table(cls), collapse="  "), "\n")

## 保存 presence-absence 稀疏矩阵（行=样本，列=显著互作）
colnames(PA) <- pairs$key
rownames(PA) <- ids
saveRDS(PA, paste0("gwas_epistasis_PA_P", i, ".rds"))
cat("PA matrix saved:", dim(PA)[1], "x", dim(PA)[2],
    " elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min\n")

## ---------- 第二步：有利上位网络（igraph） ----------
t2 <- Sys.time()
fav <- which(cls == "favorable")
nf <- length(fav)
cat("favorable edges:", nf, "\n")

## Gephi 边表（全部显著边，含类别属性）
write.csv(edge_tab[, c("snp1","snp2","class","pvalue","D","aaij","partialR2_epi")],
          paste0("gwas_epistasis_network_edges_P", i, "_csv"), row.names=FALSE)

netstat_cols <- c("node_count","edge_count","mean_neighbor","graph_diameter",
                  "mean_component_size","component_count",
                  "max_component_node_count","max_component_edge_count",
                  "max_component_mean_neighbor",
                  "max_module_size","mean_module_size","module_count")
stats <- as.data.frame(matrix(0, nrow=n, ncol=length(netstat_cols),
                              dimnames=list(ids, netstat_cols)))

## 全局有利网络
g_global <- graph_from_data_frame(edge_tab[fav, c("snp1","snp2")], directed=FALSE)
g_global <- simplify(g_global)
cat("global network: nodes =", vcount(g_global), " edges =", ecount(g_global), "\n")

## 个体水平：每个样本携带的有利互作子网络
if (nf > 0) {
  PAf <- PA[, fav, drop=FALSE]
  for (s in seq_len(n)) {
    cols <- which(PAf[s, ] != 0)
    if (length(cols) == 0) next
    e <- edge_tab[fav[cols], c("snp1","snp2")]
    g <- graph_from_data_frame(e, directed=FALSE)
    g <- simplify(g)
    stats[s, "edge_count"] <- ecount(g)
    stats[s, "node_count"] <- vcount(g)
    stats[s, "mean_neighbor"] <- mean(degree(g))
    comp <- components(g)
    stats[s, "component_count"] <- comp$no
    stats[s, "mean_component_size"] <- mean(comp$csize)
    big <- which.max(comp$csize)
    vid <- which(comp$membership == big)
    gc <- induced_subgraph(g, vid)
    stats[s, "max_component_node_count"] <- vcount(gc)
    stats[s, "max_component_edge_count"] <- ecount(gc)
    stats[s, "max_component_mean_neighbor"] <- mean(degree(gc))
    stats[s, "graph_diameter"] <- diameter(gc)
    ## 模块检测：对每个连通分量跑 fast greedy，合并所有 >=5 节点的模块
    mods <- list()
    for (cid in seq_len(comp$no)) {
      vv <- which(comp$membership == cid)
      gs <- induced_subgraph(g, vv)
      if (vcount(gs) < 5) next
      cm <- tryCatch(cluster_fast_greedy(gs), error=function(e) NULL)
      if (is.null(cm)) next
      for (m in communities(cm)) if (length(m) >= 5) mods[[length(mods)+1]] <- m
    }
    if (length(mods) > 0) {
      sz <- sapply(mods, length)
      stats[s, "module_count"] <- length(mods)
      stats[s, "max_module_size"] <- max(sz)
      stats[s, "mean_module_size"] <- mean(sz)
    }
  }
}
write.table(stats, paste0("gwas_epistasis_netstats_P", i, "_txt"),
            sep="\t", quote=FALSE, row.names=TRUE)
cat("=== DONE Network P", i, "total:",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min ===\n")
