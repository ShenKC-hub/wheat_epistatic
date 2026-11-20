# GWAS_epistasis_network_biparentalAVS.R
# 双亲群体 AVS：7个性状的上位效应网络分析（同 GWAS_epistasis_network_P*.R 口径，无协变量）
#   (1) 边分类 favorable/inferior/moderate（BiBj携带者 vs 非携带者 t检验 + D方向）
#   (2) PA 矩阵（rds）
#   (3) 有利上位网络：个体水平拓扑指标 + 模块对应表（fast greedy, >=5节点）
#   (4) Gephi 节点/边表
# 结果保存到 biparental/AVS/
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(Matrix)
library(igraph)
suppressMessages(library(data.table))
pop <- "AVS"
outdir <- file.path("biparental", pop)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
cat("===", pop, "start:", format(t0), "===\n")

## ---------- 数据提取（与 double_pop_effect_h2_<pop>.R 完全一致） ----------
df <- read.csv(paste0(pop, "_qtl_blue.csv"), header = FALSE)
df$V3 <- df$V3 * 1000000
all_block <- read.delim("all_block.bed", header = FALSE, comment.char = "#")
geno <- df[!is.na(df$V3), ]
df_odd  <- geno[seq(1, nrow(geno), by = 2), ]
df_even <- geno[(1:nrow(geno)) %% 2 == 0, ]
df_odd <- cbind(df_odd[, 1:3], end = df_even$V3, df_odd[, 4:ncol(df_odd)])

data <- NULL
for (k in 1:nrow(all_block)) {
  d2 <- df_odd[df_odd$V2 == all_block$V1[k] & df_odd$V3 <= all_block$V2[k] & df_odd$end >= all_block$V3[k], ]
  d3 <- df_odd[df_odd$V2 == all_block$V1[k] & df_odd$V3 >  all_block$V2[k] & df_odd$V3 <  all_block$V3[k], ]
  d4 <- rbind(d2, d3)
  data <- rbind(data, d4[1, ])
}
stopifnot(nrow(data) == nrow(all_block))
G_all <- t(data[, 5:ncol(data)])
colnames(G_all) <- all_block$V4
cat("samples:", ncol(df) - 3, "  blocks:", nrow(all_block), "\n")

MIN_MODULE <- 5
for (i in 1:7) {
  ti <- Sys.time()
  trait_raw <- t(df[i, 4:ncol(df)])
  snpA <- cbind(trait_raw, G_all)
  snpA <- na.omit(snpA)
  trait <- as.numeric(snpA[, 1])
  G <- snpA[, -1, drop = FALSE]
  n <- length(trait)
  ids <- sprintf("%s_S%03d", pop, seq_len(n))

  ## 效应量表（去重 A-B/B-A 重复对）
  eff <- read.delim(paste0("Epistasis_effect_", pop, "_P", i, "_txt"))
  eff$snp1 <- as.character(eff$snp1); eff$snp2 <- as.character(eff$snp2)
  key <- paste(pmin(eff$snp1, eff$snp2), pmax(eff$snp1, eff$snp2), sep = "|")
  eff <- eff[!duplicated(key), ]
  np <- nrow(eff)
  i1 <- match(eff$snp1, colnames(G))
  i2 <- match(eff$snp2, colnames(G))
  carr <- vector("list", np)
  cat("P", i, ": n =", n, " dedup pairs =", np, "\n")

  ## ---- 边分类：BiBj携带者 vs 非携带者（原始表型 Welch t 检验）+ D方向 ----
  testP <- rep(NA_real_, np); mc <- rep(NA_real_, np); mn <- rep(NA_real_, np)
  for (k in seq_len(np)) {
    mask <- G[, i1[k]] == 2 & G[, i2[k]] == 2
    rc <- which(mask)
    carr[[k]] <- rc
    r0 <- which(!mask)
    nc1 <- length(rc); nc0 <- length(r0)
    if (nc1 > 1 && nc0 > 1) {
      m1 <- mean(trait[rc]); v1 <- var(trait[rc])
      m0 <- mean(trait[r0]); v0 <- var(trait[r0])
      mc[k] <- m1; mn[k] <- m0
      den <- v1/nc1 + v0/nc0
      if (den > 0) {
        tst <- (m1 - m0) / sqrt(den)
        dfw <- den^2 / ((v1/nc1)^2/(nc1-1) + (v0/nc0)^2/(nc0-1))
        testP[k] <- 2 * pt(-abs(tst), dfw)
      }
    }
  }
  cls <- ifelse(!is.na(testP) & testP < 0.05,
                ifelse(eff$D > 0, "favorable", "inferior"), "moderate")
  edge_out <- data.frame(snp1 = eff$snp1, snp2 = eff$snp2, pvalue = eff$pvalue,
                         class = cls, testP = testP, D = eff$D, aaij = eff$aaij,
                         beta_epi = eff$beta_epi, partialR2_epi = eff$partialR2_epi,
                         n_BiBj = eff$n_BiBj, m_BiBj = eff$m_BiBj,
                         mean_non = mn, mean_carrier = mc)
  write.table(edge_out, file.path(outdir, paste0("Epistasis_edge_class_P", i, "_txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("P", i, "classes:", paste(names(table(cls)), table(cls), collapse = "  "), "\n")

  lens <- lengths(carr)
  ii <- unlist(carr, use.names = FALSE)
  jj <- rep(seq_len(np), times = lens)
  PA <- sparseMatrix(i = ii, j = jj, dims = c(n, np))
  colnames(PA) <- key[!duplicated(key)]
  rownames(PA) <- ids
  rm(carr, ii, jj)
  saveRDS(PA, file.path(outdir, paste0("Epistasis_PA_P", i, ".rds")))

  ## ---- 有利上位网络：个体拓扑 + 模块 ----
  fav <- which(cls == "favorable")
  nf <- length(fav)
  netstat_cols <- c("node_count","edge_count","mean_neighbor","graph_diameter",
                    "mean_component_size","component_count",
                    "max_component_node_count","max_component_edge_count",
                    "max_component_mean_neighbor",
                    "max_module_size","mean_module_size","module_count")
  stats <- as.data.frame(matrix(0, nrow = n, ncol = length(netstat_cols),
                                dimnames = list(ids, netstat_cols)))
  write.csv(edge_out[, c("snp1","snp2","class","pvalue","D","aaij","partialR2_epi")],
            file.path(outdir, paste0("gephi_edges_all_P", i, "_csv")), row.names = FALSE)

  if (nf > 0) {
    PAf <- PA[, fav, drop = FALSE]
    for (s in seq_len(n)) {
      cols <- which(PAf[s, ] != 0)
      if (length(cols) == 0) next
      g <- simplify(graph_from_data_frame(edge_out[fav[cols], c("snp1","snp2")], directed = FALSE))
      stats[s, "edge_count"] <- ecount(g)
      stats[s, "node_count"] <- vcount(g)
      stats[s, "mean_neighbor"] <- mean(degree(g))
      comp <- components(g)
      stats[s, "component_count"] <- comp$no
      stats[s, "mean_component_size"] <- mean(comp$csize)
      big <- which.max(comp$csize)
      gc <- induced_subgraph(g, which(comp$membership == big))
      stats[s, "max_component_node_count"] <- vcount(gc)
      stats[s, "max_component_edge_count"] <- ecount(gc)
      stats[s, "max_component_mean_neighbor"] <- mean(degree(gc))
      stats[s, "graph_diameter"] <- diameter(gc)
      mods <- list()
      for (cid in seq_len(comp$no)) {
        vv <- which(comp$membership == cid)
        gs <- induced_subgraph(g, vv)
        if (vcount(gs) < MIN_MODULE) next
        cm <- tryCatch(cluster_fast_greedy(gs), error = function(e) NULL)
        if (is.null(cm)) next
        for (m in communities(cm)) if (length(m) >= MIN_MODULE) mods[[length(mods)+1]] <- m
      }
      if (length(mods) > 0) {
        sz <- sapply(mods, length)
        stats[s, "module_count"] <- length(mods)
        stats[s, "max_module_size"] <- max(sz)
        stats[s, "mean_module_size"] <- mean(sz)
      }
    }
  }
  write.table(stats, file.path(outdir, paste0("Epistasis_netstats_P", i, "_txt")),
              sep = "\t", quote = FALSE, row.names = TRUE)

  ## ---- 全局有利网络模块对应表 + Gephi 表 ----
  g_fav <- if (nf > 0) simplify(graph_from_data_frame(edge_out[fav, c("snp1","snp2")], directed = FALSE)) else NULL
  mid <- 0
  if (nf > 0) {
    comp <- components(g_fav)
    memb <- setNames(rep(NA_character_, vcount(g_fav)), V(g_fav)$name)
    for (cid in seq_len(comp$no)) {
      vv <- V(g_fav)$name[comp$membership == cid]
      gs <- induced_subgraph(g_fav, vv)
      if (vcount(gs) < MIN_MODULE) next
      cm <- tryCatch(cluster_fast_greedy(gs), error = function(e) NULL)
      if (is.null(cm)) next
      for (m in communities(cm)) {
        if (length(m) >= MIN_MODULE) { mid <- mid + 1; memb[m] <- sprintf("%s_P%d_M%04d", pop, i, mid) }
      }
    }
    node_tab <- data.frame(HB = names(memb), module = memb, stringsAsFactors = FALSE)
    node_tab$module[is.na(node_tab$module)] <- "non_module(<5)"
    fkey <- paste(pmin(edge_out$snp1[fav], edge_out$snp2[fav]),
                  pmax(edge_out$snp1[fav], edge_out$snp2[fav]), sep = "|")
    deg <- table(fkey)
    node_tab$edge_count <- as.integer(deg[node_tab$HB])
    node_tab$edge_count[is.na(node_tab$edge_count)] <- 0L
    write.table(node_tab, file.path(outdir, paste0("Epistasis_node_module_P", i, "_txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    fm  <- match(edge_out$snp1[fav], node_tab$HB)
    fm2 <- match(edge_out$snp2[fav], node_tab$HB)
    edge_mod <- data.frame(snp1 = edge_out$snp1[fav], snp2 = edge_out$snp2[fav],
                           module_snp1 = node_tab$module[fm],
                           module_snp2 = node_tab$module[fm2],
                           D = edge_out$D[fav], aaij = edge_out$aaij[fav],
                           pvalue = edge_out$pvalue[fav])
    edge_mod$intra_module <- ifelse(edge_mod$module_snp1 == edge_mod$module_snp2,
                                    edge_mod$module_snp1, "cross_module")
    write.table(edge_mod, file.path(outdir, paste0("Epistasis_edge_module_P", i, "_txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    mod_sizes <- table(node_tab$module[node_tab$module != "non_module(<5)"])
    nmod <- length(mod_sizes)
    if (nmod > 0) {
      cc <- grDevices::col2rgb(grDevices::hsv(h = seq(0, 1, length.out = nmod + 1)[-1], s = 0.85, v = 0.95))
      col_df <- data.frame(module = names(mod_sizes), r = cc[1, ], g = cc[2, ], b = cc[3, ])
    } else col_df <- data.frame(module = character(0), r = integer(0), g = integer(0), b = integer(0))
    col_all <- rbind(col_df, data.frame(module = "non_module(<5)", r = 160, g = 160, b = 160))
    node_tab$module_size <- as.integer(mod_sizes[node_tab$module])
    node_tab$module_size[is.na(node_tab$module_size)] <- 0L
    node_tab <- cbind(node_tab, col_all[match(node_tab$module, col_all$module), c("r","g","b")])
    node_tab$Label <- node_tab$HB
    nodes_out <- node_tab[, c("HB","Label","module","module_size","edge_count","r","g","b")]
    names(nodes_out)[1] <- "Id"
    write.csv(nodes_out, file.path(outdir, paste0("gephi_nodes_P", i, "_csv")), row.names = FALSE)
    edges_out <- data.frame(Source = edge_out$snp1[fav], Target = edge_out$snp2[fav],
                            Type = "undirected", Weight = 1, module = edge_mod$intra_module)
    edges_out <- cbind(edges_out, col_all[match(edges_out$module, col_all$module), c("r","g","b")])
    write.csv(edges_out, file.path(outdir, paste0("gephi_edges_P", i, "_csv")), row.names = FALSE)
    sum_row <- data.frame(trait = paste0(pop, "_P", i), n_nodes = vcount(g_fav),
                          n_edges = ecount(g_fav), n_modules = nmod,
                          largest_module = if (nmod > 0) max(mod_sizes) else 0,
                          median_module = if (nmod > 0) round(median(mod_sizes), 1) else 0,
                          nodes_in_modules = sum(node_tab$module != "non_module(<5)"))
  } else {
    sum_row <- data.frame(trait = paste0(pop, "_P", i), n_nodes = 0, n_edges = 0,
                          n_modules = 0, largest_module = 0, median_module = 0, nodes_in_modules = 0)
  }
  write.table(sum_row, file.path(outdir, paste0("module_summary_P", i, "_txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("P", i, "done, elapsed:", round(as.numeric(difftime(Sys.time(), ti, units = "mins")), 1), "min\n")
}
cat("=== DONE", pop, "total:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ===\n")
