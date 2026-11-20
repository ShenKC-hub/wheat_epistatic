# GWAS_epistasis_network_inferior_P1.R
# 性状 P1：不利上位效应（inferior）网络分析，结果保存到 inferior/ 文件夹
#   (1) 样本×互作 presence-absence 矩阵（rds + txt）
#   (2) 个体水平 inferior 网络拓扑指标（与 favorable 版 netstats 同口径）
#   (3) 全局 inferior 网络：节点/边 -> 模块对应表（fast greedy, >=5节点）
#   (4) Gephi 节点/边表（按模块着色）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(Matrix)
library(igraph)
suppressMessages(library(data.table))
outdir <- "inferior"
i <- 1
t0 <- Sys.time()
cat("=== Inferior Network P", i, "start:", format(t0), "===\n")

## ---------- 数据读取（与 favorable 网络脚本一致） ----------
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

edge_tab <- read.delim(paste0("gwas_epistasis_edge_class_P", i, "_txt"))
inf <- edge_tab[edge_tab$class == "inferior", ]
np <- nrow(inf)
i1 <- match(inf$snp1, colnames(G))
i2 <- match(inf$snp2, colnames(G))
cat("samples:", n, " inferior edges:", np, "\n")
if (np == 0) {
  cat("=== DONE (no inferior edges) ===\n"); quit(save="no")
}

## ---------- (1) presence-absence 矩阵 ----------
PA <- Matrix(0L, nrow=n, ncol=np, sparse=TRUE)
for (k in seq_len(np)) PA[which(G[, i1[k]] == 2 & G[, i2[k]] == 2), k] <- 1L
colnames(PA) <- paste(inf$snp1, inf$snp2)
rownames(PA) <- ids
saveRDS(PA, file.path(outdir, paste0("gwas_epistasis_PA_inferior_P", i, ".rds")))
fwrite(as.data.frame(as.matrix(PA)),
       file.path(outdir, paste0("gwas_epistasis_PA_inferior_P", i, "_txt")),
       sep = "\t", quote = FALSE, row.names = TRUE)
cat("PA saved:", dim(PA)[1], "x", dim(PA)[2], "\n")

## ---------- (2) 个体水平 inferior 网络拓扑指标 ----------
netstat_cols <- c("node_count","edge_count","mean_neighbor","graph_diameter",
                  "mean_component_size","component_count",
                  "max_component_node_count","max_component_edge_count",
                  "max_component_mean_neighbor",
                  "max_module_size","mean_module_size","module_count")
stats <- as.data.frame(matrix(0, nrow=n, ncol=length(netstat_cols),
                              dimnames=list(ids, netstat_cols)))
g_global <- simplify(graph_from_data_frame(inf[, c("snp1","snp2")], directed=FALSE))
cat("global inferior network: nodes =", vcount(g_global),
    " edges =", ecount(g_global), "\n")

for (s in seq_len(n)) {
  cols <- which(PA[s, ] != 0)
  if (length(cols) == 0) next
  g <- simplify(graph_from_data_frame(inf[cols, c("snp1","snp2")], directed=FALSE))
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
write.table(stats, file.path(outdir, paste0("gwas_epistasis_netstats_inferior_P", i, "_txt")),
            sep="\t", quote=FALSE, row.names=TRUE)
cat("netstats done, elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min\n")

## ---------- (3) 全局 inferior 网络模块对应表 ----------
MIN_MODULE <- 5
comp <- components(g_global)
memb <- setNames(rep(NA_character_, vcount(g_global)), V(g_global)$name)
mid <- 0
for (cid in seq_len(comp$no)) {
  vv <- V(g_global)$name[comp$membership == cid]
  gs <- induced_subgraph(g_global, vv)
  if (vcount(gs) < MIN_MODULE) next
  cm <- tryCatch(cluster_fast_greedy(gs), error=function(e) NULL)
  if (is.null(cm)) next
  for (m in communities(cm)) {
    if (length(m) >= MIN_MODULE) {
      mid <- mid + 1
      memb[m] <- sprintf("P%d_inferior_M%04d", i, mid)
    }
  }
}
node_tab <- data.frame(HB = names(memb), module = memb, stringsAsFactors=FALSE)
node_tab$module[is.na(node_tab$module)] <- "non_module(<5)"
fkey <- apply(inf[, c("snp1","snp2")], 1, function(x) paste(sort(x), collapse="|"))
fav_u <- inf[!duplicated(fkey), ]
deg <- table(c(fav_u$snp1, fav_u$snp2))
node_tab$edge_count <- as.integer(deg[node_tab$HB])
node_tab$edge_count[is.na(node_tab$edge_count)] <- 0L
write.table(node_tab, file.path(outdir, paste0("gwas_epistasis_node_module_inferior_P", i, "_txt")),
            sep="\t", quote=FALSE, row.names=FALSE)

fm  <- match(inf$snp1, node_tab$HB)
fm2 <- match(inf$snp2, node_tab$HB)
edge_mod <- data.frame(snp1 = inf$snp1, snp2 = inf$snp2,
                       module_snp1 = node_tab$module[fm],
                       module_snp2 = node_tab$module[fm2],
                       D = inf$D, aaij = inf$aaij, pvalue = inf$pvalue,
                       stringsAsFactors=FALSE)
edge_mod$intra_module <- ifelse(edge_mod$module_snp1 == edge_mod$module_snp2,
                                edge_mod$module_snp1, "cross_module")
write.table(edge_mod, file.path(outdir, paste0("gwas_epistasis_edge_module_inferior_P", i, "_txt")),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("modules:", mid, " nodes in modules:", sum(node_tab$module != "non_module(<5)"), "\n")

## ---------- (4) Gephi 节点/边表（按模块着色） ----------
mod_sizes <- table(node_tab$module[node_tab$module != "non_module(<5)"])
nmod <- length(mod_sizes)
if (nmod > 0) {
  cc <- grDevices::col2rgb(grDevices::hsv(h = seq(0, 1, length.out = nmod + 1)[-1], s = 0.85, v = 0.95))
  col_df <- data.frame(module = names(mod_sizes), r = cc[1, ], g = cc[2, ], b = cc[3, ])
} else col_df <- data.frame(module=character(0), r=integer(0), g=integer(0), b=integer(0))
grey <- data.frame(module = "non_module(<5)", r = 160, g = 160, b = 160)
col_all <- rbind(col_df, grey)
node_tab$module_size <- as.integer(mod_sizes[node_tab$module])
node_tab$module_size[is.na(node_tab$module_size)] <- 0L
node_tab <- cbind(node_tab, col_all[match(node_tab$module, col_all$module), c("r","g","b")])
node_tab$Label <- node_tab$HB
nodes_out <- node_tab[, c("HB","Label","module","module_size","edge_count","r","g","b")]
names(nodes_out)[1] <- "Id"
write.csv(nodes_out, file.path(outdir, paste0("gephi_nodes_inferior_P", i, "_csv")),
          row.names = FALSE, quote = TRUE)
edges_out <- data.frame(Source = inf$snp1, Target = inf$snp2, Type = "undirected",
                        Weight = 1, module = edge_mod$intra_module)
edges_out <- cbind(edges_out, col_all[match(edges_out$module, col_all$module), c("r","g","b")])
write.csv(edges_out, file.path(outdir, paste0("gephi_edges_inferior_P", i, "_csv")),
          row.names = FALSE, quote = TRUE)

## 汇总行
sum_row <- data.frame(trait = paste0("P", i), n_nodes = vcount(g_global),
                      n_edges = ecount(g_global), n_modules = mid,
                      largest_module = if (nmod > 0) max(mod_sizes) else 0,
                      median_module = if (nmod > 0) round(median(mod_sizes), 1) else 0,
                      nodes_in_modules = sum(node_tab$module != "non_module(<5)"))
fwrite(sum_row, file.path(outdir, paste0("module_summary_inferior_P", i, "_txt")),
       sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== DONE Inferior Network P", i, "total:",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min ===\n")
