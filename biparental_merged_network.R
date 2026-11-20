# biparental_merged_network.R
# 23个双亲群体（排除BBC、BLB）的全部上位效应对整合去重（不区分有利/不利）
# 构建统一网络：模块对应表（fast greedy, >=5节点）+ Gephi 表
# 结果保存到 biparental/ 根目录
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(igraph); library(data.table)})
outdir <- "biparental"
t0 <- Sys.time()
pops <- setdiff(c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
                  "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK"),
                c("BBC","BLB"))
cat("pops:", length(pops), "\n")

## ---------- 整合全部上位效应对 ----------
lst <- list()
for (p in pops) {
  for (i in 1:7) {
    d <- fread(paste0("Epistasis_", p, "_P", i, "_txt"),
               select = c("snp1","snp2"))
    d$pop <- p; d$trait <- i
    lst[[length(lst)+1]] <- d
  }
}
all <- rbindlist(lst)
cat("total rows:", nrow(all), "\n")
all$snp1 <- as.character(all$snp1); all$snp2 <- as.character(all$snp2)
all$key <- paste(pmin(all$snp1, all$snp2), pmax(all$snp1, all$snp2), sep = "|")

## 去重 + 统计每条边被多少群体/群体×性状组合检测到
uk <- unique(all$key)
e1 <- sub("\\|.*", "", uk); e2 <- sub(".*\\|", "", uk)
n_pops   <- as.integer(tapply(all$pop, all$key, function(x) length(unique(x)))[uk])
n_combos <- as.integer(table(all$key)[uk])
edges <- data.frame(snp1 = e1, snp2 = e2, n_pops = n_pops, n_combos = n_combos)
cat("unique edges:", nrow(edges), "\n")

## ---------- 构建统一网络 ----------
g <- simplify(graph_from_data_frame(edges[, c("snp1","snp2")], directed = FALSE))
cat("network: nodes =", vcount(g), " edges =", ecount(g), "\n")

## ---------- 模块检测（fast greedy, 每连通分量, >=5节点） ----------
MIN_MODULE <- 5
comp <- components(g)
memb <- setNames(rep(NA_character_, vcount(g)), V(g)$name)
mid <- 0
for (cid in seq_len(comp$no)) {
  vv <- V(g)$name[comp$membership == cid]
  gs <- induced_subgraph(g, vv)
  if (vcount(gs) < MIN_MODULE) next
  cm <- tryCatch(cluster_fast_greedy(gs), error = function(e) NULL)
  if (is.null(cm)) next
  for (m in communities(cm)) {
    if (length(m) >= MIN_MODULE) { mid <- mid + 1; memb[m] <- sprintf("BP_M%04d", mid) }
  }
}
node_tab <- data.frame(HB = names(memb), module = memb, stringsAsFactors = FALSE)
node_tab$module[is.na(node_tab$module)] <- "non_module(<5)"
deg <- table(c(edges$snp1, edges$snp2))
node_tab$edge_count <- as.integer(deg[node_tab$HB])
node_tab$edge_count[is.na(node_tab$edge_count)] <- 0L
write.table(node_tab, file.path(outdir, "biparental_node_module_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fm  <- match(edges$snp1, node_tab$HB)
fm2 <- match(edges$snp2, node_tab$HB)
edge_out <- data.frame(snp1 = edges$snp1, snp2 = edges$snp2,
                       n_pops = edges$n_pops, n_combos = edges$n_combos,
                       module_snp1 = node_tab$module[fm],
                       module_snp2 = node_tab$module[fm2])
edge_out$intra_module <- ifelse(edge_out$module_snp1 == edge_out$module_snp2,
                                edge_out$module_snp1, "cross_module")
edge_out <- edge_out[order(-edge_out$n_pops), ]
write.table(edge_out, file.path(outdir, "biparental_edge_module_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("modules:", mid, " nodes in modules:", sum(node_tab$module != "non_module(<5)"), "\n")

## ---------- Gephi 表（按模块着色） ----------
mod_sizes <- table(node_tab$module[node_tab$module != "non_module(<5)"])
nmod <- length(mod_sizes)
cc <- grDevices::col2rgb(grDevices::hsv(h = seq(0, 1, length.out = nmod + 1)[-1], s = 0.85, v = 0.95))
col_df <- data.frame(module = names(mod_sizes), r = cc[1, ], g = cc[2, ], b = cc[3, ])
col_all <- rbind(col_df, data.frame(module = "non_module(<5)", r = 160, g = 160, b = 160))
node_tab$module_size <- as.integer(mod_sizes[node_tab$module])
node_tab$module_size[is.na(node_tab$module_size)] <- 0L
node_tab <- cbind(node_tab, col_all[match(node_tab$module, col_all$module), c("r","g","b")])
node_tab$Label <- node_tab$HB
nodes_out <- node_tab[, c("HB","Label","module","module_size","edge_count","r","g","b")]
names(nodes_out)[1] <- "Id"
write.csv(nodes_out, file.path(outdir, "gephi_nodes_biparental_csv"), row.names = FALSE)
edges_out <- data.frame(Source = edge_out$snp1, Target = edge_out$snp2,
                        Type = "undirected", Weight = edge_out$n_pops,
                        module = edge_out$intra_module)
edges_out <- cbind(edges_out, col_all[match(edges_out$module, col_all$module), c("r","g","b")])
write.csv(edges_out, file.path(outdir, "gephi_edges_biparental_csv"), row.names = FALSE)

sum_row <- data.frame(dataset = "biparental_23pops", n_pops_scanned = length(pops),
                      n_nodes = vcount(g), n_edges = ecount(g), n_modules = nmod,
                      largest_module = max(mod_sizes), median_module = round(median(mod_sizes), 1),
                      nodes_in_modules = sum(node_tab$module != "non_module(<5)"))
write.table(sum_row, file.path(outdir, "module_summary_biparental_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== DONE biparental merged network, elapsed:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ===\n")
