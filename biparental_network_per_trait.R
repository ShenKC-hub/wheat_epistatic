# biparental_network_per_trait.R
# 23个双亲群体（排除BBC、BLB）按性状分别整合上位效应对（不区分有利/不利），
# 每个性状构建一个网络（共7个）：模块对应表 + Gephi 表
# 结果保存到 biparental/ 根目录，文件名含性状编号
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(igraph); library(data.table)})
outdir <- "biparental"
t0 <- Sys.time()
pops <- setdiff(c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
                  "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK"),
                c("BBC","BLB"))
MIN_MODULE <- 5
summary_all <- NULL

for (i in 1:7) {
  ti <- Sys.time()
  lst <- list()
  for (p in pops) {
    d <- fread(paste0("Epistasis_", p, "_P", i, "_txt"), select = c("snp1","snp2"))
    d$pop <- p
    lst[[length(lst)+1]] <- d
  }
  all <- rbindlist(lst)
  all$snp1 <- as.character(all$snp1); all$snp2 <- as.character(all$snp2)
  all$key <- paste(pmin(all$snp1, all$snp2), pmax(all$snp1, all$snp2), sep = "|")

  uk <- unique(all$key)
  e1 <- sub("\\|.*", "", uk); e2 <- sub(".*\\|", "", uk)
  n_pops <- as.integer(tapply(all$pop, all$key, function(x) length(unique(x)))[uk])
  edges <- data.frame(snp1 = e1, snp2 = e2, n_pops = n_pops)
  cat("P", i, ": rows", nrow(all), " unique edges", nrow(edges), "\n")

  g <- simplify(graph_from_data_frame(edges[, c("snp1","snp2")], directed = FALSE))

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
      if (length(m) >= MIN_MODULE) { mid <- mid + 1; memb[m] <- sprintf("BP_P%d_M%04d", i, mid) }
    }
  }
  node_tab <- data.frame(HB = names(memb), module = memb, stringsAsFactors = FALSE)
  node_tab$module[is.na(node_tab$module)] <- "non_module(<5)"
  deg <- table(c(edges$snp1, edges$snp2))
  node_tab$edge_count <- as.integer(deg[node_tab$HB])
  node_tab$edge_count[is.na(node_tab$edge_count)] <- 0L
  write.table(node_tab, file.path(outdir, paste0("biparental_node_module_P", i, "_txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  fm  <- match(edges$snp1, node_tab$HB)
  fm2 <- match(edges$snp2, node_tab$HB)
  edge_out <- data.frame(snp1 = edges$snp1, snp2 = edges$snp2, n_pops = edges$n_pops,
                         module_snp1 = node_tab$module[fm],
                         module_snp2 = node_tab$module[fm2])
  edge_out$intra_module <- ifelse(edge_out$module_snp1 == edge_out$module_snp2,
                                  edge_out$module_snp1, "cross_module")
  edge_out <- edge_out[order(-edge_out$n_pops), ]
  write.table(edge_out, file.path(outdir, paste0("biparental_edge_module_P", i, "_txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("P", i, ": nodes", vcount(g), " edges", ecount(g), " modules", mid,
      " nodes_in_modules", sum(node_tab$module != "non_module(<5))"), "\n")

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
  write.csv(nodes_out, file.path(outdir, paste0("gephi_nodes_biparental_P", i, "_csv")),
            row.names = FALSE)
  edges_out <- data.frame(Source = edge_out$snp1, Target = edge_out$snp2,
                          Type = "undirected", Weight = edge_out$n_pops,
                          module = edge_out$intra_module)
  edges_out <- cbind(edges_out, col_all[match(edges_out$module, col_all$module), c("r","g","b")])
  write.csv(edges_out, file.path(outdir, paste0("gephi_edges_biparental_P", i, "_csv")),
            row.names = FALSE)

  summary_all <- rbind(summary_all,
    data.frame(trait = paste0("P", i), n_pops_scanned = length(pops),
               n_rows_raw = nrow(all), n_nodes = vcount(g), n_edges = ecount(g),
               n_modules = nmod,
               largest_module = if (nmod > 0) max(mod_sizes) else 0,
               median_module = if (nmod > 0) round(median(mod_sizes), 1) else 0,
               nodes_in_modules = sum(node_tab$module != "non_module(<5)"),
               edges_2pops = sum(edge_out$n_pops >= 2),
               elapsed_min = round(as.numeric(difftime(Sys.time(), ti, units = "mins")), 1)))
}
write.table(summary_all, file.path(outdir, "module_summary_biparental_per_trait_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_all)
cat("=== DONE per-trait biparental networks, total elapsed:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ===\n")
