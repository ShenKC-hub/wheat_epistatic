setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages(library(igraph))
MIN_MODULE <- 5

for (i in 1:7) {
  t0 <- Sys.time()
  e <- read.delim(paste0("gwas_epistasis_edge_class_P", i, "_txt"))
  fav <- e[e$class == "favorable", c("snp1","snp2")]
  if (nrow(fav) == 0) next
  g <- simplify(graph_from_data_frame(fav, directed=FALSE))

  ## 逐连通分量做 fast greedy，全局统一编号模块
  comp <- components(g)
  memb <- setNames(rep(NA_character_, vcount(g)), V(g)$name)
  mid <- 0
  for (cid in seq_len(comp$no)) {
    vv <- V(g)$name[comp$membership == cid]
    gs <- induced_subgraph(g, vv)
    if (vcount(gs) < MIN_MODULE) next
    cm <- tryCatch(cluster_fast_greedy(gs), error=function(e) NULL)
    if (is.null(cm)) next
    for (m in communities(cm)) {
      if (length(m) >= MIN_MODULE) {
        mid <- mid + 1
        memb[m] <- sprintf("P%d_M%04d", i, mid)
      }
    }
  }
  ## 节点度：参与的有利边数（去重）
  fkey <- apply(fav, 1, function(x) paste(sort(x), collapse="|"))
  fav_u <- fav[!duplicated(fkey), ]
  deg <- table(c(fav_u$snp1, fav_u$snp2))
  node_tab <- data.frame(HB = names(memb), module = memb, stringsAsFactors=FALSE)
  node_tab$edge_count <- as.integer(deg[node_tab$HB])
  node_tab$edge_count[is.na(node_tab$edge_count)] <- 0L
  node_tab$module[is.na(node_tab$module)] <- "non_module(<5)"
  write.table(node_tab, paste0("gwas_epistasis_node_module_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  ## 边 -> 模块对应表（favorable 边；两端同模块记入该模块，跨模块标注 cross_module）
  fm <- match(fav$snp1, node_tab$HB); fm2 <- match(fav$snp2, node_tab$HB)
  edge_tab <- data.frame(snp1 = fav$snp1, snp2 = fav$snp2,
                         module_snp1 = node_tab$module[fm],
                         module_snp2 = node_tab$module[fm2],
                         stringsAsFactors=FALSE)
  edge_tab$intra_module <- ifelse(edge_tab$module_snp1 == edge_tab$module_snp2,
                                  edge_tab$module_snp1, "cross_module")
  write.table(edge_tab, paste0("gwas_epistasis_edge_module_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("P", i, ": nodes", vcount(g), "edges", ecount(g),
      "| modules", mid, "| intra-module edges",
      sum(edge_tab$intra_module != "cross_module"),
      "| elapsed", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "s\n")
}
cat("ALL DONE\n")
