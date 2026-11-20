setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
MIN_MODULE <- 5
summary_tab <- NULL

for (i in 1:7) {
  nt <- read.delim(paste0("gwas_epistasis_node_module_P", i, "_txt"))
  e  <- read.delim(paste0("gwas_epistasis_edge_class_P", i, "_txt"))
  fav <- e[e$class == "favorable", c("snp1","snp2")]
  key <- apply(fav, 1, function(x) paste(sort(x), collapse="|"))
  fav <- fav[!duplicated(key), ]

  ## 模块大小与颜色（每性状独立取色环，non_module 固定灰色）
  mod_sizes <- table(nt$module[nt$module != "non_module(<5)"])
  nmod <- length(mod_sizes)
  cols <- grDevices::hsv(h = seq(0, 1, length.out = nmod + 1)[-1], s = 0.85, v = 0.95)
  cols <- grDevices::col2rgb(cols)
  col_df <- data.frame(module = names(mod_sizes),
                       r = cols[1, ], g = cols[2, ], b = cols[3, ],
                       stringsAsFactors = FALSE)
  grey <- data.frame(module = c("non_module(<5)", "cross_module"),
                     r = c(160, 120), g = c(160, 120), b = c(160, 120))
  col_all <- rbind(col_df, grey)
  nt <- cbind(nt, col_all[match(nt$module, col_all$module), c("r","g","b")])
  nt$module_size <- as.integer(mod_sizes[nt$module])
  nt$module_size[is.na(nt$module_size)] <- 0L
  nt$Label <- nt$HB
  nodes_out <- nt[, c("Id" = "HB", "Label", "module", "module_size", "edge_count", "r","g","b")]
  names(nodes_out)[1] <- "Id"
  write.csv(nodes_out, paste0("gephi_nodes_P", i, "_csv"), row.names = FALSE, quote = TRUE)

  ## 边表（去重后的有利边），按模块着色
  em <- match(fav$snp1, nt$HB); em2 <- match(fav$snp2, nt$HB)
  edgemod <- ifelse(nt$module[em] == nt$module[em2], nt$module[em], "cross_module")
  edges_out <- data.frame(Source = fav$snp1, Target = fav$snp2, Type = "undirected",
                          Weight = 1, module = edgemod)
  edges_out <- cbind(edges_out, col_all[match(edgemod, col_all$module), c("r","g","b")])
  write.csv(edges_out, paste0("gephi_edges_P", i, "_csv"), row.names = FALSE, quote = TRUE)

  summary_tab <- rbind(summary_tab,
    data.frame(trait = paste0("P", i),
               n_nodes = nrow(nt),
               n_edges_favorable = nrow(fav),
               n_modules = nmod,
               largest_module = max(mod_sizes),
               median_module = round(median(mod_sizes), 1),
               nodes_in_modules = sum(nt$module != "non_module(<5)")))
}
write.csv(summary_tab, "epistasis_module_summary_csv", row.names = FALSE)
print(summary_tab)
cat("ALL DONE\n")
