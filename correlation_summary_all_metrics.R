# correlation_summary_all_metrics.R
# 四类指标与表型的相关性统一汇总表：
#   1) 优异单倍型数目   (gwas_haplotype_counts_P*_txt: n_favorable)
#   2) 优异上位效应数目 (gwas_epistasis_netstats_P*_txt: edge_count, favorable)
#   3) 模块数目         (gwas_module_stats_P*_txt: n_modules_covered)
#   4) 模块平均完整度   (gwas_module_stats_P*_txt: mean_completeness)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(dplyr); library(tidyr)})

tr7 <- paste0("P", 1:7)
phe <- read.delim("phe2.txt", header = FALSE, comment.char = "#")
phe$V1 <- as.character(phe$V1)
names(phe) <- c("ID", tr7)

cor_one <- function(x, y) {
  ct <- cor.test(x, y, method = "pearson")
  st <- cor.test(x, y, method = "spearman", exact = FALSE)
  data.frame(n = length(x),
             pearson_r = unname(ct$estimate), pearson_p = ct$p.value,
             spearman_rho = unname(st$estimate), spearman_p = st$p.value)
}

rows <- list()
for (i in 1:7) {
  t <- tr7[i]
  y <- phe[[t]]
  names(y) <- phe$ID          # 给表型向量补名字，便于按ID取子集

  ## 1) 优异单倍型数目
  hc <- read.delim(paste0("gwas_haplotype_counts_P", i, "_txt"))
  m <- match(hc$ID, names(y))
  keep <- !is.na(m) & !is.na(y[m])
  rows[[length(rows) + 1]] <- cbind(trait = t, metric = "n_favorable_haplotypes",
                                    cor_one(hc$n_favorable[keep], y[m[keep]]))

  ## 2) 优异上位效应数目（有利网络 edge_count）
  st <- read.delim(paste0("gwas_epistasis_netstats_P", i, "_txt"), row.names = 1)
  m2 <- match(rownames(st), names(y))
  keep2 <- !is.na(m2) & !is.na(y[m2])
  rows[[length(rows) + 1]] <- cbind(trait = t, metric = "n_favorable_epistasis_edges",
                                    cor_one(st$edge_count[keep2], y[m2[keep2]]))

  ## 3) 模块数目  4) 模块平均完整度
  ms <- read.delim(paste0("gwas_module_stats_P", i, "_txt"))
  m3 <- match(ms$ID, names(y))
  keep3 <- !is.na(m3) & !is.na(y[m3])
  rows[[length(rows) + 1]] <- cbind(trait = t, metric = "n_modules_covered",
                                    cor_one(ms$n_modules_covered[keep3], y[m3[keep3]]))
  rows[[length(rows) + 1]] <- cbind(trait = t, metric = "mean_module_completeness",
                                    cor_one(ms$mean_completeness[keep3], y[m3[keep3]]))
}
res <- do.call(rbind, rows)
res$metric <- factor(res$metric, levels = c("n_favorable_haplotypes",
                                            "n_favorable_epistasis_edges",
                                            "n_modules_covered",
                                            "mean_module_completeness"))
res <- as.data.frame(res)
res <- res[order(res$trait, res$metric), ]

## 长表 + 宽表两种输出
write.table(res, "phenotype_correlation_all_metrics_long_txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
wide <- res %>%
  transmute(trait, metric,
            pearson_r = round(pearson_r, 3), pearson_p = signif(pearson_p, 2),
            spearman_rho = round(spearman_rho, 3), spearman_p = signif(spearman_p, 2), n) %>%
  pivot_wider(names_from = metric,
              values_from = c(pearson_r, pearson_p, spearman_rho, spearman_p, n),
              names_glue = "{metric}_{.value}") %>%
  as.data.frame()
wide <- wide[order(wide$trait), ]
write.table(wide, "phenotype_correlation_all_metrics_wide_txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("===== 相关性统一汇总（长表） =====\n")
print(res %>% mutate(across(where(is.numeric), ~round(.x, 4))), row.names = FALSE, digits = 3)
cat("\nDONE\n")
