# trait_cor_heatmap.R
# 681份自然群体 P1-P7 农艺性状(第2-8列)的 Pearson 相关热图
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages(library(pheatmap))

phe <- read.delim("phe2.txt", header = FALSE, comment.char = "#")
X <- phe[, 2:8]
colnames(X) <- paste0("P", 1:7)

r <- cor(X, use = "pairwise.complete.obs")
write.csv(signif(r, 3), "trait_correlation_matrix_csv", row.names = TRUE)

## 显著性（n=646 完整观测对）
label <- matrix(sprintf("%.2f", r), nrow = nrow(r))
for (ii in seq_len(nrow(r))) for (jj in seq_len(ncol(r))) {
  if (ii == jj) next
  nobs <- sum(complete.cases(X[, ii], X[, jj]))
  tstat <- r[ii, jj] * sqrt((nobs - 2) / (1 - r[ii, jj]^2))
  pv <- 2 * pt(-abs(tstat), nobs - 2)
  if (pv < 0.001) label[ii, jj] <- paste0(label[ii, jj], "***")
  else if (pv < 0.01) label[ii, jj] <- paste0(label[ii, jj], "**")
  else if (pv < 0.05) label[ii, jj] <- paste0(label[ii, jj], "*")
}
diag(label) <- "1"

pheatmap(r,
         display_numbers = label,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#B2182B", "#FFFFFF", "#2C64AE"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color = "grey80",
         main = "Phenotypic correlations among P1-P7 (n = 681)",
         fontsize = 11, fontsize_number = 10,
         legend = TRUE, filename = "trait_cor_heatmap.pdf")
pheatmap(r,
         display_numbers = label,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#B2182B", "#FFFFFF", "#2C64AE"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color = "grey80",
         main = "Phenotypic correlations among P1-P7 (n = 681)",
         fontsize = 11, fontsize_number = 10,
         legend = TRUE, filename = "trait_cor_heatmap.png",
         width = 7, height = 6)
cat("cor range:", round(min(r), 3), "to", round(max(r), 3), "\nDONE\n")
