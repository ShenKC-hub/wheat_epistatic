# h2_epi_heatmap.R
# 25个双亲群体 x 7个性状 的上位效应遗传力(h2_epi)热图，标注两位有效数字
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages(library(pheatmap))

pops <- c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
          "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK")
traits <- paste0("P", 1:7)

mat <- matrix(NA_real_, nrow = length(pops), ncol = length(traits),
              dimnames = list(pops, traits))
for (p in pops) for (i in 1:7) {
  f <- paste0("Epistasis_h2_", p, "_P", i, "_txt")
  if (file.exists(f)) {
    d <- read.delim(f)
    mat[p, i] <- d$h2_epi[1]
  }
}

## 数值矩阵（保留两位有效数字）与标签矩阵
mat_sig <- signif(mat, 2)
label_mat <- matrix(ifelse(is.na(mat), "", sprintf("%.2g", mat)),
                    nrow = nrow(mat), dimnames = dimnames(mat))
mat_t <- t(mat)
label_t <- t(label_mat)
write.csv(mat_sig, "epistasis_h2_epi_matrix_25pops_csv", row.names = TRUE)

## 热图（横坐标=群体，纵坐标=性状，白至#2C64AE；NA留白）
pheatmap(mat_t,
         display_numbers = label_t,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#FFFFFF", "#2C64AE"))(100),
         na_col = "grey90",
         border_color = "grey80",
         angle_col = 45,
         main = "Epistatic heritability (h2_epi) in 25 biparental populations",
         fontsize = 10, fontsize_number = 9,
         legend = TRUE, filename = "h2_epi_heatmap_25pops_T.pdf")
pheatmap(mat_t,
         display_numbers = label_t,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#FFFFFF", "#2C64AE"))(100),
         na_col = "grey90",
         border_color = "grey80",
         angle_col = 45,
         main = "Epistatic heritability (h2_epi) in 25 biparental populations",
         fontsize = 10, fontsize_number = 9,
         legend = TRUE, filename = "h2_epi_heatmap_25pops_T.png", width = 9, height = 4)

cat("h2_epi 概览:\n")
print(round(mat_sig, 3))
cat("NA cells:", sum(is.na(mat)), "\n")
cat("DONE\n")
