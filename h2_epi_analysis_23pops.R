# h2_epi_analysis_23pops.R
# 双亲群体(排除BBC/BLB)上位效应遗传力：热图 + 分布图 + 显著上位对数与遗传力的相关分析
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(pheatmap); library(ggplot2); library(data.table); library(dplyr)})

pops <- c("AVS","BBR","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
          "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK")
traits <- paste0("P", 1:7)
stopifnot(length(pops) == 23)

## ---------- 数据收集 ----------
mat   <- matrix(NA_real_, 23, 7, dimnames = list(pops, traits))
npair <- matrix(NA_real_, 23, 7, dimnames = list(pops, traits))
for (p in pops) for (i in 1:7) {
  fh <- paste0("Epistasis_h2_", p, "_P", i, "_txt")
  fp <- paste0("Epistasis_", p, "_P", i, "_txt")
  if (file.exists(fh)) mat[p, i] <- read.delim(fh)$h2_epi[1]
  if (file.exists(fp)) {
    d <- fread(fp, select = c("snp1","snp2"))
    k <- paste(pmin(as.character(d$snp1), as.character(d$snp2)),
               pmax(as.character(d$snp1), as.character(d$snp2)), sep = "|")
    npair[p, i] <- length(unique(k))          # 去重后的上位效应对数(全部类别)
  }
}

dat <- expand.grid(trait = traits, pop = pops)
dat$h2_epi <- as.vector(t(mat))
dat$n_pairs <- as.vector(t(npair))
write.csv(dat, "biparental_23pops_h2_epi_npairs_csv", row.names = FALSE)

## ---------- 1) 热图（横=群体，纵=性状，白至#2C64AE，两位有效数字） ----------
label_mat <- matrix(ifelse(is.na(mat), "", sprintf("%.2g", mat)),
                    nrow = nrow(mat), dimnames = dimnames(mat))
pdf("h2_epi_heatmap_23pops.pdf", width = 11, height = 4.5)
pheatmap(t(mat), display_numbers = t(label_mat),
         cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45,
         color = colorRampPalette(c("#FFFFFF", "#2C64AE"))(100),
         na_col = "grey90", border_color = "grey80",
         main = "Epistatic heritability (h2_epi) in 23 biparental populations",
         fontsize = 10, fontsize_number = 9)
dev.off()
png("h2_epi_heatmap_23pops.png", width = 11, height = 4.5, res = 150, units = "in")
pheatmap(t(mat), display_numbers = t(label_mat),
         cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45,
         color = colorRampPalette(c("#FFFFFF", "#2C64AE"))(100),
         na_col = "grey90", border_color = "grey80",
         main = "Epistatic heritability (h2_epi) in 23 biparental populations",
         fontsize = 10, fontsize_number = 9)
dev.off()

## ---------- 2) 遗传力分布图：总体直方图+密度 & 分性状箱线图 ----------
dd <- dat[!is.na(dat$h2_epi), ]
p_hist <- ggplot(dd, aes(x = h2_epi)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05,
                 fill = "#2C64AE", colour = "grey30", alpha = 0.75) +
  geom_density(colour = "#DF0101", linewidth = 0.9) +
  geom_vline(xintercept = median(dd$h2_epi), linetype = "dashed", colour = "grey25") +
  annotate("text", x = median(dd$h2_epi), y = Inf, hjust = -0.05, vjust = 2,
           label = sprintf("median = %.2f", median(dd$h2_epi)), size = 3.4) +
  labs(title = "Distribution of epistatic heritability across 23 biparental populations",
       x = "h2_epi", y = "Density") +
  theme_bw() + theme(panel.grid.minor = element_blank())

p_box <- ggplot(dd, aes(x = reorder(trait, h2_epi, median), y = h2_epi, fill = trait)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.5, shape = 21, fill = "grey90") +
  labs(title = "h2_epi by trait", x = "Trait", y = "h2_epi", fill = "Trait") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())

pdf("h2_epi_distribution_hist_23pops.pdf", width = 7, height = 5); print(p_hist); dev.off()
pdf("h2_epi_distribution_box_23pops.pdf", width = 7, height = 5); print(p_box); dev.off()
png("h2_epi_distribution_hist_23pops.png", width = 7, height = 5, res = 150, units = "in")
print(p_hist); dev.off()
png("h2_epi_distribution_box_23pops.png", width = 7, height = 5, res = 150, units = "in")
print(p_box); dev.off()

## ---------- 3) 显著上位对数 vs 遗传力：散点 + 线性相关分析 ----------
dd2 <- dd[dd$n_pairs > 0, ]
dd2$log10_pairs <- log10(dd2$n_pairs)
ct <- cor.test(dd2$log10_pairs, dd2$h2_epi, method = "pearson")
lmfit <- lm(h2_epi ~ log10_pairs, dd2)
r2 <- summary(lmfit)$r.squared
pv <- pf(summary(lmfit)$fstatistic[1], 1, summary(lmfit)$fstatistic[2], lower.tail = FALSE)
cat(sprintf("线性相关(log10对数 vs h2_epi): Pearson r=%.3f, P=%.2e, R^2=%.3f, n=%d\n",
            unname(ct$estimate), ct$p.value, r2, nrow(dd2)))
ct_raw <- cor.test(dd2$n_pairs, dd2$h2_epi, method = "pearson")
cat(sprintf("原始尺度: r=%.3f, P=%.2e\n", unname(ct_raw$estimate), ct_raw$p.value))

p_sc <- ggplot(dd2, aes(x = n_pairs, y = h2_epi)) +
  geom_point(shape = 21, fill = "#2C64AE", colour = "grey25",
             size = 2.4, alpha = 0.75, stroke = 0.4) +
  geom_smooth(method = "lm", colour = "#DF0101", fill = "#DF0101",
              alpha = 0.12, linewidth = 0.8) +
  scale_x_log10() +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.6,
           label = sprintf("log10 scale: r = %.2f, P = %.1e, R\u00b2 = %.2f, n = %d",
                           unname(ct$estimate), ct$p.value, r2, nrow(dd2))) +
  labs(title = "Detected epistasis pairs vs epistatic heritability (all classes pooled)",
       x = "Number of significant HB-HB pairs (log10 scale)",
       y = "h2_epi") +
  theme_bw() + theme(panel.grid.minor = element_blank())

pdf("h2_epi_vs_npairs_scatter.pdf", width = 8, height = 6); print(p_sc); dev.off()
png("h2_epi_vs_npairs_scatter.png", width = 8, height = 6, res = 150, units = "in")
print(p_sc); dev.off()
cat("DONE\n")