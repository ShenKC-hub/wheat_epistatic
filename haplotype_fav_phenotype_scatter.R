# haplotype_fav_phenotype_scatter.R
# 每个性状：材料携带的有利单倍型数目 vs 表型值 散点图 + 拟合光滑曲线(loess)与线性回归
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr)})

phe <- read.delim("phe2.txt", header = FALSE, comment.char = "#")
phe$ID <- as.character(phe$V1)
phe$V1 <- NULL
names(phe) <- c(paste0("P", 1:7), "ID")
phe_l <- reshape(phe, direction = "long", varying = paste0("P", 1:7),
                 v.names = "phenotype", times = paste0("P", 1:7),
                 timevar = "trait", idvar = "ID")
phe_l$trait <- factor(phe_l$trait, levels = paste0("P", 1:7))

dat <- list()
cor_tab <- NULL
for (i in 1:7) {
  c1 <- read.delim(paste0("gwas_haplotype_counts_P", i, "_txt"))
  c1$trait <- paste0("P", i)
  dat[[i]] <- c1
  d <- merge(c1, phe_l, by = c("ID", "trait"))
  d <- d[!is.na(d$phenotype), ]
  ct <- cor.test(d$n_favorable, d$phenotype, method = "pearson")
  st <- cor.test(d$n_favorable, d$phenotype, method = "spearman", exact = FALSE)
  cor_tab <- rbind(cor_tab, data.frame(
    trait = paste0("P", i), n = nrow(d),
    pearson_r = unname(ct$estimate), pearson_p = ct$p.value,
    spearman_rho = unname(st$estimate), spearman_p = st$p.value))
  dat[[i]]$phenotype <- NA  # 防止重复行; 实际绘图用合并后的 d
}
cor_tab$label <- sprintf("r = %.2f, P = %.1e\nrho = %.2f, P = %.1e",
                         cor_tab$pearson_r, cor_tab$pearson_p,
                         cor_tab$spearman_rho, cor_tab$spearman_p)
write.table(cor_tab, "gwas_haplotype_fav_phenotype_cor_txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## 合并所有性状的绘图数据
plot_dat <- NULL
for (i in 1:7) {
  c1 <- read.delim(paste0("gwas_haplotype_counts_P", i, "_txt"))
  c1$trait <- paste0("P", i)
  d <- merge(c1, phe_l, by = c("ID", "trait"))
  d <- d[!is.na(d$phenotype), c("ID", "trait", "phenotype", "n_favorable")]
  plot_dat <- rbind(plot_dat, d)
}
plot_dat$trait <- factor(plot_dat$trait, levels = paste0("P", 1:7))
lab_df <- cor_tab[, c("trait", "label")]
names(lab_df)[1] <- "trait"

p <- ggplot(plot_dat, aes(x = n_favorable, y = phenotype)) +
  geom_point(shape = 21, fill = "#2C64AE", colour = "grey25",
             alpha = 0.55, size = 1.6, stroke = 0.3) +
  geom_smooth(method = "loess", se = TRUE, colour = "#DF0101",
              fill = "#DF0101", alpha = 0.12, linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey20",
              linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~trait, nrow = 2, scales = "free") +
  geom_text(data = lab_df, aes(label = label), x = Inf, y = Inf,
            hjust = 1.05, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  labs(title = "Favorable haplotype count vs phenotype for each trait",
       x = "Number of favorable haplotypes carried",
       y = "Phenotype value") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

pdf("gwas_haplotype_fav_phenotype_scatter.pdf", width = 12, height = 7); print(p); dev.off()
png("gwas_haplotype_fav_phenotype_scatter.png", width = 12, height = 7,
    units = "in", res = 150); print(p); dev.off()
cat("correlation summary:\n")
print(cor_tab, row.names = FALSE, digits = 3)
cat("DONE\n")
