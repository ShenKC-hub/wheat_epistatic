# haplotype_class_and_counts.R
# 1) 每个性状：按 p_snp/alpha_snp 将单倍型分为 有利/不利/无显著加性效应
# 2) 每个材料携带的有利/不利单倍型数目（genotype==2 计携带）
# 3) 按育种进程分组统计数目变化（参考 epistasis_netstats_by_group.R）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr)})

P_CUT <- 0.05

## ---------- 数据读取（与 GWAS_epistasis_effect_h2_P*.R 一致） ----------
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")
pca <- read.delim("pca.txt", header=FALSE, comment.char="#", sep="\t")
phe <- merge(phe, pca, by="V1")
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
phe <- merge(phe, genotype, by.x="V1", by.y="ID")
snpA <- na.omit(phe[, c(1, 12:ncol(phe))])   # 跳过7性状+3PC，仅ID与单倍型
ids <- as.character(snpA$V1)
G <- as.matrix(sapply(snpA[, -1], as.numeric))
cat("samples:", nrow(G), " SNPs:", ncol(G), "\n")

grp <- read.delim("sample_group.txt", header=FALSE,
                  col.names=c("ID", "group"), comment.char="#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace", "Introduced", "Pre-1960",
               "1961-1980", "1981-2000", "2001-2020")

counts_all <- list()
class_summ <- NULL

for (i in 1:7) {
  eff <- read.delim(paste0("gwas_Epistasis_effect_P", i, "_txt"))
  eff$snp1 <- as.character(eff$snp1); eff$snp2 <- as.character(eff$snp2)

  ## ---- 1) 单倍型分类：取该单倍型最显著出现的 (p, alpha) 判定 ----
  occ <- rbind(
    data.frame(HB = eff$snp1, p = eff$p_snp1, alpha = eff$alpha_snp1),
    data.frame(HB = eff$snp2, p = eff$p_snp2, alpha = eff$alpha_snp2))
  agg <- occ %>% group_by(HB) %>%
    summarise(n_occurrences = n(),
              n_sig = sum(p < P_CUT),
              min_p = min(p),
              alpha_at_min_p = alpha[which.min(p)],
              .groups = "drop") %>% as.data.frame()
  agg$class <- with(agg, ifelse(min_p < P_CUT & alpha_at_min_p > 0, "favorable",
                         ifelse(min_p < P_CUT & alpha_at_min_p < 0, "unfavorable",
                                "no_significant_additive")))
  write.table(agg, paste0("gwas_haplotype_class_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  class_summ <- rbind(class_summ,
    data.frame(trait = paste0("P", i), table(agg$class)))
  cat("P", i, ":", paste(names(table(agg$class)), table(agg$class), collapse = "  "), "\n")

  ## ---- 2) 每个材料携带的有利/不利单倍型数目 ----
  inG <- agg$HB %in% colnames(G)
  if (any(!inG)) cat("  WARNING:", sum(!inG), "haplotypes not in genotype matrix\n")
  fav <- agg$HB[agg$class == "favorable" & inG]
  unf <- agg$HB[agg$class == "unfavorable" & inG]
  n_fav <- as.integer(rowSums(G[, fav, drop = FALSE] == 2))
  n_unf <- as.integer(rowSums(G[, unf, drop = FALSE] == 2))
  cnt <- data.frame(ID = ids, trait = paste0("P", i),
                    n_favorable = n_fav, n_unfavorable = n_unf,
                    n_fav_total = length(fav), n_unfav_total = length(unf))
  write.table(cnt, paste0("gwas_haplotype_counts_P", i, "_txt"),
              sep="\t", quote=FALSE, row.names=FALSE)
  counts_all[[i]] <- cnt
}
counts_all <- do.call(rbind, counts_all)
counts_all <- merge(counts_all, grp, by = "ID", all.x = TRUE)
counts_all$trait <- factor(counts_all$trait, levels = paste0("P", 1:7))
counts_all$group <- factor(counts_all$group, levels = grp_order)
write.csv(class_summ, "gwas_haplotype_class_summary_csv", row.names = FALSE)

## ---- 3) 育种进程分组统计（均值 ± 标准差 + 柱状图，同参考代码风格） ----
summ <- counts_all %>%
  filter(!is.na(group)) %>%
  group_by(trait, group) %>%
  summarise(n = n(),
            fav_mean = mean(n_favorable), fav_sd = sd(n_favorable),
            unf_mean = mean(n_unfavorable), unf_sd = sd(n_unfavorable),
            .groups = "drop") %>% as.data.frame()
summ$fav_ms <- sprintf("%.2f±%.2f", summ$fav_mean, summ$fav_sd)
summ$unf_ms <- sprintf("%.2f±%.2f", summ$unf_mean, summ$unf_sd)
out_wide <- summ %>%
  transmute(trait, group, n, fav_mean, fav_sd, fav_mean_sd = fav_ms,
            unf_mean, unf_sd, unf_mean_sd = unf_ms) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
write.table(out_wide, "gwas_haplotype_counts_by_group.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("saved: gwas_haplotype_counts_by_group.txt\n")
print(out_wide, right = FALSE)

pd <- position_dodge(0.8)
p1 <- ggplot(summ, aes(x = group, y = fav_mean, fill = trait)) +
  geom_bar(stat = "identity", width = 0.75,
           colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = fav_mean - fav_sd, ymax = fav_mean + fav_sd),
                width = 0.25, linewidth = 0.3) +
  facet_wrap(~trait, nrow = 2, scales = "free_y") +
  labs(title = "Favorable haplotypes carried per accession by breeding era",
       x = "Breeding era", y = "Mean n_favorable (±SD)", fill = "Trait") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     strip.text = element_text(face = "bold"))
p2 <- ggplot(summ, aes(x = group, y = unf_mean, fill = trait)) +
  geom_bar(stat = "identity", width = 0.75,
           colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = unf_mean - unf_sd, ymax = unf_mean + unf_sd),
                width = 0.25, linewidth = 0.3) +
  facet_wrap(~trait, nrow = 2, scales = "free_y") +
  labs(title = "Unfavorable haplotypes carried per accession by breeding era",
       x = "Breeding era", y = "Mean n_unfavorable (±SD)", fill = "Trait") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     strip.text = element_text(face = "bold"))
pdf("gwas_haplotype_favorable_by_group.pdf", width = 10, height = 6); print(p1); dev.off()
pdf("gwas_haplotype_unfavorable_by_group.pdf", width = 10, height = 6); print(p2); dev.off()
png("gwas_haplotype_favorable_by_group.png", width = 10, height = 6, units = "in", res = 150); print(p1); dev.off()
png("gwas_haplotype_unfavorable_by_group.png", width = 10, height = 6, units = "in", res = 150); print(p2); dev.off()
cat("DONE\n")
