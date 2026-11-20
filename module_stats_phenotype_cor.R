# module_stats_phenotype_cor.R
# 每个性状：n_modules_covered / mean_completeness / overall_completeness 与表型的相关性
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr)})

phe <- read.delim("phe2.txt", header = FALSE, comment.char = "#")
phe$ID <- as.character(phe$V1); phe$V1 <- NULL
names(phe) <- c(paste0("P", 1:7), "ID")
phe_l <- reshape(phe, direction = "long", varying = paste0("P", 1:7),
                 v.names = "phenotype", times = paste0("P", 1:7),
                 timevar = "trait", idvar = "ID")
phe_l$trait <- factor(phe_l$trait, levels = paste0("P", 1:7))

dat <- list(); cor_tab <- NULL
for (i in 1:7) {
  st <- read.delim(paste0("gwas_module_stats_P", i, "_txt"))
  st$trait <- paste0("P", i)
  d <- merge(st, phe_l, by = c("ID", "trait"))
  d <- d[!is.na(d$phenotype), ]
  dat[[i]] <- d
  for (v in c("n_modules_covered", "mean_completeness", "overall_completeness")) {
    ct <- cor.test(d[[v]], d$phenotype, method = "pearson")
    stt <- cor.test(d[[v]], d$phenotype, method = "spearman", exact = FALSE)
    cor_tab <- rbind(cor_tab, data.frame(
      trait = paste0("P", i), variable = v, n = nrow(d),
      pearson_r = unname(ct$estimate), pearson_p = ct$p.value,
      spearman_rho = unname(stt$estimate), spearman_p = stt$p.value))
  }
}
write.table(cor_tab, "gwas_module_phenotype_cor_txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

dat_all <- do.call(rbind, dat)
dat_all$trait <- factor(dat_all$trait, levels = paste0("P", 1:7))

## ---- 图A: n_modules_covered vs 表型 ----
pA <- ggplot(dat_all, aes(x = n_modules_covered, y = phenotype)) +
  geom_point(shape = 21, fill = "#2C64AE", colour = "grey25",
             alpha = 0.55, size = 1.6, stroke = 0.3) +
  geom_smooth(method = "loess", se = TRUE, colour = "#DF0101",
              fill = "#DF0101", alpha = 0.12, linewidth = 0.8) +
  facet_wrap(~trait, nrow = 2, scales = "free") +
  labs(title = "Number of modules covered vs phenotype",
       x = "n_modules_covered", y = "Phenotype value") +
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     panel.grid.minor = element_blank())

## ---- 图B: mean_completeness vs 表型 ----
pB <- ggplot(dat_all, aes(x = mean_completeness, y = phenotype)) +
  geom_point(shape = 21, fill = "#04B404", colour = "grey25",
             alpha = 0.55, size = 1.6, stroke = 0.3) +
  geom_smooth(method = "loess", se = TRUE, colour = "#DF0101",
              fill = "#DF0101", alpha = 0.12, linewidth = 0.8) +
  facet_wrap(~trait, nrow = 2, scales = "free") +
  labs(title = "Mean module completeness vs phenotype",
       x = "Mean completeness of covered modules", y = "Phenotype value") +
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     panel.grid.minor = element_blank())

## ---- 图C: overall_completeness vs 表型 ----
pC <- ggplot(dat_all, aes(x = overall_completeness, y = phenotype)) +
  geom_point(shape = 21, fill = "#A945FF", colour = "grey25",
             alpha = 0.55, size = 1.6, stroke = 0.3) +
  geom_smooth(method = "loess", se = TRUE, colour = "#DF0101",
              fill = "#DF0101", alpha = 0.12, linewidth = 0.8) +
  facet_wrap(~trait, nrow = 2, scales = "free") +
  labs(title = "Overall module completeness vs phenotype",
       x = "Overall completeness (edges carried / total module edges)",
       y = "Phenotype value") +
  theme_bw() + theme(strip.text = element_text(face = "bold"),
                     panel.grid.minor = element_blank())

pdf("gwas_module_coverage_phenotype_scatter.pdf", width = 12, height = 7); print(pA); dev.off()
pdf("gwas_module_completeness_phenotype_scatter.pdf", width = 12, height = 7); print(pB); dev.off()
pdf("gwas_module_overall_completeness_phenotype_scatter.pdf", width = 12, height = 7); print(pC); dev.off()
png("gwas_module_coverage_phenotype_scatter.png", width = 12, height = 7, units = "in", res = 150); print(pA); dev.off()
png("gwas_module_completeness_phenotype_scatter.png", width = 12, height = 7, units = "in", res = 150); print(pB); dev.off()
png("gwas_module_overall_completeness_phenotype_scatter.png", width = 12, height = 7, units = "in", res = 150); print(pC); dev.off()

cat("\n===== 相关性汇总 (Pearson / Spearman) =====\n")
ctab <- cor_tab
ctab$pearson_p <- signif(ctab$pearson_p, 2); ctab$spearman_p <- signif(ctab$spearman_p, 2)
print(ctab, row.names = FALSE, digits = 3)
cat("DONE\n")
