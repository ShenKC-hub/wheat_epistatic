# TOC_MSD_boxplot_P7.R
# SZY、BBC 群体：1A(TOC)/2D(MSD) 位点基因型 vs P7 表型箱线图 + 均值表
# 以及 1A×2D 基因型组合箱线图 + 均值表（基因型/表型获取同 double_pop.R）
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "TOC_MSD"
suppressMessages({library(ggplot2); library(dplyr)})

pops <- c("SZY", "BBC")
targets <- data.frame(chr = c("1A", "2D"),
                      start = c(369678240, 9078043),
                      end   = c(369712915, 9080881),
                      block = c("1A_block_20384(TOC)", "2D_block_177(MSD)"),
                      stringsAsFactors = FALSE)

plot_list_g <- list(); mean_list_g <- list()
plot_list_c <- list(); mean_list_c <- list()

for (p in pops) {
  df <- read.csv(paste0(p, "_qtl_blue.csv"), header = FALSE)
  df$V3 <- df$V3 * 1000000
  geno <- df[!is.na(df$V3), ]
  df_odd  <- geno[seq(1, nrow(geno), by = 2), ]
  df_even <- geno[(1:nrow(geno)) %% 2 == 0, ]
  df_odd <- cbind(df_odd[, 1:3], end = df_even$V3, df_odd[, 4:ncol(df_odd)])

  gvec <- list()
  for (k in 1:nrow(targets)) {
    d2 <- df_odd[df_odd$V2 == targets$chr[k] & df_odd$V3 <= targets$start[k] &
                 df_odd$end >= targets$end[k], ]
    d3 <- df_odd[df_odd$V2 == targets$chr[k] & df_odd$V3 >  targets$start[k] &
                 df_odd$V3 <  targets$end[k], ]
    d4 <- rbind(d2, d3)
    gvec[[k]] <- as.numeric(d4[1, 5:ncol(d4)])
  }
  g1 <- gvec[[1]]; g2 <- gvec[[2]]
  y <- as.numeric(t(df[7, 4:ncol(df)]))   # P7 = 第7行(性状)
  keep <- !is.na(y) & !is.na(g1) & !is.na(g2)
  y <- y[keep]; g1 <- g1[keep]; g2 <- g2[keep]
  n <- length(y)

  g1f <- factor(g1, levels = c(0, 1, 2), labels = c("AA(0)", "AB(1)", "BB(2)"))
  g2f <- factor(g2, levels = c(0, 1, 2), labels = c("AA(0)", "AB(1)", "BB(2)"))
  combo <- paste(g1f, g2f, sep = " × ")
  combo <- factor(combo, levels = c("AA(0) × AA(0)", "AA(0) × AB(1)", "AA(0) × BB(2)",
                                    "AB(1) × AA(0)", "AB(1) × AB(1)", "AB(1) × BB(2)",
                                    "BB(2) × AA(0)", "BB(2) × AB(1)", "BB(2) × BB(2)"))

  d1 <- data.frame(pop = p, site = "1A (TOC)", genotype = g1f, P7 = y)
  d2 <- data.frame(pop = p, site = "2D (MSD)", genotype = g2f, P7 = y)
  dc <- data.frame(pop = p, genotype = combo, P7 = y)

  plot_list_g[[p]] <- bind_rows(d1, d2)
  mean_list_g[[p]] <- bind_rows(
    d1 %>% group_by(pop, site, genotype) %>%
      summarise(n = n(), mean = mean(P7), sd = sd(P7), .groups = "drop"),
    d2 %>% group_by(pop, site, genotype) %>%
      summarise(n = n(), mean = mean(P7), sd = sd(P7), .groups = "drop"))
  mean_list_c[[p]] <- dc %>% group_by(pop, genotype) %>%
    summarise(n = n(), mean = mean(P7), sd = sd(P7), .groups = "drop")
  plot_list_c[[p]] <- dc

  tb1 <- table(g1f); tb2 <- table(g2f)
  cat(p, ": n =", n, " 1A分布:", paste(names(tb1), as.character(tb1), collapse = " "),
      "| 2D分布:", paste(names(tb2), as.character(tb2), collapse = " "), "\n")
}
plot_g <- bind_rows(plot_list_g)
mean_g <- bind_rows(mean_list_g)
plot_c <- bind_rows(plot_list_c)
mean_c <- bind_rows(mean_list_c)

## ---- 均值表 ----
mean_g_out <- mean_g %>%
  transmute(pop, site, genotype, n, mean = round(mean, 4), sd = round(sd, 4)) %>%
  mutate(mean_sd = sprintf("%.3f±%.3f", mean, sd))
mean_c_out <- mean_c %>%
  transmute(pop, genotype, n, mean = round(mean, 4), sd = round(sd, 4)) %>%
  mutate(mean_sd = sprintf("%.3f±%.3f", mean, sd)) %>%
  arrange(pop, genotype)
write.table(mean_g_out, file.path(outdir, "P7_single_site_means_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mean_c_out, file.path(outdir, "P7_combo_means_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 单位点箱线图 ----
pA <- ggplot(plot_g, aes(x = genotype, y = P7, fill = site)) +
  geom_boxplot(outlier.size = 1, width = 0.7, linewidth = 0.35) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6,
               fill = "white", colour = "black") +
  facet_grid(pop ~ site, scales = "free_y") +
  labs(title = "P7 phenotype by genotype at 1A(TOC) and 2D(MSD) in SZY and BBC",
       x = "Genotype", y = "P7 phenotype") +
  theme_bw() + theme(legend.position = "none",
                     strip.text = element_text(face = "bold"))
pdf(file.path(outdir, "P7_single_site_boxplot.pdf"), width = 9, height = 7); print(pA); dev.off()
png(file.path(outdir, "P7_single_site_boxplot.png"), width = 9, height = 7, units = "in", res = 150)
print(pA); dev.off()

## ---- 组合箱线图 ----
pC <- ggplot(plot_c, aes(x = genotype, y = P7, fill = pop)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4,
               fill = "white", colour = "black") +
  facet_wrap(~pop, nrow = 1, scales = "free") +
  labs(title = "P7 phenotype by 1A × 2D genotype combination",
       x = "1A genotype × 2D genotype", y = "P7 phenotype") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     strip.text = element_text(face = "bold"))
pdf(file.path(outdir, "P7_combo_boxplot.pdf"), width = 12, height = 5.5); print(pC); dev.off()
png(file.path(outdir, "P7_combo_boxplot.png"), width = 12, height = 5.5, units = "in", res = 150)
print(pC); dev.off()

cat("\n== 单位点均值表 ==\n"); print(mean_g_out, row.names = FALSE, digits = 3)
cat("\n== 组合均值表 ==\n"); print(mean_c_out, row.names = FALSE, digits = 3)
cat("DONE\n")
