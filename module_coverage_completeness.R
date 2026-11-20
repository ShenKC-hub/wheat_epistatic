# module_coverage_completeness.R
# 以全局有利网络模块(gwas_epistasis_edge_module_P*_txt)为固定框架：
#   1) 每个材料覆盖的模块数 n_modules_covered
#   2) 每材料×每模块的完整度 = 携带边数/模块总边数
# 仅统计有利(favorable)边；非模块(<5节点小分量)边不统计
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

grp <- read.delim("sample_group.txt", header = FALSE,
                  col.names = c("ID", "group"), comment.char = "#")
grp$ID <- as.character(grp$ID)
grp_order <- c("Landrace", "Introduced", "Pre-1960",
               "1961-1980", "1981-2000", "2001-2020")

stats_all <- list()
for (i in 1:7) {
  PA <- readRDS(paste0("gwas_epistasis_PA_P", i, ".rds"))
  em <- read.delim(paste0("gwas_epistasis_edge_module_P", i, "_txt"))
  em$snp1 <- as.character(em$snp1); em$snp2 <- as.character(em$snp2)

  ## 去重（A-B/B-A 重复对）并仅保留真实模块内的边
  em$key <- paste(pmin(em$snp1, em$snp2), pmax(em$snp1, em$snp2), sep = "|")
  em <- em[!duplicated(em$key), ]
  keep <- grepl(paste0("^P", i, "_M"), em$intra_module)
  em <- em[keep, ]
  mod <- em$intra_module
  Nm  <- table(mod)                       # 每模块总边数(去重)
  mods <- names(Nm)

  ## PA列按key匹配(两种朝向)
  lut <- setNames(seq_along(colnames(PA)), colnames(PA))
  idx <- lut[paste(em$snp1, em$snp2)]
  rv  <- lut[paste(em$snp2, em$snp1)]
  idx[is.na(idx)] <- rv[is.na(idx)]
  stopifnot(!any(is.na(idx)))
  PAf <- PA[, idx, drop = FALSE]

  ## 每材料×每模块的携带边数与完整度
  C <- sapply(mods, function(m) as.integer(rowSums(PAf[, mod == m, drop = FALSE])))
  rownames(C) <- rownames(PA)
  Comp <- sweep(C, 2, as.integer(Nm[mods]), "/")

  ## 长表输出
  long <- as.data.frame(table(
    ID = factor(rep(rownames(C), ncol(C)), levels = rownames(C)),
    module = factor(rep(mods, each = nrow(C)), levels = mods)))
  names(long) <- c("ID", "module")
  long$n_edges_module   <- as.integer(Nm[long$module])
  long$n_edges_carried  <- as.vector(C)
  long$completeness     <- round(as.vector(Comp), 4)
  long$trait <- paste0("P", i)
  write.table(long, paste0("gwas_module_completeness_P", i, "_txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## 材料汇总表
  n_cov <- rowSums(C > 0)
  mean_comp <- rowSums(Comp * (C > 0)) / pmax(n_cov, 1)
  mean_comp[n_cov == 0] <- NA
  st <- data.frame(ID = rownames(C), trait = paste0("P", i),
                   n_modules_total = length(mods),
                   n_modules_covered = n_cov,
                   total_module_edges = sum(Nm),
                   n_edges_carried_in_modules = rowSums(C),
                   mean_completeness = round(mean_comp, 4))
  st$overall_completeness <- round(st$n_edges_carried_in_modules / st$total_module_edges, 4)
  write.table(st, paste0("gwas_module_stats_P", i, "_txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  stats_all[[i]] <- st
  cat("P", i, ": modules =", length(mods), " module edges =", sum(Nm),
      " covered(中位) =", median(n_cov), " mean_completeness(中位) =",
      round(median(mean_comp, na.rm = TRUE), 3), "\n")
}
stats_all <- do.call(rbind, stats_all)
stats_all <- merge(stats_all, grp, by = "ID", all.x = TRUE)
stats_all$trait <- factor(stats_all$trait, levels = paste0("P", 1:7))
stats_all$group <- factor(stats_all$group, levels = grp_order)

## ---- 育种进程分组统计 ----
summ <- stats_all %>%
  filter(!is.na(group)) %>%
  group_by(trait, group) %>%
  summarise(n = n(),
            cov_mean = mean(n_modules_covered), cov_sd = sd(n_modules_covered),
            comp_mean = mean(mean_completeness, na.rm = TRUE),
            comp_sd = sd(mean_completeness, na.rm = TRUE),
            overall_comp_mean = mean(overall_completeness),
            .groups = "drop") %>% as.data.frame()
summ$cov_ms  <- sprintf("%.1f±%.1f", summ$cov_mean, summ$cov_sd)
summ$comp_ms <- sprintf("%.3f±%.3f", summ$comp_mean, summ$comp_sd)
out_wide <- summ %>%
  transmute(trait, group, n,
            cov_mean, cov_sd, cov_mean_sd = cov_ms,
            comp_mean, comp_sd, comp_mean_sd = comp_ms,
            overall_comp_mean) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
write.table(out_wide, "gwas_module_stats_by_group_txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

pd <- position_dodge(0.8)
p1 <- ggplot(summ, aes(x = group, y = cov_mean, fill = trait)) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = cov_mean - cov_sd, ymax = cov_mean + cov_sd),
                width = 0.25, linewidth = 0.3) +
  facet_wrap(~trait, nrow = 2, scales = "free_y") +
  labs(title = "Number of modules covered by accessions across breeding eras",
       x = "Breeding era", y = "Mean n_modules_covered (±SD)", fill = "Trait") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     strip.text = element_text(face = "bold"))
p2 <- ggplot(summ, aes(x = group, y = comp_mean, fill = trait)) +
  geom_bar(stat = "identity", width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = comp_mean - comp_sd, ymax = comp_mean + comp_sd),
                width = 0.25, linewidth = 0.3) +
  facet_wrap(~trait, nrow = 2, scales = "free_y") +
  labs(title = "Module completeness of accessions across breeding eras",
       x = "Breeding era", y = "Mean module completeness (±SD)", fill = "Trait") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     strip.text = element_text(face = "bold"))

pdf("gwas_module_coverage_by_group.pdf", width = 11, height = 7); print(p1); dev.off()
pdf("gwas_module_completeness_by_group.pdf", width = 11, height = 7); print(p2); dev.off()
png("gwas_module_coverage_by_group.png", width = 11, height = 7, units = "in", res = 150)
print(p1); dev.off()
png("gwas_module_completeness_by_group.png", width = 11, height = 7, units = "in", res = 150)
print(p2); dev.off()
cat("DONE\n")
