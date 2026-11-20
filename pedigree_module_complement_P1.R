# pedigree_module_complement_P1.R  (v2: 支持3亲本杂交)
# 系谱材料在 P1 总网络模块中的覆盖与完整度分析
# X1第一代: 23WZH-536 × 23WZH-562 × 23WZH-498(三亲本) → 23WZH-558
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "pedigree_module_P1"
dir.create(outdir, showWarnings = FALSE)

## ---- 基因型(全部681材料, 不做表型过滤) ----
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
G <- as.matrix(apply(genotype, 2, as.numeric))
G <- G[, -1]   # 第一列(名为"ID")是表头行残留
snames <- unlist(strsplit(readLines("gwas_all_block_genotype.txt", n=1), "\t"))[10:690]
rownames(G) <- snames
cat("samples:", nrow(G), " SNPs:", ncol(G), "\n")

## ---- P1 总网络模块(有利边, 去重, 真实模块) ----
em <- read.delim("gwas_epistasis_edge_module_P1_txt")
em$snp1 <- as.character(em$snp1); em$snp2 <- as.character(em$snp2)
em$key <- paste(pmin(em$snp1, em$snp2), pmax(em$snp1, em$snp2), sep = "|")
em <- em[!duplicated(em$key), ]
em <- em[grepl("^P1_M", em$intra_module), ]
mod <- em$intra_module
Nm <- table(mod); mods <- names(Nm)
cat("P1 总网络: 模块数 =", length(mods), " 模块内边数 =", sum(Nm), "\n")

## ---- 16份系谱材料 ----
mats <- c("23WZH-536","23WZH-562","23WZH-498","24GYL-181","23WZH-558","23HN-WZH-552","23WZH-559",
          "23HN-WZH-308","23WZH-263","23WZH-438",
          "23WZH-456","23WZH-626","23HN-BBC-157","23WZH-397",
          "24GYL-176","24GYL-177","23GZM-295")
srow <- match(mats, rownames(G))
stopifnot(!any(is.na(srow)))
carried <- sapply(srow, function(s) as.integer(G[s, em$snp1] == 2 & G[s, em$snp2] == 2))
colnames(carried) <- mats
Cm <- sapply(mods, function(m) colSums(carried[mod == m, , drop = FALSE]))
rownames(Cm) <- mats
Comp <- sweep(Cm, 2, as.integer(Nm[mods]), "/")
tot_edges_module <- sum(Nm)

## ---- 材料统计表 ----
stat <- data.frame(
  material = mats,
  n_edges_in_modules = as.integer(rowSums(Cm)),
  n_modules_covered  = as.integer(rowSums(Cm > 0)),
  n_modules_total    = length(mods),
  mean_completeness  = round(rowSums(Comp * (Cm > 0)) / pmax(rowSums(Cm > 0), 1), 4),
  overall_completeness = round(rowSums(Cm) / tot_edges_module, 4))
stat$mean_completeness[stat$n_modules_covered == 0] <- NA
write.table(stat, file.path(outdir, "pedigree_materials_module_stats_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(stat, row.names = FALSE, digits = 3)

## ---- 系谱组合(parents 用 | 分隔, 支持任意亲本数) ----
crosses <- data.frame(
  cross     = c("X1_layer1", "X1_layer2", "X2", "X3_off1", "X3_off2", "X5"),
  parents   = c("23WZH-536|23WZH-562|23WZH-498",
                "23WZH-558|23HN-WZH-552",
                "23HN-WZH-308|23WZH-263",
                "23WZH-456|23WZH-626",
                "23WZH-456|23WZH-626",
                "24GYL-176|24GYL-177"),
  offspring = c("23WZH-558","23WZH-559","23WZH-438",
                "23HN-BBC-157","23WZH-397","23GZM-295"),
  stringsAsFactors = FALSE)

sum_rows <- list(); mod_rows <- list()
for (k in seq_len(nrow(crosses))) {
  ps   <- strsplit(crosses$parents[k], "|", fixed = TRUE)[[1]]
  o    <- crosses$offspring[k]
  covP <- lapply(ps, function(pp) Cm[pp, ] > 0)
  compP<- lapply(ps, function(pp) Comp[pp, ])
  covO <- Cm[o, ] > 0; compO <- Comp[o, ]
  union <- Reduce(`|`, covP)
  maxpar <- Reduce(pmax, compP)
  ## 模块类别
  cat_mod <- rep("not_covered_by_all", length(mods))
  cat_mod[covO & !union] <- "de_novo_module"
  cat_mod[covO & union]  <- "some_parents_and_offspring"
  cat_mod[covO & Reduce(`&`, covP)] <- "all_parents_and_offspring"
  inc <- covO & union & (compO > maxpar + 1e-9)
  cat_mod[inc] <- "complementary_increase"
  ## 汇总
  pdet <- paste(sapply(seq_along(ps), function(j)
    sprintf("%s(%d模块/%d边)", ps[j], sum(covP[[j]]), sum(Cm[ps[j], ]))),
    collapse = "; ")
  sum_rows[[k]] <- data.frame(
    cross = crosses$cross[k], offspring = o, n_parents = length(ps),
    parents_detail = pdet,
    n_modules_union = sum(union), n_modules_offspring = sum(covO),
    n_modules_offspring_not_in_union = sum(covO & !union),
    n_modules_in_union_not_offspring = sum(union & !covO),
    n_edges_union = sum(Reduce(`|`, lapply(ps, function(pp) carried[, pp] == 1))),
    n_edges_offspring = sum(Cm[o, ]),
    n_modules_complementary_increase = sum(inc),
    n_modules_de_novo = sum(covO & !union))
  ## 模块明细(动态亲本列)
  mt <- data.frame(cross = crosses$cross[k], offspring = o, module = mods,
                   n_edges_module = as.integer(Nm[mods]),
                   comp_offspring = round(compO, 4),
                   covered_offspring = as.integer(covO),
                   n_parents_covered = as.integer(Reduce(`+`, covP)),
                   comp_max_parents = round(maxpar, 4))
  for (j in 1:3) {   # 最多3亲本, 缺失补NA
    mt[[paste0("comp_parent", j)]] <- if (j <= length(ps)) round(compP[[j]], 4) else NA_real_
    mt[[paste0("covered_parent", j)]] <- if (j <= length(ps)) as.integer(covP[[j]]) else NA_integer_
  }
  mt$parents <- crosses$parents[k]
  mt$category <- cat_mod
  mod_rows[[k]] <- mt
}
sum_tab <- do.call(rbind, sum_rows)
mod_tab <- do.call(rbind, mod_rows)
write.table(sum_tab, file.path(outdir, "pedigree_cross_summary_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mod_tab, file.path(outdir, "pedigree_cross_modules_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== 系谱组合汇总 =====\n")
print(sum_tab, row.names = FALSE, digits = 3)

inc_all <- mod_tab[mod_tab$category == "complementary_increase", ]
cat("\n互补提高完整度的模块总数:", nrow(inc_all), "\n")
if (nrow(inc_all) > 0) {
  u <- inc_all[!duplicated(paste(inc_all$cross, inc_all$module)), ]
  print(u[, c("cross","offspring","module","n_edges_module",
              "comp_parent1","comp_parent2","comp_parent3",
              "comp_offspring","comp_max_parents")], row.names = FALSE, digits = 3)
}

## 与 gwas_module_stats_P1_txt 一致性校验(有表型的材料)
old <- read.delim("gwas_module_stats_P1_txt")
chk <- merge(stat, old[, c("ID","n_modules_covered","mean_completeness")],
             by.x = "material", by.y = "ID")
cat("\n与旧结果(有表型材料)一致性: 覆盖模块数一致",
    sum(chk$n_modules_covered.x == chk$n_modules_covered.y), "/", nrow(chk), "\n")
cat("DONE\n")
