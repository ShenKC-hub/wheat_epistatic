# overlap_181_295_P1.R
# P1 上位效应网络中 24GYL-181 与 23GZM-295 的模块及边重叠分析
# (EZX 双亲群体 = 24GYL-181 × 23GZM-295; 框架与 pedigree_module_complement_P1.R 一致)
suppressMessages({library(ggplot2); library(dplyr)})
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "pedigree_module_P1"

## ---- 基因型与P1网络(与系谱脚本一致) ----
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
G <- as.matrix(apply(genotype, 2, as.numeric))
G <- G[, -1]
snames <- unlist(strsplit(readLines("gwas_all_block_genotype.txt", n=1), "\t"))[10:690]
rownames(G) <- snames

em <- read.delim("gwas_epistasis_edge_module_P1_txt")
em$snp1 <- as.character(em$snp1); em$snp2 <- as.character(em$snp2)
em$key <- paste(pmin(em$snp1, em$snp2), pmax(em$snp1, em$snp2), sep = "|")
em <- em[!duplicated(em$key), ]
em <- em[grepl("^P1_M", em$intra_module), ]
mod <- em$intra_module
Nm <- table(mod); mods <- names(Nm)

A <- "24GYL-181"; B <- "23GZM-295"
srow <- match(c(A, B), rownames(G))
carA <- as.integer(G[srow[1], em$snp1] == 2 & G[srow[1], em$snp2] == 2)
carB <- as.integer(G[srow[2], em$snp1] == 2 & G[srow[2], em$snp2] == 2)

## ---- 边水平重叠 ----
nA <- sum(carA); nB <- sum(carB)
nboth <- sum(carA == 1 & carB == 1)
nunion <- sum(carA == 1 | carB == 1)
jac_e <- nboth / nunion
ovl_e <- nboth / min(nA, nB)          # 重叠系数
cat(sprintf("边重叠: %s携带%d, %s携带%d, 共同%d, 并集%d, Jaccard=%.3f, 重叠系数=%.3f\n",
            A, nA, B, nB, nboth, nunion, jac_e, ovl_e))

## ---- 模块水平重叠 ----
Cm <- cbind(A = as.vector(tapply(carA, mod, sum)), B = as.vector(tapply(carB, mod, sum)))
colnames(Cm) <- c(A, B); rownames(Cm) <- mods
covA <- Cm[, A] > 0; covB <- Cm[, B] > 0
mboth <- sum(covA & covB); munion <- sum(covA | covB)
jac_m <- mboth / munion
cat(sprintf("模块重叠: %s覆盖%d, %s覆盖%d, 共同%d, 并集%d, Jaccard=%.3f\n",
            A, sum(covA), B, sum(covB), mboth, munion, jac_m))

## ---- 逐模块明细 ----
eA <- as.integer(Cm[, A])
eB <- as.integer(Cm[, B])
eBoth <- sapply(mods, function(m) sum(carA[mod == m] == 1 & carB[mod == m] == 1))
eUni  <- sapply(mods, function(m) sum(carA[mod == m] == 1 | carB[mod == m] == 1))
tab <- data.frame(module = mods, n_edges_module = as.integer(Nm[mods]),
                  edges_181 = as.integer(eA), edges_295 = as.integer(eB),
                  edges_both = as.integer(eBoth), edges_union = as.integer(eUni))
tab$jaccard_edges = round(ifelse(tab$edges_union > 0, tab$edges_both / tab$edges_union, NA), 3)
tab$overlap_coeff = round(ifelse(tab$edges_both > 0, tab$edges_both / pmin(tab$edges_181, tab$edges_295), NA), 3)
tab$comp_181 = round(eA / tab$n_edges_module, 4)
tab$comp_295 = round(eB / tab$n_edges_module, 4)
tab$comp_union = round(eUni / tab$n_edges_module, 4)   # EZX后代可达到的潜力上限
tab$covered_181 = as.integer(covA)
tab$covered_295 = as.integer(covB)
## 结合EZX子代实际均值(来自EZX_progeny_module_completeness_P1_txt)
l <- read.delim(file.path("EZX_progeny_module_P1", "EZX_progeny_module_completeness_P1_txt"))
ms <- l %>% group_by(module) %>%
  summarise(ezx_mean_comp_eff = round(mean(completeness_effective, na.rm = TRUE), 4),
            ezx_mean_comp_raw = round(mean(completeness_raw, na.rm = TRUE), 4),
            .groups = "drop")
tab <- merge(tab, ms, by = "module", all.x = TRUE)
tab <- tab[order(-tab$n_edges_module), ]
write.table(tab, file.path(outdir, "overlap_181_295_modules_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 共同携带的边明细 ----
shared <- em[carA == 1 & carB == 1, c("snp1","snp2","intra_module","key")]
shared$module <- shared$intra_module
write.table(shared[, c("snp1","snp2","module")],
            file.path(outdir, "overlap_181_295_shared_edges_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== 逐模块明细(按边数降序) =====\n")
print(tab[, c("module","n_edges_module","edges_181","edges_295","edges_both",
              "jaccard_edges","comp_181","comp_295","comp_union","ezx_mean_comp_eff")],
      row.names = FALSE, digits = 3)
cat("\nDONE\n")
