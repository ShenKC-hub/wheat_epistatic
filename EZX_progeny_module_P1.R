# EZX_progeny_module_P1.R
# EZX双亲群体(24GYL-181 × 23GZM-295)：将子代区间基因型转换为
# gwas_all_block_genotype.txt 框架的单倍型基因型, 计算 P1 上位效应网络
# 每个子代的模块覆盖数与各模块完整度
# 规则: 区间编码 0=23GZM-295型, 2=24GYL-181型, 1=杂合(忽略该区间相关的上位效应)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "EZX_progeny_module_P1"
dir.create(outdir, showWarnings = FALSE)
suppressMessages(library(data.table))

P295 <- "23GZM-295"; P181 <- "24GYL-181"

## ---- 1) 自然群体单倍型矩阵(父母本携带信息) ----
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
G <- as.matrix(apply(genotype, 2, as.numeric))
G <- G[, -1]
snames <- unlist(strsplit(readLines("gwas_all_block_genotype.txt", n=1), "\t"))[10:690]
rownames(G) <- snames

## ---- 2) P1 上位效应网络: 去重后的模块内有利边 ----
em <- read.delim("gwas_epistasis_edge_module_P1_txt")
em$snp1 <- as.character(em$snp1); em$snp2 <- as.character(em$snp2)
em$key <- paste(pmin(em$snp1, em$snp2), pmax(em$snp1, em$snp2), sep = "|")
em <- em[!duplicated(em$key), ]
em <- em[grepl("^P1_M", em$intra_module), ]
mod <- em$intra_module
Nm <- table(mod); mods <- names(Nm)
haps <- unique(c(em$snp1, em$snp2))          # 网络涉及的单倍型
cat("P1 网络: 边", nrow(em), " 节点单倍型", length(haps),
    " 模块", length(mods), " 模块内边", sum(Nm), "\n")

## ---- 3) 单倍型 -> 物理位置 -> EZX区间编码 ----
ann <- read.delim("all_block_annotated.txt", header = TRUE)
ann$block <- sub("_number\\.txt.*$", "", sub("^X", "", ann$block))
hap2block <- sub("_number\\.txt.*$", "", sub("^X", "", haps))
m <- match(hap2block, ann$block)
stopifnot(!any(is.na(m)))
hap_chr  <- ann$chr[m]
hap_mid  <- (ann$start[m] + ann$end[m]) / 2

w <- fread("/home/workspace/results/merge/EZX/second/adj/EZX_genotype_width.txt")
setkey(w, chr, region)
progs <- setdiff(names(w), c("chr", "region"))
w[, idx := seq_len(.N)]
cat("EZX 子代数:", length(progs), " 区间数:", nrow(w), "\n")

hap_idx <- integer(length(haps))
for (cc in unique(hap_chr)) {
  ii <- which(hap_chr == cc)
  wr <- w[chr == cc, region]
  hap_idx[ii] <- findInterval(hap_mid[ii], wr)   # 0=在第一个region之前(NA处理)
}
hap_i0 <- which(hap_idx == 0)
if (length(hap_i0)) cat("警告:", length(hap_i0), "个单倍型位置在EZX区间覆盖之外\n")
segM <- as.matrix(w[, ..progs])                   # 区间数 × 子代
seg_of_hap <- segM[hap_idx, , drop = FALSE]       # 单倍型 × 子代 (0/1/2)
dimnames(seg_of_hap) <- list(haps, progs)

## ---- 4) 每子代每边的携带判定 ----
p295 <- G[P295, haps] == 2      # 父母本携带(0/2)
p181 <- G[P181, haps] == 2
seg1 <- seg_of_hap[em$snp1, , drop = FALSE]   # 边 × 子代
seg2 <- seg_of_hap[em$snp2, , drop = FALSE]

p295_v <- as.integer(p295[em$snp1])   # 边水平: 3A...每条边一端单倍型是否被295携带
p181_v <- as.integer(p181[em$snp1])
q295_v <- as.integer(p295[em$snp2])
q181_v <- as.integer(p181[em$snp2])
c1 <- p295_v * (seg1 == 0) + p181_v * (seg1 == 2)
c2 <- q295_v * (seg2 == 0) + q181_v * (seg2 == 2)
c1[seg1 == 1] <- NA_integer_          # 杂合区间忽略
c2[seg2 == 1] <- NA_integer_
evaluable <- !is.na(c1) & !is.na(c2)
carried <- evaluable & (c1 == 1) & (c2 == 1)
dimnames(carried) <- dimnames(evaluable) <- list(edge = seq_len(nrow(em)), prog = progs)

## ---- 5) 每子代 × 每模块 汇总 ----
N_M    <- as.integer(Nm[mods])
n_ign  <- sapply(mods, function(mm) colSums(!evaluable[mod == mm, , drop = FALSE]))
n_eval <- sapply(mods, function(mm) colSums(evaluable[mod == mm, , drop = FALSE]))
n_car  <- sapply(mods, function(mm) colSums(carried[mod == mm, , drop = FALSE]))
rownames(n_ign) <- rownames(n_eval) <- rownames(n_car) <- progs
Comp_raw  <- sweep(n_car, 2, N_M, "/")
Comp_eff  <- sweep(n_car, 2, N_M, "/")
Comp_eff[n_eval == 0] <- NA

n_cov <- rowSums(n_car > 0)
stat <- data.frame(
  progeny = progs,
  n_edges_module_total = sum(N_M),
  n_edges_evaluable = as.integer(colSums(evaluable)),
  n_edges_ignored_hetero = as.integer(colSums(!evaluable)),
  n_edges_carried = as.integer(colSums(carried)),
  n_modules_total = length(mods),
  n_modules_covered = as.integer(n_cov),
  mean_completeness_raw = round(rowSums(Comp_raw * (n_car > 0)) / pmax(n_cov, 1), 4),
  mean_completeness_effective = round(rowSums(Comp_eff * (n_car > 0), na.rm = TRUE) / pmax(n_cov, 1), 4),
  overall_completeness_raw = round(rowSums(n_car) / sum(N_M), 4))
write.table(stat, file.path(outdir, "EZX_progeny_module_stats_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 6) 长表: 子代 × 模块 ----
p_idx <- rep(progs, times = length(mods))
m_idx <- factor(rep(mods, each = length(progs)), levels = mods)
long <- data.frame(
  progeny = p_idx, module = m_idx,
  n_edges_module = as.integer(Nm[m_idx]),
  n_edges_ignored_hetero = as.integer(n_ign[cbind(match(p_idx, rownames(n_ign)),
                                                    match(as.character(m_idx), colnames(n_ign)))]),
  n_edges_evaluable = as.integer(n_eval[cbind(match(p_idx, rownames(n_eval)),
                                              match(as.character(m_idx), colnames(n_eval)))]),
  n_edges_carried = as.integer(n_car[cbind(match(p_idx, rownames(n_car)),
                                           match(as.character(m_idx), colnames(n_car)))]))
long$completeness_raw = round(long$n_edges_carried / long$n_edges_module, 4)
long$completeness_effective = ifelse(long$n_edges_evaluable == 0, NA,
    round(long$n_edges_carried / long$n_edges_evaluable, 4))
long <- long[order(long$progeny, long$module), ]
write.table(long, file.path(outdir, "EZX_progeny_module_completeness_P1_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== 汇总 =====\n")
cat("子代数:", length(progs), "\n")
cat("n_edges_evaluable 范围:", range(stat$n_edges_evaluable), "\n")
cat("n_edges_ignored_hetero 范围:", range(stat$n_edges_ignored_hetero), "\n")
cat("n_edges_carried 范围:", range(stat$n_edges_carried), "\n")
cat("n_modules_covered 范围:", range(stat$n_modules_covered), "\n")
cat("mean_completeness_effective 中位:", median(stat$mean_completeness_effective, na.rm = TRUE), "\n")
cat("DONE\n")
