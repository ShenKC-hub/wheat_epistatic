# qinben_progeny_favorable_hap.R
# 参考 qinben_progeny_module_all.R 的子代基因型获取逻辑:
#   子代区间基因型(0=亲本1型,2=亲本2型,1=杂合忽略) + 亲本在自然群体中的单倍型携带(0/2)
#   -> 子代是否携带某单倍型: seg==0 且 亲本1携带, 或 seg==2 且 亲本2携带
# 统计: 每个双亲群体每个子代, 对应不同性状(P1-P7)携带的 favorable 单倍型数目
#       (单倍型分类来自 gwas_haplotype_class_P*_txt 最后一列 class)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(data.table); library(dplyr)})
root <- "qinben_progeny_favorable_hap"
dir.create(root, showWarnings = FALSE)

## ---- 自然群体单倍型矩阵 (681 × 16102, 0/2, 行=材料含亲本) ----
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
G <- apply(genotype, 2, as.numeric)
rownames(G) <- unlist(strsplit(readLines("gwas_all_block_genotype.txt", n=1), "\t"))[10:690]
rm(genotype)

## ---- block 坐标注释 (单倍型 -> 物理位置) ----
ann <- read.delim("all_block_annotated.txt", header=TRUE)
ann$block <- sub("_number\\.txt.*$", "", sub("^X", "", ann$block))
ann <- ann[!duplicated(ann$block), ]

## ---- 22 个双亲群体 ----
pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
pops$p0 <- as.character(pops$p0); pops$p2 <- as.character(pops$p2)
stopifnot(all(c(pops$p0, pops$p2) %in% rownames(G)))

## ---- 每性状的单倍型分类 (favorable/unfavorable) ----
cls <- list(); all_haps <- character(0)
for (i in 1:7) {
  h <- read.delim(paste0("gwas_haplotype_class_P", i, "_txt"))
  h$HB <- as.character(h$HB)
  cls[[i]] <- h
  all_haps <- union(all_haps, h$HB)
}
all_haps <- unique(all_haps)
inG <- all_haps %in% colnames(G)
cat("union haps:", length(all_haps), "  not in G:", sum(!inG), "(dropped)\n")
all_haps <- all_haps[inG]

summary_all <- NULL
long_all <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]; P0 <- pops$p0[r]; P2 <- pops$p2[r]
  pdir <- file.path(root, pop); dir.create(pdir, showWarnings=FALSE)

  ## ---- 子代区间基因型 (来自 {POP}_genotype_width.txt) ----
  w <- fread(paste0(pop, "_genotype_width.txt"))
  setkey(w, chr, region)
  progs <- setdiff(names(w), c("chr","region"))
  w[, idx := seq_len(.N)]

  ## ---- 单倍型 -> block -> 区间索引 (与 qinben_progeny_module_all.R 一致) ----
  hap2block <- sub("_number\\.txt.*$", "", sub("^X", "", all_haps))
  m <- match(hap2block, ann$block)
  stopifnot(!any(is.na(m)))
  hap_chr <- ann$chr[m]; hap_mid <- (ann$start[m] + ann$end[m]) / 2
  hap_idx <- integer(length(all_haps))
  for (cc in unique(hap_chr)) {
    ii <- which(hap_chr == cc)
    hap_idx[ii] <- findInterval(hap_mid[ii], w[chr==cc, region])
  }
  segM <- as.matrix(w[, ..progs])
  seg_of_hap <- segM[hap_idx, , drop=FALSE]          # 单倍型 × 子代 (0/1/2)
  dimnames(seg_of_hap) <- list(all_haps, progs)
  rm(w, segM)

  ## ---- 亲本是否携带各单倍型 (自然群体 0/2) ----
  p0c <- setNames(G[P0, all_haps] == 2, all_haps)    # 亲本1携带
  p2c <- setNames(G[P2, all_haps] == 2, all_haps)    # 亲本2携带

  for (i in 1:7) {
    h <- cls[[i]]
    fav <- h$HB[h$class == "favorable"   & h$HB %in% all_haps]
    unf <- h$HB[h$class == "unfavorable" & h$HB %in% all_haps]

    ## 子代携带判定矩阵 (单倍型 × 子代)
    seg_f <- seg_of_hap[fav, , drop=FALSE]
    seg_u <- seg_of_hap[unf, , drop=FALSE]
    carry_f <- (seg_f == 0 & p0c[fav]) | (seg_f == 2 & p2c[fav])
    carry_u <- (seg_u == 0 & p0c[unf]) | (seg_u == 2 & p2c[unf])

    n_fav  <- as.integer(colSums(carry_f, na.rm=TRUE))
    n_unf  <- as.integer(colSums(carry_u, na.rm=TRUE))
    n_hetf <- as.integer(colSums(seg_f == 1, na.rm=TRUE))  # 杂合区间忽略的 favorable 单倍型数
    n_hetu <- as.integer(colSums(seg_u == 1, na.rm=TRUE))

    out <- data.frame(
      pop = pop,
      progeny = progs,
      trait = paste0("P", i),
      n_favorable = n_fav,
      n_favorable_total = length(fav),
      n_favorable_ignored_hetero = n_hetf,
      n_unfavorable = n_unf,
      n_unfavorable_total = length(unf),
      n_unfavorable_ignored_hetero = n_hetu)
    write.table(out, file.path(pdir, paste0(pop, "_progeny_favorable_P", i, "_txt")),
                sep="\t", quote=FALSE, row.names=FALSE)
    long_all <- rbind(long_all, out)

    summary_all <- rbind(summary_all, data.frame(
      pop = pop, trait = paste0("P", i),
      n_progeny = length(progs),
      n_fav_total = length(fav), n_unf_total = length(unf),
      mean_fav = round(mean(n_fav), 2), median_fav = as.integer(median(n_fav)),
      min_fav = min(n_fav), max_fav = max(n_fav),
      mean_unf = round(mean(n_unf), 2), median_unf = as.integer(median(n_unf))))
  }
  cat(sprintf("%s done (progeny=%d)\n", pop, length(progs)))
  rm(seg_of_hap)
}
write.table(summary_all, file.path(root, "summary_favorable_by_pop_trait_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
long_all <- long_all[order(long_all$pop, long_all$progeny, long_all$trait), ]
write.table(long_all, file.path(root, "all_progeny_favorable_counts_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: summary_favorable_by_pop_trait_txt / all_progeny_favorable_counts_txt\n")
cat("DONE\n")