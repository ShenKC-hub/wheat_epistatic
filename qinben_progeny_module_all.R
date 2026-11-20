# qinben_progeny_module_all.R
# 22个双亲群体 × P1-P7：
#   (1) 参照 overlap_181_295_P1.R：两个亲本在自然群体有利上位网络中的
#       边/模块重叠（互补性），逐模块完整度 comp_A/comp_B/comp_union
#   (2) 参照 EZX_progeny_module_P1.R：子代区间基因型(0=父本1型,2=父本2型,1=杂合忽略)
#       -> 每子代每模块的携带边数与完整度(raw/effective)
# 规则与参考代码一致: qinben_group2.txt 第2列亲本=0, 第3列亲本=2
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(data.table); library(dplyr)})
root <- "qinben_progeny_module"
dir.create(root, showWarnings = FALSE)

## ---- 自然群体单倍型矩阵（681 × 16102, 0/2）----
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <- t(genotype[, -c(1:2, 4:9)])
colnames(genotype) <- genotype[1, ]; genotype <- genotype[-1, ]
G <- apply(genotype, 2, as.numeric)
rownames(G) <- unlist(strsplit(readLines("gwas_all_block_genotype.txt", n=1), "\t"))[10:690]
rm(genotype)

## ---- block 坐标注释（单倍型 -> 物理位置）----
ann <- read.delim("all_block_annotated.txt", header=TRUE)
ann$block <- sub("_number\\.txt.*$", "", sub("^X", "", ann$block))
ann <- ann[!duplicated(ann$block), ]

## ---- 22 个群体 ----
pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
pops$p0 <- as.character(pops$p0); pops$p2 <- as.character(pops$p2)
stopifnot(all(c(pops$p0, pops$p2) %in% rownames(G)))

summary_all <- NULL
for (i in 1:7) {
  ## ---- P{i} 网络: 去重后的模块内边 ----
  em <- read.delim(paste0("gwas_epistasis_edge_module_P", i, "_txt"))
  em$snp1 <- as.character(em$snp1); em$snp2 <- as.character(em$snp2)
  em$key <- paste(pmin(em$snp1, em$snp2), pmax(em$snp1, em$snp2), sep="|")
  em <- em[!duplicated(em$key), ]
  em <- em[grepl(paste0("^P", i, "_M"), em$intra_module), ]
  mod <- em$intra_module
  Nm <- table(mod); mods <- names(Nm)
  haps <- unique(c(em$snp1, em$snp2))
  hap2block <- sub("_number\\.txt.*$", "", sub("^X", "", haps))
  m <- match(hap2block, ann$block)
  stopifnot(!any(is.na(m)))
  hap_chr <- ann$chr[m]; hap_mid <- (ann$start[m] + ann$end[m]) / 2
  cat(sprintf("P%d: edges=%d nodes=%d modules=%d\n", i, nrow(em), length(haps), length(mods)))

  for (r in seq_len(nrow(pops))) {
    pop <- pops$pop[r]; P0 <- pops$p0[r]; P2 <- pops$p2[r]
    pdir <- file.path(root, pop); dir.create(pdir, showWarnings=FALSE)

    ## ---- (1) 亲本在边/模块水平的重叠(互补) ----
    carA <- as.integer(G[P0, em$snp1] == 2 & G[P0, em$snp2] == 2)
    carB <- as.integer(G[P2, em$snp1] == 2 & G[P2, em$snp2] == 2)
    nA <- sum(carA); nB <- sum(carB)
    nboth <- sum(carA==1 & carB==1); nunion <- sum(carA==1 | carB==1)
    jac_e <- nboth/nunion; ovl_e <- nboth/max(min(nA,nB), 1)

    Cm <- cbind(A=as.vector(tapply(carA, mod, sum)), B=as.vector(tapply(carB, mod, sum)))
    colnames(Cm) <- c(P0, P2); rownames(Cm) <- mods
    covA <- Cm[,1] > 0; covB <- Cm[,2] > 0
    mboth <- sum(covA & covB); munion <- sum(covA | covB)

    eA <- as.integer(Cm[,1]); eB <- as.integer(Cm[,2])
    eBoth <- sapply(mods, function(mm) sum(carA[mod==mm]==1 & carB[mod==mm]==1))
    eUni  <- sapply(mods, function(mm) sum(carA[mod==mm]==1 | carB[mod==mm]==1))
    tab <- data.frame(module=mods, n_edges_module=as.integer(Nm[mods]),
                      edges_p0=eA, edges_p2=eB, edges_both=eBoth, edges_union=eUni)
    tab$jaccard_edges <- round(ifelse(tab$edges_union>0, tab$edges_both/tab$edges_union, NA), 3)
    tab$overlap_coeff <- round(ifelse(tab$edges_both>0, tab$edges_both/pmin(tab$edges_p0,tab$edges_p2), NA), 3)
    tab$comp_p0   <- round(eA/tab$n_edges_module, 4)
    tab$comp_p2   <- round(eB/tab$n_edges_module, 4)
    tab$comp_union<- round(eUni/tab$n_edges_module, 4)   # 子代潜力上限
    tab$covered_p0 <- as.integer(covA); tab$covered_p2 <- as.integer(covB)
    tab <- tab[order(-tab$n_edges_module), ]
    write.table(tab, file.path(pdir, paste0(pop, "_overlap_modules_P", i, "_txt")),
                sep="\t", quote=FALSE, row.names=FALSE)

    shared <- em[carA==1 & carB==1, c("snp1","snp2","intra_module")]
    write.table(shared, file.path(pdir, paste0(pop, "_shared_edges_P", i, "_txt")),
                sep="\t", quote=FALSE, row.names=FALSE)

    ## ---- (2) 子代模块覆盖与完整度 ----
    w <- fread(paste0(pop, "_genotype_width.txt"))
    setkey(w, chr, region)
    progs <- setdiff(names(w), c("chr","region"))
    w[, idx := seq_len(.N)]
    hap_idx <- integer(length(haps))
    for (cc in unique(hap_chr)) {
      ii <- which(hap_chr == cc)
      hap_idx[ii] <- findInterval(hap_mid[ii], w[chr==cc, region])
    }
    segM <- as.matrix(w[, ..progs])
    seg_of_hap <- segM[hap_idx, , drop=FALSE]
    dimnames(seg_of_hap) <- list(haps, progs)

    p0c <- G[P0, haps] == 2      # 父本1携带(0/2)
    p2c <- G[P2, haps] == 2      # 父本2携带
    seg1 <- seg_of_hap[em$snp1, , drop=FALSE]
    seg2 <- seg_of_hap[em$snp2, , drop=FALSE]
    c1 <- p0c[em$snp1]*(seg1==0) + p2c[em$snp1]*(seg1==2)
    c2 <- p0c[em$snp2]*(seg2==0) + p2c[em$snp2]*(seg2==2)
    c1[seg1==1] <- NA_integer_; c2[seg2==1] <- NA_integer_
    evaluable <- !is.na(c1) & !is.na(c2)
    carried <- evaluable & (c1==1) & (c2==1)
    dimnames(carried) <- dimnames(evaluable) <- list(edge=seq_len(nrow(em)), prog=progs)

    N_M <- as.integer(Nm[mods])
    n_ign  <- sapply(mods, function(mm) colSums(!evaluable[mod==mm, , drop=FALSE]))
    n_eval <- sapply(mods, function(mm) colSums( evaluable[mod==mm, , drop=FALSE]))
    n_car  <- sapply(mods, function(mm) colSums( carried[mod==mm, , drop=FALSE]))
    rownames(n_ign) <- rownames(n_eval) <- rownames(n_car) <- progs
    Comp_raw <- sweep(n_car, 2, N_M, "/")
    Comp_eff <- sweep(n_car, 2, N_M, "/"); Comp_eff[n_eval==0] <- NA
    n_cov <- rowSums(n_car > 0)

    stat <- data.frame(
      progeny=progs,
      n_edges_module_total=sum(N_M),
      n_edges_evaluable=as.integer(colSums(evaluable)),
      n_edges_ignored_hetero=as.integer(colSums(!evaluable)),
      n_edges_carried=as.integer(colSums(carried)),
      n_modules_total=length(mods),
      n_modules_covered=as.integer(n_cov),
      mean_completeness_raw=round(rowSums(Comp_raw*(n_car>0))/pmax(n_cov,1), 4),
      mean_completeness_effective=round(rowSums(Comp_eff*(n_car>0), na.rm=TRUE)/pmax(n_cov,1), 4),
      overall_completeness_raw=round(rowSums(n_car)/sum(N_M), 4))
    write.table(stat, file.path(pdir, paste0(pop, "_progeny_stats_P", i, "_txt")),
                sep="\t", quote=FALSE, row.names=FALSE)

    p_idx <- rep(progs, times=length(mods))
    m_idx <- factor(rep(mods, each=length(progs)), levels=mods)
    ii_p <- match(p_idx, progs); ii_m <- match(as.character(m_idx), mods)
    long <- data.frame(
      progeny=p_idx, module=m_idx,
      n_edges_module=as.integer(Nm[m_idx]),
      n_edges_ignored_hetero=as.integer(n_ign[cbind(ii_p, ii_m)]),
      n_edges_evaluable=as.integer(n_eval[cbind(ii_p, ii_m)]),
      n_edges_carried=as.integer(n_car[cbind(ii_p, ii_m)]))
    long$completeness_raw <- round(long$n_edges_carried/long$n_edges_module, 4)
    long$completeness_effective <- ifelse(long$n_edges_evaluable==0, NA,
        round(long$n_edges_carried/long$n_edges_evaluable, 4))
    long <- long[order(long$progeny, long$module), ]
    write.table(long, file.path(pdir, paste0(pop, "_progeny_completeness_P", i, "_txt")),
                sep="\t", quote=FALSE, row.names=FALSE)

    summary_all <- rbind(summary_all, data.frame(
      pop=pop, trait=paste0("P", i),
      n_progeny=length(progs),
      edges_p0=nA, edges_p2=nB, edges_both=nboth, edges_union=nunion,
      jaccard_edges=round(jac_e,3), overlap_coeff=round(ovl_e,3),
      modules_p0=sum(covA), modules_p2=sum(covB),
      modules_both=mboth, modules_union=munion,
      jaccard_modules=round(mboth/max(munion,1),3),
      prog_median_carried=as.integer(median(stat$n_edges_carried)),
      prog_median_cov=as.integer(median(stat$n_modules_covered)),
      prog_median_comp_eff=round(median(stat$mean_completeness_effective, na.rm=TRUE),4),
      prog_median_comp_raw=round(median(stat$mean_completeness_raw, na.rm=TRUE),4)))
    cat(sprintf("  %s: edges %d/%d both=%d J=%.3f | prog median comp_eff=%.3f\n",
                pop, nA, nB, nboth, jac_e,
                median(stat$mean_completeness_effective, na.rm=TRUE)))
  }
  saveRDS(NULL, tempfile())
}
write.table(summary_all, file.path(root, "summary_all_pops_traits_txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\nsaved: qinben_progeny_module/summary_all_pops_traits_txt\n")
cat("DONE\n")