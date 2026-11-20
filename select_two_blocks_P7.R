# select_two_blocks_P7.R
# 筛选 1A_block_20384 与 2D_block_177 在 P7 性状的上位效应对
# (1) 681份自然群体  (2) 25个双亲群体
# 保留这些互作对中已计算的两位点加性效应(alpha/p)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
outdir <- "select_1A20384_2D177_P7"
dir.create(outdir, showWarnings = FALSE)
TARGETS <- c("1A_block_20384", "2D_block_177")

## 自然群体单倍型ID -> block ID
block_of <- function(x) sub("_number\\.txt-Hap.*$", "", sub("^X", "", x))

## ============ (1) 681份自然群体 ============
eff <- read.delim("gwas_Epistasis_effect_P7_txt")
eff$snp1 <- as.character(eff$snp1); eff$snp2 <- as.character(eff$snp2)
eff$key <- paste(pmin(eff$snp1, eff$snp2), pmax(eff$snp1, eff$snp2), sep = "|")
eff <- eff[!duplicated(eff$key), ]
eff$b1 <- block_of(eff$snp1); eff$b2 <- block_of(eff$snp2)
hit <- eff$b1 %in% TARGETS | eff$b2 %in% TARGETS
nat <- eff[hit, ]

## 补充类别标签
ec <- read.delim("gwas_epistasis_edge_class_P7_txt")
ec$key <- paste(pmin(as.character(ec$snp1), as.character(ec$snp2)),
                pmax(as.character(ec$snp1), as.character(ec$snp2)), sep = "|")
m <- match(nat$key, ec$key)
nat$class <- ec$class[m]; nat$testP <- ec$testP[m]

## 将加性效应映射到目标block列（涉及两个目标block时各占一列）
nat$target_block <- apply(cbind(nat$b1, nat$b2), 1,
                          function(x) paste(TARGETS[TARGETS %in% x], collapse = ";"))
nat$alpha_1A_block_20384 <- NA; nat$p_1A_block_20384 <- NA
nat$alpha_2D_block_177   <- NA; nat$p_2D_block_177   <- NA
for (r in seq_len(nrow(nat))) {
  if (nat$b1[r] == "1A_block_20384") { nat$alpha_1A_block_20384[r] <- nat$alpha_snp1[r]; nat$p_1A_block_20384[r] <- nat$p_snp1[r] }
  if (nat$b2[r] == "1A_block_20384") { nat$alpha_1A_block_20384[r] <- nat$alpha_snp2[r]; nat$p_1A_block_20384[r] <- nat$p_snp2[r] }
  if (nat$b1[r] == "2D_block_177")   { nat$alpha_2D_block_177[r]   <- nat$alpha_snp1[r]; nat$p_2D_block_177[r]   <- nat$p_snp1[r] }
  if (nat$b2[r] == "2D_block_177")   { nat$alpha_2D_block_177[r]   <- nat$alpha_snp2[r]; nat$p_2D_block_177[r]   <- nat$p_snp2[r] }
}
nat_out <- nat[, c("snp1","snp2","target_block","pvalue","class","testP",
                   "beta_epi","se_epi","aaij","D","partialR2_epi",
                   "alpha_1A_block_20384","p_1A_block_20384",
                   "alpha_2D_block_177","p_2D_block_177",
                   "n_BiBj","m_BiBj")]
write.table(nat_out, file.path(outdir, "natural_681_P7_two_blocks_pairs_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("自然群体: 命中互作对(去重) =", nrow(nat_out), "\n")
print(table(nat_out$target_block))

## ============ (2) 25个双亲群体 ============
pops <- c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
          "JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK")
bp_list <- list()
for (p in pops) {
  f_eff <- paste0("Epistasis_effect_", p, "_P7_txt")
  f_cls <- file.path("biparental", p, "Epistasis_edge_class_P7_txt")
  if (!file.exists(f_eff)) next
  e <- read.delim(f_eff)
  e$snp1 <- as.character(e$snp1); e$snp2 <- as.character(e$snp2)
  e$key <- paste(pmin(e$snp1, e$snp2), pmax(e$snp1, e$snp2), sep = "|")
  e <- e[!duplicated(e$key), ]
  hit <- e$snp1 %in% TARGETS | e$snp2 %in% TARGETS
  if (!any(hit)) next
  s <- e[hit, ]
  s$pop <- p
  if (file.exists(f_cls)) {
    cc <- read.delim(f_cls)
    cc$key <- paste(pmin(as.character(cc$snp1), as.character(cc$snp2)),
                    pmax(as.character(cc$snp1), as.character(cc$snp2)), sep = "|")
    m <- match(s$key, cc$key)
    s$class <- cc$class[m]; s$testP <- cc$testP[m]
  }
  s$alpha_1A_block_20384 <- NA; s$p_1A_block_20384 <- NA
  s$alpha_2D_block_177   <- NA; s$p_2D_block_177   <- NA
  for (r in seq_len(nrow(s))) {
    if (s$snp1[r] == "1A_block_20384") { s$alpha_1A_block_20384[r] <- s$alpha_snp1[r]; s$p_1A_block_20384[r] <- s$p_snp1[r] }
    if (s$snp2[r] == "1A_block_20384") { s$alpha_1A_block_20384[r] <- s$alpha_snp2[r]; s$p_1A_block_20384[r] <- s$p_snp2[r] }
    if (s$snp1[r] == "2D_block_177")   { s$alpha_2D_block_177[r]   <- s$alpha_snp1[r]; s$p_2D_block_177[r]   <- s$p_snp1[r] }
    if (s$snp2[r] == "2D_block_177")   { s$alpha_2D_block_177[r]   <- s$alpha_snp2[r]; s$p_2D_block_177[r]   <- s$p_snp2[r] }
  }
  bp_list[[p]] <- s
}
bp <- do.call(rbind, bp_list)
bp$target_block <- apply(cbind(bp$snp1, bp$snp2), 1,
                         function(x) paste(TARGETS[TARGETS %in% x], collapse = ";"))
bp_out <- bp[, c("pop","snp1","snp2","target_block","pvalue","class","testP",
                 "beta_epi","se_epi","aaij","D","partialR2_epi",
                 "alpha_1A_block_20384","p_1A_block_20384",
                 "alpha_2D_block_177","p_2D_block_177",
                 "n_BiBj","m_BiBj")]
write.table(bp_out, file.path(outdir, "biparental_25pops_P7_two_blocks_pairs_txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n双亲群体: 命中互作对(去重) =", nrow(bp_out), "  涉及群体数 =", length(unique(bp_out$pop)), "\n")
cat("\n各群体命中数:\n")
print(table(bp_out$pop))
cat("\n与两个目标block配对的伙伴block(Top):\n")
partner <- c(bp_out$snp1[bp_out$snp2 %in% TARGETS], bp_out$snp2[bp_out$snp1 %in% TARGETS])
print(sort(table(partner), decreasing = TRUE)[1:10])
cat("\n显著(P<0.05)的组合:\n")
sig <- bp_out[bp_out$pvalue < 0.05, c("pop","snp1","snp2","pvalue","beta_epi","aaij")]
if (nrow(sig) > 0) print(sig, row.names = FALSE, digits = 3) else cat("无\n")
cat("DONE\n")
