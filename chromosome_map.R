# chromosome_map.R  (v2: 单行横版)
# 21条染色体一行排列；GWAS(圆点)/QTL(方块)/已知基因(三角+名称)紧凑三列，竖线分隔
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")

genome <- read.delim("../../../mygenome.txt", header = TRUE)
colnames(genome) <- c("chr", "start", "end", "num")
genome <- genome[genome$chr != "Un", ]
bl <- read.delim("all_block_annotated.txt", header = TRUE, check.names = FALSE)
bl <- bl[bl$chr != "Un", ]

trait_col <- c("#2121D9", "#9999FF", "#DF0101", "#04B404",
               "#FFFB23", "#FF9326", "#A945FF")
gene_col  <- "#000000"          # 基因名黑色
SEP_COL   <- "grey55"
chrs <- paste0(rep(1:7, each = 3), c("A", "B", "D"))

## 槽内坐标（每条染色体占1个横向单位，整体xlim=0..21）
BAR_HW  <- 0.015                       # 染色体半宽（约为旧版1/10）
X_GW    <- c(0.075, 0.097)             # GWAS 2子列，几乎紧挨
SEP1    <- 0.125                       # GWAS|QTL 分割线
X_QTL   <- 0.145 + (0:5) * 0.022       # QTL 6子列，几乎紧挨
SEP2    <- 0.288                       # QTL|基因 分割线
X_GENE  <- 0.315
X_TEXT  <- 0.345
CEX_PT  <- 0.22
CEX_G   <- 0.26
CEX_TXT <- 0.69
GAP     <- 7.2                         # 基因名纵向最小间距(Mb)

draw_map <- function() {
par(mar = c(4, 4.5, 3.5, 2))
plot(NA, xlim = c(0, 21), ylim = c(0, 960), axes = FALSE, xlab = "", ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)

## 顶部图例
text(10.5, 935, "Chromosomal distribution of HB-based GWAS signals, QTL and known genes",
     cex = 1.3, font = 2, xpd = TRUE)
legend(x = 0.2, y = 892, bty = "n", ncol = 3, xpd = TRUE, cex = 1.0,
       legend = c("GWAS (circle)", "QTL (square)", "Known gene (triangle)"),
       pch = c(21, 22, 17), pt.cex = 1.2,
       col = c("grey30", "grey30", gene_col), pt.bg = c("grey30", "grey30", gene_col))
legend(x = 5.2, y = 892, bty = "n", ncol = 7, title = "Trait", xpd = TRUE, cex = 1.0,
       legend = paste0("P", 1:7), pch = 21, pt.cex = 1.2, pt.bg = trait_col, col = trait_col)

## y轴
axis(2, at = seq(0, 800, by = 100), las = 1, cex.axis = 0.75, tck = -0.008, lwd = 0.8)
mtext("Physical position (Mb)", side = 2, line = 2.6, cex = 0.9)

for (k in seq_along(chrs)) {
  chr <- chrs[k]
  len <- genome$end[genome$chr == chr]
  blc <- bl[bl$chr == chr, ]
  x0  <- k - 1                          # 槽左边界

  ## 染色体（细条）
  rect(x0 + 0.045 - BAR_HW, 0, x0 + 0.045 + BAR_HW, len,
       col = "grey93", border = "grey35", lwd = 1)

  ## 分割竖线
  segments(x0 + SEP1, 0, x0 + SEP1, len, col = SEP_COL, lwd = 0.5)
  segments(x0 + SEP2, 0, x0 + SEP2, len, col = SEP_COL, lwd = 0.5)

  gw_idx <- as.matrix(blc[, paste0("GWAS_P", 1:7)])
  qt_idx <- as.matrix(blc[, paste0("QTL_P", 1:7)])
  gene_rows <- integer(0)

  for (r in seq_len(nrow(blc))) {
    y <- (blc$start[r] + blc$end[r]) / 2 / 1e6
    gw <- which(gw_idx[r, ] == 1)
    if (length(gw)) points(x0 + X_GW[seq_along(gw)], rep(y, length(gw)),
                           pch = 21, cex = CEX_PT, lwd = 0.4,
                           col = trait_col[gw], bg = trait_col[gw])
    qt <- which(qt_idx[r, ] == 1)
    if (length(qt)) points(x0 + X_QTL[seq_along(qt)], rep(y, length(qt)),
                           pch = 22, cex = CEX_PT, lwd = 0.4,
                           col = trait_col[qt], bg = trait_col[qt])
    if (blc$known_gene[r] != "") {
      points(x0 + X_GENE, y, pch = 17, cex = CEX_G, col = gene_col)
      gene_rows <- c(gene_rows, r)
    }
  }

  ## 基因名纵向避让（黑色，超出末端整体上移）
  if (length(gene_rows) > 0) {
    ys <- (blc$start[gene_rows] + blc$end[gene_rows]) / 2 / 1e6
    o <- order(ys); ys_s <- ys[o]; rows_s <- gene_rows[o]
    ylab <- ys_s
    for (q in 2:length(ylab)) if (ylab[q] - ylab[q-1] < GAP) ylab[q] <- ylab[q-1] + GAP
    ov <- ylab[length(ylab)] - (len + 5)
    if (ov > 0) ylab <- ylab - ov
    for (q in seq_along(rows_s)) {
      if (abs(ylab[q] - ys_s[q]) > 1.5)
        segments(x0 + X_GENE + 0.012, ys_s[q], x0 + X_TEXT - 0.015, ylab[q],
                 col = "grey70", lwd = 0.3)
      text(x0 + X_TEXT, ylab[q], blc$known_gene[rows_s[q]], cex = CEX_TXT,
           col = gene_col, adj = c(0, 0.5))
    }
  }

  ## 染色体名
  text(x0 + 0.12, len + 14, chr, cex = 0.85, font = 2, adj = c(0, 0))
}
}

pdf("chromosome_map.pdf", width = 32, height = 12.75)
draw_map()
dev.off()
png("chromosome_map.png", width = 32, height = 12.75, res = 120, units = "in")
draw_map()
dev.off()
