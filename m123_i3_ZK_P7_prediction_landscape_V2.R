# m123_i3_ZK_P7_prediction_landscape_V2.R
# V2版: 在原图3(标准化坐标)基础上,
#   面板1: 边数 x 完整度   (模块数固定在均值)
#   面板2: 边数 x 模块数   (完整度固定在均值)
#   面板3: 模块数 x 完整度 (边数固定在均值)
# 三面板组合为一图; 标准化坐标与原图一致, 含等高线与观测子代点
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr)})

POP <- "ZK"; TR <- "P7"; METRICS <- c("n_edges_carried","n_modules_covered",
                                      "mean_completeness_effective")
LAB <- c(n_edges_carried="Number of favorable epistasis pairs (edges)",
         n_modules_covered="Number of modules covered",
         mean_completeness_effective="Mean completeness (effective)")
names(LAB) <- METRICS

sp <- read.delim("shuangqin_phe2.txt", header=TRUE, check.names=FALSE)
names(sp)[1:2] <- c("group","ID")
norm_id <- function(x) gsub("[_.]", "-", sub("^X(?=[0-9])", "", x, perl=TRUE))

st <- read.delim(file.path("qinben_progeny_module", POP,
                           paste0(POP, "_progeny_stats_", TR, "_txt")))
st$ID <- norm_id(st$progeny)
d <- merge(st, sp[sp$group==POP, c("ID", TR)], by="ID")
d <- d[, c(METRICS, TR)]
d <- d[complete.cases(d), ]
cat("n =", nrow(d), "\n")

## 标准化 + 交互项 (与 m123_i3_ZK_P7_plots.R 完全一致)
dz <- as.data.frame(scale(d))
names(dz) <- c(METRICS, "y")
dz$int_em_mo  <- dz$n_edges_carried * dz$n_modules_covered
dz$int_em_cpe <- dz$n_edges_carried * dz$mean_completeness_effective
dz$int_mo_cpe <- dz$n_modules_covered * dz$mean_completeness_effective
dz$int_all    <- dz$n_edges_carried * dz$n_modules_covered * dz$mean_completeness_effective
terms7 <- c(METRICS, "int_em_mo","int_em_cpe","int_mo_cpe","int_all")
fit <- lm(as.formula(paste("y ~", paste(terms7, collapse="+"))), data=dz)
s <- summary(fit)
R2 <- s$r.squared; mp <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                            lower.tail=FALSE)
cat(sprintf("m123_i3: R2=%.4f  P=%.3g\n", R2, mp))

## ---- 每对指标的预测景观(第三指标固定在均值=0, 标准化坐标) ----
TRIPLES <- list(c("n_edges_carried","mean_completeness_effective","n_modules_covered"),
                c("n_edges_carried","n_modules_covered","mean_completeness_effective"),
                c("n_modules_covered","mean_completeness_effective","n_edges_carried"))
make_landscape <- function(tv) {
  xv <- tv[1]; yv <- tv[2]; fixedv <- tv[3]
  x_seq <- seq(min(dz[[xv]]), max(dz[[xv]]), length.out=60)
  y_seq <- seq(min(dz[[yv]]), max(dz[[yv]]), length.out=60)
  gr <- data.frame(n_edges_carried=0, n_modules_covered=0,
                   mean_completeness_effective=0)[rep(1, length(x_seq)*length(y_seq)), ]
  gr[[xv]]     <- rep(x_seq, each=length(y_seq))
  gr[[yv]]     <- rep(y_seq, times=length(x_seq))
  gr[[fixedv]] <- 0        # 均值(标准化尺度) -> 其主效应及交互贡献为0
  gr$int_em_mo  <- gr$n_edges_carried * gr$n_modules_covered
  gr$int_em_cpe <- gr$n_edges_carried * gr$mean_completeness_effective
  gr$int_mo_cpe <- gr$n_modules_covered * gr$mean_completeness_effective
  gr$int_all    <- gr$n_edges_carried * gr$n_modules_covered * gr$mean_completeness_effective
  gr$pred <- predict(fit, newdata=gr)
  gr$x <- gr[[xv]]; gr$y <- gr[[yv]]
  gr$pair <- sprintf("%s (x-axis)\n%s (y-axis)\n%s fixed at mean",
                     LAB[xv], LAB[yv], LAB[fixedv])
  gr
}
land <- do.call(rbind, lapply(TRIPLES, make_landscape))
land$pair <- factor(land$pair, levels=unique(land$pair))

pts_stack <- do.call(rbind, lapply(seq_along(TRIPLES), function(k) {
  xv <- TRIPLES[[k]][1]; yv <- TRIPLES[[k]][2]
  data.frame(pair=levels(land$pair)[k], x=dz[[xv]], y=dz[[yv]])
}))
stopifnot(nrow(pts_stack)==nrow(d)*3)

p3 <- ggplot(land, aes(x=x, y=y, fill=pred)) +
  geom_tile() +
  geom_contour(aes(x=x, y=y, z=pred), inherit.aes=FALSE,
               colour="grey30", linewidth=0.3, bins=8) +
  geom_point(data=pts_stack, aes(x=x, y=y), inherit.aes=FALSE,
             alpha=0.25, size=1.1, colour="grey20") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, name=sprintf("Predicted %s\n(standardized)", TR)) +
  facet_wrap(~pair, nrow=1) +
  scale_x_continuous(breaks=-2:2) + scale_y_continuous(breaks=-2:2) +
  labs(x="Indicator (standardized)", y="Indicator (standardized)",
       title=sprintf("%s population, %s trait: predicted phenotype landscape (m123_i3 model)\neach panel: joint effect of a pair of indicators, third indicator fixed at mean; points = observed progeny",
                     POP, TR)) +
  theme_bw() + theme(strip.text=element_text(face="bold", size=8))
ggsave("m123_i3_ZK_P7_prediction_landscape_V2.pdf", p3, width=13, height=4.5)

cat("saved: m123_i3_ZK_P7_prediction_landscape_V2.pdf\nDONE\n")