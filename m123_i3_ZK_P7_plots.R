# m123_i3_ZK_P7_plots.R
# 以 ZK 群体 P7 性状为例, 展示 m123_i3 模型(3主效应+3两两交互+1三元交互)
#   图1: 拟合值 vs 实际值 散点图 (整体关联强度, 标准化量纲)
#   图2: 两两交互切片图 (边数/模块数/完整度 两两组合, 标准化表型)
#   图3: 预测景观热图 (标准化坐标, 按模块数 低/均值/高分面 -> 展示三元交互)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

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

## 标准化 + 交互项 (与 parent_complementarity_h2_gain_v2.R 完全一致)
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
fitted_v <- fitted(fit)
r_fit <- cor(fitted_v, dz$y)
cat(sprintf("m123_i3: R2=%.4f  P=%.3g  fitted-observed r=%.3f\n", R2, mp, r_fit))

## ============ 图1: 拟合值 vs 实际值 (标准化) ============
p1 <- ggplot(NULL, aes(x=fitted_v, y=dz$y)) +
  geom_point(shape=21, fill="#4DAF4A", colour="grey30", alpha=0.6,
             size=2.2, stroke=0.3) +
  geom_smooth(method="lm", se=TRUE, colour="#E41A1C", fill="#E41A1C",
              alpha=0.12, linewidth=0.6) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey40") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.6, size=4.2,
           fontface="italic",
           label=sprintf("r = %.3f\nR2 = %.3f,  model P = %.1e\nn = %d",
                         r_fit, R2, mp, nrow(d))) +
  labs(x=sprintf("Predicted %s (standardized, m123_i3 model)", TR),
       y=sprintf("%s phenotype (standardized)", TR),
       title=sprintf("%s population, %s trait: m123_i3 model\n3 main effects + 3 two-way interactions + 1 three-way interaction",
                     POP, TR)) +
  theme_bw()
ggsave("m123_i3_ZK_P7_fitted_observed.pdf", p1, width=6.5, height=6.5)

## ============ 图2: 两两交互切片图 ============
## x = 变量A(原始值), y = P7(标准化), 按变量B三分位分组画回归线
mk_slice <- function(vx, vg, labx, labg) {
  dd <- data.frame(x=d[[vx]], y=dz$y,
                   g=cut(d[[vg]], breaks=quantile(d[[vg]], probs=c(0,1/3,2/3,1),
                                                  na.rm=TRUE),
                         include.lowest=TRUE, labels=c("Low","Mid","High")))
  dd <- dd[complete.cases(dd), ]
  nn <- table(dd$g)
  dd$g <- factor(as.character(dd$g),
                 levels=c("Low","Mid","High"),
                 labels=sprintf("%s (n=%d)", c("Low","Mid","High"), nn[c("Low","Mid","High")]))
  ggplot(dd, aes(x=x, y=y, colour=g)) +
    geom_point(alpha=0.35, size=1.4, show.legend=FALSE) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.8) +
    labs(x=labx, y=sprintf("%s phenotype (standardized)", TR),
         colour=labg, subtitle=sprintf("Interaction: %s x %s", vx, vg)) +
    theme_bw() + theme(legend.position="bottom", legend.title=element_blank(),
                       legend.text=element_text(size=8))
}
s1 <- mk_slice("mean_completeness_effective","n_edges_carried",
               LAB["mean_completeness_effective"], "Edges")
s2 <- mk_slice("mean_completeness_effective","n_modules_covered",
               LAB["mean_completeness_effective"], "Modules")
s3 <- mk_slice("n_modules_covered","n_edges_carried",
               LAB["n_modules_covered"], "Edges")
if (requireNamespace("patchwork", quietly=TRUE)) {
  ggsave("m123_i3_ZK_P7_interaction_slices.pdf",
         patchwork::wrap_plots(s1, s2, s3, nrow=1), width=13, height=4.8)
} else {
  library(cowplot)
  ggsave("m123_i3_ZK_P7_interaction_slices.pdf",
         plot_grid(s1, s2, s3, nrow=1), width=13, height=4.8)
}

## ============ 图3: 预测景观热图 ============
## 网格: 边数 x 完整度(标准化), 模块数固定在 低(-1SD)/中(均值)/高(+1SD)
e_seq <- seq(min(dz$n_edges_carried), max(dz$n_edges_carried), length.out=60)
c_seq <- seq(min(dz$mean_completeness_effective),
             max(dz$mean_completeness_effective), length.out=60)
m_lv <- c("Low (-1SD)"=-1, "Mid (mean)"=0, "High (+1SD)"=1)
grid <- expand.grid(n_edges_carried=e_seq,
                    mean_completeness_effective=c_seq,
                    m_val=unname(m_lv))
grid$n_modules_covered <- grid$m_val
grid$int_em_mo  <- grid$n_edges_carried * grid$n_modules_covered
grid$int_em_cpe <- grid$n_edges_carried * grid$mean_completeness_effective
grid$int_mo_cpe <- grid$n_modules_covered * grid$mean_completeness_effective
grid$int_all    <- grid$n_edges_carried * grid$n_modules_covered * grid$mean_completeness_effective
grid$pred <- predict(fit, newdata=grid)
grid$mod_lab <- factor(names(m_lv)[match(grid$m_val, unname(m_lv))],
                       levels=names(m_lv))

p3 <- ggplot(grid, aes(x=n_edges_carried, y=mean_completeness_effective,
                       fill=pred)) +
  geom_tile() +
  geom_contour(aes(x=n_edges_carried, y=mean_completeness_effective, z=pred),
               inherit.aes=FALSE, colour="grey30", linewidth=0.3, bins=8) +
  geom_point(data=data.frame(dz, check.names=FALSE),
             aes(x=n_edges_carried, y=mean_completeness_effective, fill=NULL),
             inherit.aes=FALSE, alpha=0.25, size=1.1, colour="grey20") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, name=sprintf("Predicted %s\n(standardized)", TR)) +
  facet_wrap(~mod_lab) +
  scale_x_continuous(breaks=-2:2) + scale_y_continuous(breaks=-2:2) +
  labs(x=LAB["n_edges_carried"], y=LAB["mean_completeness_effective"],
       title=sprintf("%s population, %s trait: predicted phenotype landscape (m123_i3 model)\nmodules fixed at low/mean/high (-1SD/mean/+1SD); points = observed progeny",
                     POP, TR)) +
  theme_bw() + theme(strip.text=element_text(face="bold"))
ggsave("m123_i3_ZK_P7_prediction_landscape.pdf", p3, width=11, height=4.5)

cat("saved: m123_i3_ZK_P7_fitted_observed.pdf\n",
    "      m123_i3_ZK_P7_interaction_slices.pdf\n",
    "      m123_i3_ZK_P7_prediction_landscape.pdf\nDONE\n")