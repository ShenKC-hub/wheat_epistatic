# n_favorable_single_R2.R
# n_favorable 单个指标对表型的变异解释度: 每群体×性状做单变量线性回归
# 表型 ~ n_favorable, R2 = Pearson r 的平方, 按用户约定用 Pearson
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages({library(ggplot2); library(dplyr)})

pops <- read.delim("qinben_group2.txt", header=FALSE,
                   col.names=c("pop","p0","p2"), comment.char="#")
tr7 <- paste0("P", 1:7)
sp <- read.delim("shuangqin_phe2.txt", header=TRUE, check.names=FALSE)
names(sp)[1:2] <- c("group","ID")
norm_id <- function(x) gsub("[_.]", "-", sub("^X(?=[0-9])", "", x, perl=TRUE))

res <- NULL
for (r in seq_len(nrow(pops))) {
  pop <- pops$pop[r]
  for (i in 1:7) {
    f <- file.path("qinben_progeny_favorable_hap", pop,
                   paste0(pop, "_progeny_favorable_P", i, "_txt"))
    if (!file.exists(f)) next
    d <- read.delim(f)
    d$ID <- norm_id(d$progeny)
    d <- merge(d, sp[sp$group==pop, c("ID", tr7[i])], by="ID")
    d <- d[, c("n_favorable", tr7[i])]
    d <- d[complete.cases(d), ]
    if (nrow(d) < 30) next
    ct <- cor.test(d$n_favorable, d[[tr7[i]]], method="pearson")
    res <- rbind(res, data.frame(
      pop=pop, trait=tr7[i], n=nrow(d),
      r=round(ct$estimate,3), p=signif(ct$p.value,3),
      R2=round(ct$estimate^2, 4)))
  }
}
res <- res[order(res$pop, res$trait), ]
write.table(res, "n_favorable_single_R2_txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("saved: n_favorable_single_R2_txt\n")

summ <- res %>% group_by(trait) %>%
  summarise(n_models=n(), n_sig=sum(p<0.05),
            mean_R2=round(mean(R2),4), median_R2=round(median(R2),4),
            min_R2=round(min(R2),4), max_R2=round(max(R2),4),
            .groups="drop") %>% as.data.frame()
write.table(summ, "n_favorable_single_R2_summary_txt",
            sep="\t", quote=FALSE, row.names=FALSE)
print(summ, right=FALSE)

p <- ggplot(res, aes(x=trait, y=R2)) +
  geom_boxplot(fill="#8DD3C7", outlier.size=0.9, width=0.6) +
  geom_jitter(width=0.12, size=0.5, alpha=0.4, colour="grey25") +
  labs(title="Variance in phenotype explained by n_favorable alone (single-metric R2)",
       subtitle="Each point = one population; R2 from Pearson correlation of n_favorable with phenotype",
       x="Trait", y="Single-metric R2") +
  theme_bw()
ggsave("n_favorable_single_R2_boxplot.pdf", p, width=8, height=6)
cat("saved: n_favorable_single_R2_boxplot.pdf\nDONE\n")