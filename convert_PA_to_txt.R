setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
suppressMessages(library(data.table))
for (i in 1:7) {
  t0 <- Sys.time()
  PA <- readRDS(paste0("gwas_epistasis_PA_P", i, ".rds"))
  fwrite(as.data.frame(as.matrix(PA)), paste0("gwas_epistasis_PA_P", i, "_txt"),
         sep = "\t", quote = FALSE, row.names = TRUE)
  cat("P", i, "done:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min\n")
}
cat("ALL DONE\n")
