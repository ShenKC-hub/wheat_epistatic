setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(MatrixEpistasis)
name<-c("AVS","BBC","BBR","BLB","DH","EZX","GCS","GYGZ","GZH","GZX","GZY",
"JA","JL","JM","JZ","KJH","SKH","SZM","SZS","SZY","XLJ","XNG","ZB","ZJ","ZK")

for (j in 1:length(name)) {
  df <- read.csv(paste0(name[j],"_qtl_blue.csv"),header = F)
df$V3<-df$V3*1000000

all_block <- read.delim("all_block.bed", header=FALSE, comment.char="#")

geno<-df[!(is.na(df$V3))  ,]
df_odd <- geno[seq(1, nrow(geno), by = 2), ]
df_even <- geno[ (1:nrow(geno)) %% 2 == 0, ]

df_odd <- df_odd %>%
  mutate(end = df_even$V3, .after = 3) 
data<-NULL
for (i in 1:nrow(all_block)) {
  df_odd2<-filter(df_odd,V2==all_block$V1[i],V3<=all_block$V2[i],end>=all_block$V3[i])
  df_odd3<-filter(df_odd,V2==all_block$V1[i],V3>all_block$V2[i],V3<all_block$V3[i])
  df_odd4<-rbind(df_odd2,df_odd3)
  data<-rbind(data,df_odd4[1,])
  }

for (i in 1:7) {
  trait<-t(df[i,4:ncol(df)])
  snpA = t(data[,5:ncol(data)])
  colnames(snpA)<-all_block$V4
  snpA<-cbind(trait,snpA)
  snpA<-na.omit(snpA)
  trait<-snpA[,1]
  snpA<-snpB<-snpA[,-1]

 # covariate <- prcomp(snpA)$x[,1:5]
  MatrixEpistasis_main( snpA=snpA, snpB=snpB, trait=trait, #covariate=covariate ,
                      pvalThreshold=1e-5, outputFileName=paste0("Epistasis_",name[j],"_P", i,"_txt"))
}
}

