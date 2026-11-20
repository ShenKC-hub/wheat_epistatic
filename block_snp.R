#!/usr/bin/env Rscript
#setwd("/home/workspace/results/merge/qinben/split_vcf/second/1A_block")

args=commandArgs(TRUE)
library(tidyverse)
library(readr)
intput <- args[1]


path<-c("/home/workspace/results/merge/AVS/adj/",
        
        "/home/workspace/results/merge/BBR/second/adj/",
        
        "/home/workspace/results/merge/EZX/second/adj/",
        
        "/home/workspace/results/merge/GYGZ/adj/",
        
        "/home/workspace/results/merge/GZX/new/adj/",
        
        "/home/workspace/results/merge/JA/second/adj/",
        
        "/home/workspace/results/merge/JM/adj/",
        
        "/home/workspace/results/merge/KJH/second/adj/",
        
        "/home/workspace/results/merge/SKH/second/adj/",
        
        "/home/workspace/results/merge/SZS/second/adj/",
        
        "/home/workspace/results/merge/XNG/second/adj/",
        
        "/home/workspace/results/merge/ZJ/second/adj/",
        
        "/home/workspace/results/merge/BLB/second/adj/",
        
        "/home/workspace/results/merge/DH/second/adj/",
        
        "/home/workspace/results/merge/GCS/second/adj/",
        
        "/home/workspace/results/merge/GZH/adj/",
        
        "/home/workspace/results/merge/GZY/second/adj/",
        
        "/home/workspace/results/merge/JL/second/adj/",
        
        "/home/workspace/results/merge/JZ/second/adj/",
        
        "/home/workspace/results/merge/SZM/second/adj/",
        
        "/home/workspace/results/merge/SZY/second/adj/",
        
        "/home/workspace/results/merge/XLJ/second/adj/",
        
        "/home/workspace/results/merge/ZB/second/adj/",
        
        "/home/workspace/results/merge/ZK/second/adj/",
        
        "/home/workspace/results/merge/BBC/second/adj/")

name1<-c("AVS",  "BBR",  "EZX",  "GYGZ",  "GZX", "JA", "JM",  "KJH", "SKH",  "SZS", "XNG",  "ZJ",
         "BLB",  "DH", "GCS", "GZH", "GZY", "JL", "JZ", "SZM", "SZY",  "XLJ", "ZB", "ZK", "BBC")

data5<-read.table(intput, sep = "\t", header = T)

all_block_ID<- read.delim("/home/workspace/results/merge/qinben/split_vcf/second/all_block_ID.txt", header=T)
all_block_ID$BP1<-as.numeric(all_block_ID$BP1)
all_block_ID$BP2<-as.numeric(all_block_ID$BP2)

data6<-data7<-NULL
a1<-0
for (i in 1:nrow(data5)) {
  block1<-filter(all_block_ID,CHR==data5[i,2],BP1>=data5[i,10],BP2<=data5[i,11])
  
  a<-which(name1==data5[i,8])
  if(a!=a1){snp_vcf<-read.table(paste0(path[a],"new_SNP_",name1[a],".vcf"), sep = "\t", header = F)}
  a1=a
  snp_vcf2<-filter(snp_vcf,V1==data5[i,2],V2>=data5[i,10],V2<=data5[i,11])
  
  gwas<-read.table(paste0("/home/workspace/results/merge/qinben/GWAS/output2/chr",data5[i,2],"_",as.numeric(gsub("P","",data5[i,7]))+1,".GEMMA.assoc.txt"), sep = "\t", header = T)
  gwas1<-filter(gwas,rs%in%snp_vcf2$V3)
  
  if(nrow(block1)>0 & nrow(gwas1)>0){
    gwas1<-merge(gwas1,snp_vcf2,by.x="rs",by.y="V3")
    gwas1<-mutate(gwas1,beta=if_else(allele1==V4,-1*beta,beta))
    data5[i,13]<-mean(gwas1$beta)
    expanded_df <- do.call(rbind, lapply(1:nrow(block1), function(x) data5[i,]))  
    block2<-data.frame(expanded_df,block1)
    data6<-rbind(data6,block2)
    
    
    expanded_df <- do.call(rbind, lapply(1:nrow(gwas1), function(x) data5[i,])) 
    block3<-data.frame(expanded_df,gwas1)
    data7<-rbind(data7,block3)
  }
}


write.table(data6,paste0("qtl_block_",intput),sep = '\t',
            row.names = FALSE, col.names = T, quote = FALSE)

write.table(data7,paste0("qtl_snp_",intput),sep = '\t',
            row.names = FALSE, col.names = T, quote = FALSE)
