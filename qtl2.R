library(qtl)
library(tidyverse)
setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")

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

data2<-NULL
shuangqin_phe<-read.table("shuangqin_phe2.txt", sep = "\t", header = TRUE)  

for (i in 1:25) {
  genotype<-read.table( paste0(path[i],name1[i],"_genotype_width.txt"), sep = "\t", header = TRUE)  
  rownames(genotype)<-paste0("chr_",genotype$chr,"_",round(genotype$region))
  genotype$region<-genotype$region/1000000
  genotype2<-t(genotype)%>%as.data.frame()%>%
    rownames_to_column(.,"ID")
  genotype2$ID<-gsub("\\.","-",genotype2$ID)
  genotype2$ID<-gsub("_","-",genotype2$ID)
  genotype2$ID<-gsub("^X(\\d+.*)$", "\\1",genotype2$ID)
  phe<-filter(shuangqin_phe,group==name1[i])
  
  phe<-phe[,colSums(is.na(phe)) < (0.25 * nrow(phe))]
  a<-ncol(phe)-2
  phe<-merge(phe,genotype2,by="ID",all = T)
  phe<-phe[rowSums(is.na(phe)) < (0.25 * ncol(phe)), ]  
  
  phe<-t(phe)%>%as.data.frame()
  phe<-rownames_to_column(phe,"ID")
  
  
  phe<-phe[3:nrow(phe),c(1,(ncol(phe)-1),ncol(phe),2:(ncol(phe)-2))]
  phe[1:a,2:3]<-""
  write.table(phe, paste0(name1[i],"_qtl_blue.csv"), sep = ",",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  

  
  ffdata<-read.cross("csvr",".",paste0(name1[i],"_qtl_blue.csv"),crosstype = "f2",genotypes = c("0","1","2"),alleles=c("0","2"))
  
  ffdata<-calc.genoprob(ffdata,step=0,error.prob=0.001)
  p1.hk<-scanone(ffdata,method="hk",pheno.col=1:a)
  
  operm<-scanone(ffdata,pheno.col=1:a,method="hk",n.perm=50)
  thr<-summary(operm,alpha=c(0.01,0.05,0.1))%>%as.matrix()
  thr1<-rbind(thr)%>%as.data.frame()%>%
    rownames_to_column(.,"ID")
  
  write.table(thr1,paste0(name1[i],"qtl_thr.txt") ,sep = '\t',
              row.names = FALSE,col.names = T, quote = FALSE)
  p1.hk<-p1.hk%>%rownames_to_column(.,"ID")
  write.table(p1.hk,paste0(name1[i],"qtl_lod.txt") ,sep = '\t',
              row.names = FALSE, col.names = T, quote = FALSE)
  
  
  for (j in 1:a) {
    pdf(paste0(name1[i],"_p",j,"_qtl.pdf"), width=10, height=4)  
    plot(p1.hk[,c(2,3,j+3)],main=paste0(name1[i],"_p",j,"_qtl"))
    abline(h=thr1[3,j+1], col="black", lwd=2, lty=2)
    dev.off()  
    
    #p1.hk<-read.table("JMqtl_lod.txt", sep = "\t", header = T)
    data1<-p1.hk[,c(1:3,j+3)]
    colnames(data1)<-c("marker","chr","pos","value")
    data1<-data1[data1$value>=thr1[3,j+1],]
    if(nrow(data1)>0){
      data1$phe<-paste0("p",j)
      data1$pop<-name1[i]
      data1$thr<-thr1[3,j+1]
      
      data2<-rbind(data2,data1)
    }
    
  }
}

write.table(data2,"all_qtl.txt",sep = '\t',
            row.names = FALSE, col.names = T, quote = FALSE)
