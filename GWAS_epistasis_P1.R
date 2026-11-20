setwd("/home/workspace/results/merge/BBC/second/adj/new_qtl")
library(MatrixEpistasis)
library(tidyverse)
genotype <- read.delim("gwas_all_block_genotype.txt", header=FALSE, comment.char="#")
genotype <-genotype[,-c(1:2,4:9)]

genotype <-t(genotype)

phe <- read.delim("phe2.txt", header=FALSE, comment.char="#")

pca <- read.delim("pca.txt", header=FALSE, comment.char="#",sep = "\t")

phe<-merge(phe,pca,by="V1")

colnames(genotype)<-genotype[1,]
genotype<-genotype[-1,]

phe<-merge(phe,genotype,by.x="V1",by.y="ID")

  i=1
  snpA<-phe[,c(i+1,9:ncol(phe))]
  snpA<-na.omit(snpA)
  trait<-snpA[,1]
  covariate <- snpA[,2:4]
  
  snpA <- snpA[,-c(1:4)]
  snpA[] <- lapply(snpA, function(x) as.numeric(x))
  snpB<-snpA<-as.matrix(snpA)
  
  
  MatrixEpistasis_main( snpA=snpA, snpB=snpB, trait=trait, covariate=covariate ,
                        pvalThreshold=1e-5, outputFileName=paste0("gwas_Epistasis_P",i,"_txt"))
