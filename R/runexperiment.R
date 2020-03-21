# Test performance of SCODE and TIGRESS in different settings
# Jean-Philippe Vert
# 6/21/2018


# Ensure reproducibility
set.seed(4395)

# Load packages and functions
source("myscode.R")
library(ROCR)
library(tigress)
library(zinbwave)

# SCODE parameters
repnum=100
maxite=100
# TIGRESS parameters
nstepsLARS<-100
nsplit<-1000
# ZinbWave parameters
zinbK <- 10

# Main loop: two datasets
for (dataset in c(2,3)) {

  ## Data preparation

  # Load data
  datapath <- paste("../SCODE-master/data",dataset,"/",sep="")
  fdata <- paste(datapath,"data.txt",sep="") # file of expression matrix
  ftime <- paste(datapath,"time.txt",sep="") # file of pseudotimes
  fA <- paste(datapath,"A.txt",sep="") # file of gold standard network
  X <- as.matrix(read.table(fdata, sep="\t")) # expression matrix. Rows=genes, col=experiments
  rownames(X) = seq(dim(X)[1]) # give row names (necessary for TIGRESS)
  ntf <- dim(X)[1] # Number of genes (total)
  pseudotime <- read.table(ftime, sep="\t")[,2] # vector of pseudotimes
  trueA <- as.matrix(read.table(fA, sep="\t")) # Gold standard network. Col=regulator, Row=target
  diag(trueA) <- 0 # remove self-regulation
  diagind <- seq(1,ntf*ntf,ntf+1) # indices of the non-diagonal elements in trueA

  # Create smaller dataset with only TF having at least one edge
  keepTF <- apply(trueA,2,sum) + apply(trueA,1,sum) >0 # TF with at least one regulation
  Xsmall <- X[keepTF,]
  ntfsmall <- sum(keepTF)
  trueAsmall <- trueA[keepTF,keepTF]
  diagindsmall <- seq(1,ntfsmall*ntfsmall,ntfsmall+1)

  # Replace 0 by ZinbWave estimates in the expression matrix
  Z <- zinbwave(SummarizedExperiment(round(exp(X)-1)),imputedValues=TRUE,K=zinbK)
  Xzinb <- log(Z@assays[[3]]+1)
  write.table(Xzinb, file=paste("data",dataset,"zinb.txt",sep=""), sep='\t',col.names=FALSE,row.names=FALSE)

  ## Test different algorithms on different datasets
  auc=NULL # a list to store results

  # 1. Train SCODE on all data
  Ascode <- myscode(X,pseudotime,maxite=maxite,repnum=repnum)

  # Perf on all TFs
  pred <- prediction(abs(Ascode[-diagind]), trueA[-diagind])
  auc <- rbind(auc,c("SCODE train all test all",performance(pred, measure = "auc")@y.values))

  # Perf on small set only
  Ascodesmall <- Ascode[keepTF,keepTF]
  pred <- prediction(abs(Ascodesmall[-diagindsmall]), trueAsmall[-diagindsmall])
  auc <- rbind(auc,c("SCODE train all test small",performance(pred, measure = "auc")@y.values))

  # 2. Train SCODE on selected TF only
  Ascodesmall <- myscode(Xsmall,pseudotime,maxite=maxite,repnum=repnum)

  # Perf (on small dataset only)
  pred <- prediction(abs(Ascodesmall[-diagindsmall]), trueAsmall[-diagindsmall])
  auc <- rbind(auc,c("SCODE train small test small",performance(pred, measure = "auc")@y.values))

  # 3. Train SCODE on all data ZINB
  Ascode <- myscode(Xzinb,pseudotime,maxite=maxite,repnum=repnum)

  # Perf on all TFs
  pred <- prediction(abs(Ascode[-diagind]), trueA[-diagind])
  auc <- rbind(auc,c("SCODE ZINB train all test all",performance(pred, measure = "auc")@y.values))

  # Perf on small set only
  Ascodesmall <- Ascode[keepTF,keepTF]
  pred <- prediction(abs(Ascodesmall[-diagindsmall]), trueAsmall[-diagindsmall])
  auc <- rbind(auc,c("SCODE ZINB train all test small",performance(pred, measure = "auc")@y.values))

  # 4. Train SCODE ZINB on selected TF only
  Ascodesmall <- myscode(Xzinb[keepTF,],pseudotime,maxite=maxite,repnum=repnum)

  # Perf (on small dataset only)
  pred <- prediction(abs(Ascodesmall[-diagindsmall]), trueAsmall[-diagindsmall])
  auc <- rbind(auc,c("SCODE ZINB train small test small",performance(pred, measure = "auc")@y.values))

  # 5. Train TIGRESS on all TF
  predTigress <- tigress(t(X),nstepsLARS=nstepsLARS,nsplit=nsplit)

  # Perf on all TF
  aucTigress <- lapply(predTigress,function(pr) {
    pred <- prediction(pr[-diagind], t(trueA)[-diagind])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigress))) {
    auc <- rbind(auc,c(paste("TIGRESS train all test all L=",i,sep=""),aucTigress[[i]][[1]][1]))
  }
  # Perf on small TF set
  aucTigressSmall <- lapply(predTigress,function(pr) {
    pred <- prediction(pr[keepTF,keepTF][-diagindsmall], t(trueAsmall)[-diagindsmall])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigressSmall))) {
    auc <- rbind(auc,c(paste("TIGRESS train all test small L=",i,sep=""),aucTigressSmall[[i]][[1]][1]))
  }

  # 6. Train TIGRESS on small data
  predTigressSmall <- tigress(t(X[keepTF,]),nstepsLARS=nstepsLARS,nsplit=nsplit)
  aucTigressSmallSmall <- lapply(predTigressSmall,function(pr) {
    pred <- prediction(pr[-diagindsmall], t(trueAsmall)[-diagindsmall])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigressSmallSmall))) {
    auc <- rbind(auc,c(paste("TIGRESS train small test small L=",i,sep=""),aucTigressSmallSmall[[i]][[1]][1]))
  }

  # 7. Train TIGRESS on all TF
  predTigress <- tigress(t(Xzinb),nstepsLARS=nstepsLARS,nsplit=nsplit)

  # Perf on all TF
  aucTigress <- lapply(predTigress,function(pr) {
    pred <- prediction(pr[-diagind], t(trueA)[-diagind])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigress))) {
    auc <- rbind(auc,c(paste("TIGRESS ZINB train all test all L=",i,sep=""),aucTigress[[i]][[1]][1]))
  }
  # Perf on small TF set
  aucTigressSmall <- lapply(predTigress,function(pr) {
    pred <- prediction(pr[keepTF,keepTF][-diagindsmall], t(trueAsmall)[-diagindsmall])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigressSmall))) {
    auc <- rbind(auc,c(paste("TIGRESS ZINB train all test small L=",i,sep=""),aucTigressSmall[[i]][[1]][1]))
  }

  # 8. Train TIGRESS ZINB on small data
  predTigressSmall <- tigress(t(Xzinb[keepTF,]),nstepsLARS=nstepsLARS,nsplit=nsplit)
  aucTigressSmallSmall <- lapply(predTigressSmall,function(pr) {
    pred <- prediction(pr[-diagindsmall], t(trueAsmall)[-diagindsmall])
    performance(pred, measure = "auc")@y.values
  })
  for (i in seq(length(aucTigressSmallSmall))) {
    auc <- rbind(auc,c(paste("TIGRESS ZINB train small test small L=",i,sep=""),aucTigressSmallSmall[[i]][[1]][1]))
  }

  write.table(auc, file=paste("auc",dataset,".txt",sep=""), sep='\t',col.names=FALSE,row.names=FALSE)
  save.image(file=paste("data",dataset,".Rdata",sep=""))
}
