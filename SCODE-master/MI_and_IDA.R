workingDir = "C:/Users/pierr/OneDrive/Documents/Biologie_Vert/Code_bio/Code GRISLI v2/SCODE-master";
setwd(workingDir);
# Script to run IDA (pcalg) and mutual inference-based network inference-based
# the following parameters have to be changed by hand
fdata <- 'data3/datamatrix.txt'
tfnum <- 100 #100
cnum <- 758 #373
X <- t(as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum])
# W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
# Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
# WZ <- matrix(nrow=tfnum, ncol=cnum)

library(infotheo)
Xdisc<-discretize(X)
I <- mutinformation(Xdisc,method= "mm")
write.table(I,"data3/A_MI.txt", row.names=F, col.names=F, sep="\t")

library(pcalg)
suffStat <- list(C = cor(X), n = cnum)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=tfnum, alpha = 0.01)
## Supppose that we know the true CPDAG and covariance matrix
IDA.c1<- idaFast(1,seq(1,tfnum), cov(X),
                    pc.fit@graph)
covX<-cov(X)
PredictedA_min<-matrix(0, tfnum, tfnum)
PredictedA_mean<-matrix(0, tfnum, tfnum)
PredictedA_number<-matrix(0, tfnum, tfnum)
for(i in 1:tfnum){
  temp<-abs(idaFast(i,seq(1,tfnum),covX,
                          pc.fit@graph))
  PredictedA_min[,i]<-apply(temp, 1, FUN=min)
  PredictedA_mean[,i]<-rowMeans(temp)
  PredictedA_number[,i]<-rowSums(array(as.logical(temp),dim(temp)))
}
# PredictedA_mean_stored<-PredictedA_mean
write.table(PredictedA_min,"A_IDA_min.txt", row.names=F, col.names=F, sep="\t")
write.table(PredictedA_mean,"A_IDA_mean.txt", row.names=F, col.names=F, sep="\t")

# # data(gmG)
# ## Simulate the true DAG
# set.seed(123)
# p <- 7
# myDAG <- randomDAG(p, prob = 0.2) ## true DAG
# myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
# ## simulate Gaussian data from the true DAG
# n <- 10000
# dat <- rmvDAG(n, myDAG)
# ## estimate CPDAG and PDAG -- see help(pc)
# suffStat <- list(C = cor(dat), n = n)
# pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01)
# ## Supppose that we know the true CPDAG and covariance matrix
# (l.ida.cpdag <- ida(2,7, cov(dat),
#                     myCPDAG, method = "local", type = "cpdag"))