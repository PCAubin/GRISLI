library(MASS)

args <- commandArgs(trailingOnly = T)
fdata <- args[1]
ftime <- args[2]
direc <- args[3]
tfnum <- as.numeric(args[4])
pnum <- as.numeric(args[5])
cnum <- as.numeric(args[6])
maxite <- as.numeric(args[7])
repnum <- as.numeric(args[8])

# fdata <- 'data2/datamatrix.txt'
# ftime <- 'data2/pseudotime.txt'
# direc <- out
# tfnum <- 100
# pnum <- 4
# cnum <- 373
# maxite <- 50
# repnum <- 2

maxB <- 2.0
minB <- -10.0

# dir.create(direc)
#system(paste("mkdirec", direc, sep=" "))

X <- as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum]
W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
WZ <- matrix(nrow=tfnum, ncol=cnum)

#read pseudo-time and normalize pseudo-time
pseudotime <- as.vector(as.matrix(read.table(ftime, sep="\t")))#[1:cnum,2]
pseudotime <- pseudotime/max(pseudotime)

new_B <- rep(0, pnum)
old_B <- rep(0, pnum)

meanA<-0
for(k in 1:repnum ){
#initialization
RSS <- Inf
for(i in 1:pnum){
	new_B[i] <- runif(1, min=minB, max=maxB)
	old_B[i] <- new_B[i]
}

#function to sample Z
sample_Z <- function(){
	for(i in 1:pnum){
		for(j in 1:cnum){
			Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
		}
	}
}

#optimize W and B iteratively
for(ite in 1:maxite){
	#sampling B
	target <- floor(runif(1, min=1, max=pnum+1))
	new_B[target] <- runif(1, min=minB, max=maxB)

	#for last calculation
	if(ite == maxite){
		for(i in 1:pnum){
			new_B[i] <- old_B[i]
		}
	}

	#sample Z from new B
 	sample_Z()

	#regression
	for(i in 1:tfnum){
		X.lm <- lm(X[i,] ~ t(Z)-1)
		for(j in 1:pnum){
	    	W[i,j] <- X.lm$coefficients[j]
	  	}
		WZ[i,] <- W[i,] %*% Z
	}

	#RSS
	tmp_RSS <- sum((X-WZ)**2)
	if(tmp_RSS < RSS){
		RSS <- tmp_RSS
	}
	else{
		new_B[target] <- old_B[target]
	}
}

#output RSS
#write.table(RSS, paste(direc,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")

#output W
#write.table(W, paste(direc,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")

#infer A
B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
for(i in 1:pnum){
	B[i,i] <- new_B[i]
}
invW <- ginv(W)
A <- W %*% B %*% invW

meanA <- meanA +A
}
meanA<-meanA/repnum
write.table(meanA,"A.txt", row.names=F, col.names=F, sep="\t")
# prmatrix(meanA, rowlab=rep("",dim(meanA)[1]), collab=rep("",dim(meanA)[2]))
# dim(meanA)
#output A and B
# write.table(A, paste(direc,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
# write.table(A, paste(direc,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
#write.table(B, paste(direc,"/B.txt",sep=""), row.names=F, col.names=F, sep="\t")
