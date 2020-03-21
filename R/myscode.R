#' Re-implementation of the SCODE method as an R function.
#'
#' Original: https://github.com/hmatsu1226/SCODE/blob/master/SCODE.R
#'
#'
#' @param X a G x C matrix of expression data
#' @param pseudotime a vector of length C of pseudo-times
#' @param pnum the number of z
#' @param maxite the number of iterations for convergence
myscode <-
  function(X, pseudotime, targetlist=colnames(expdata), pnum=4, maxite=50, repnum=1, maxB=2.0, minB=-10.0)
  {
    require(MASS)
    tfnum <- dim(X)[1]
    cnum <- dim(X)[2]
    pseudotime <- pseudotime/max(pseudotime)

    meanA <- matrix(rep(0,tfnum*tfnum), nrow=tfnum, ncol=tfnum)
    for(irep in 1:repnum ){

    W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
    Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
    WZ <- matrix(nrow=tfnum, ncol=cnum)

    new_B <- rep(0, pnum)
    old_B <- rep(0, pnum)

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
    #infer A
    B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
    for(i in 1:pnum){
      B[i,i] <- new_B[i]
    }
    invW <- ginv(W)
    meanA <- meanA + W %*% B %*% invW
    }
    meanA/repnum
  }
