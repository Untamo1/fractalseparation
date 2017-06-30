
frobnorm <- function(X){
  # crossprod(Conj(Z.matrix),Z.matrix) = Conj(t(Z.matrix)) %*% Z.matrix
  tmp <- sqrt(sum(diag(tcrossprod(X,Conj(X)))))
  #I like to add these parts for myself just to make sure that there is not too much computational error
  if(sum(Im(tmp))>1e-14){ 
    print("Error in frobnorm")
  }
  return(Re(tmp))
}


cdet <- function(X){
  return(prod(eigen(X, only.values=TRUE)$values))
}

pMatrix.min <- function(A){
  cost <- t(apply(A^2, 1, sum) - 2 * A + 1)
  vec <- c(solve_LSAP(cost))
  list(A=A[vec,], pvec=vec)
}


MD_fun <- function(W.hat,A){
  G <- W.hat %*% A
  RowNorms <- sqrt(rowSums(abs(G)^2))
  G.0 <- sweep(abs(G),1,RowNorms, "/")
  G.tilde <- G.0^2
  p <- nrow(A)
  Pmin <- pMatrix.min(G.tilde)
  G.tilde.p <- Pmin$A
  md <- sqrt(p - sum(diag(G.tilde.p)))/sqrt(p-1)
  return(md)
}


msqrt <- function(X){ #matrix squareroot
  #Eigendecomposition
  eigen(X)$vectors %*% (diag(eigen(X)$values)^(1/2) %*%Conj(t(eigen(X)$vectors))) 
}

msqrt2 <- function(A){
  tmp <- svd(A)
  result <- tmp$u %*% diag((tmp$d)^(1/2)) %*% Conj(t(tmp$v))
  return(result)
}

cov1 <- function(X){
  cent <- sweep(X,2,colMeans(X),"-")
  return(1/(dim(X)[1]-1)*crossprod(cent,Conj(cent)))
}

#Named cov2 since package JADE has a cov4() function
cov2 <- function(X){ #cov4 matrix
  p <- dim(X)[2]
  n <- dim(X)[1]
  cent <- sweep(X,2,colMeans(X),"-")
  S1 <- cov1(X)   
  SQ1 <- solve(msqrt(S1))
  Z=matrix(data=NA, nrow=n, ncol=p) 
  for(k in 1:n){
    Z[k,] <- SQ1 %*% cent[k,]
  }
  K = matrix(0,nrow=p,ncol=p)
  for(k in 1:n){
    a <- Z[k,] %*% t(Conj(Z[k,])) %*%  Z[k,] %*% t(Conj(Z[k,]))
    K <- K + a
  }
  msqrt(S1) %*% K %*% msqrt(S1) /(n*(p+2))
}

acov <- function(X,tau){ #autocovariance matrix
  n <- dim(X)[1]
  Z <- sweep(X,2,colMeans(X),"-")
  tmp1 <- crossprod(Z[1:(n-tau),],Conj(Z[(1+tau):n,])) 
  AC <- 1/(2*(n-tau)) * (tmp1 + t(Conj(tmp1)))
  return(AC)
}


NAMUSE <- function(D,t){
  #return the gamma matrix estimate
  n <- dim(D)[1]
  p <- dim(D)[2]
  cent <- sweep(D,2,colMeans(D),"-")
  S1 <- cov1(D)
  COV.sqrt.i <- solve(msqrt(S1))
  #Cov is affine equivariant, lose less precision
  # when calculated from original
  Z <- tcrossprod(cent,COV.sqrt.i) #whitening
  S2 <- acov(Z,t)
  U2 <- eigen(S2,symmetric=TRUE)$vectors
  G<- crossprod(Conj(U2),COV.sqrt.i)
  DT <- Z %*% Conj(U2)
  L <- list(Gamma=G,Data=DT)
  return(L)
}


NFOBI <- function(D){
  #return the gamma matrix estimate
  n <- dim(D)[1]
  p <- dim(D)[2]
  cent <- sweep(D,2,colMeans(D),"-")
  S1 <- cov1(D)
  COV.sqrt.i <- solve(msqrt(S1))
  #Cov is affine equivariant, lose less precision
  # when calculated from original
  Z <- tcrossprod(cent,COV.sqrt.i) #whitening
  S2 <- cov2(Z)
  U2 <- eigen(S2,symmetric=TRUE)$vectors
  G<- crossprod(Conj(U2),COV.sqrt.i)
  DT <- Z %*% Conj(U2)
  L <- list(Gamma=G,Data=DT)
  return(L)
}

NSOBI <- function(D){
  #return the gamma matrix estimate
  n <- dim(D)[1]
  p <- dim(D)[2]
  cent <- sweep(D,2,colMeans(D),"-")
  S1 <- cov1(D)
  COV.sqrt.i <- solve(msqrt(S1))
  #Cov is affine equivariant, lose less precision
  # when calculated from original
  Z <- tcrossprod(cent,COV.sqrt.i) #whitening
  S2 <- cov2(Z)
  U2 <- eigen(S2,symmetric=TRUE)$vectors
  G<- crossprod(Conj(U2),COV.sqrt.i)
  DT <- Z %*% Conj(U2)
  L <- list(Gamma=G,Data=DT)
  return(L)
}




