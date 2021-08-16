# .onUnload <- function (libpath) {
#   library.dynam.unload("glasso", libpath)
# }
#
# .onUnload <- function (libpath) {
#   library.dynam.unload("SFGen", libpath)
# }

norm_vec_inf <- function(X){ return(max(abs(X))) }
tr_mat <- function(X){  return(sum(diag(X)))  }
norm_vec <- function(x) sqrt(sum(x^2))

Gene_cov<-function(p){
  sigma <- runif(p-1,0.5,1)
  covmat0 <- diag(1,p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      temp <- exp(-sum(sigma[i:(j-1)]/2))
      covmat0[i,j] <- temp
      covmat0[j,i] <- temp
    }
  }
  return(covmat0)
}

rmn <- function(M, Srow, Scol){
  m=dim(Srow)[1]
  n=dim(Scol)[1]
  tmp=eigen(Srow)
  Srow.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
  tmp=eigen(Scol)
  Scol.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
  Z=matrix(rnorm(m * n), m, n)
  Srow.h %*% Z %*% Scol.h + M
}

rmn_rep <- function(N, M, Srow, Scol){
  m=dim(Srow)[1]
  n=dim(Scol)[1]
  tmp=eigen(Srow)
  Srow.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
  tmp=eigen(Scol)
  Scol.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
  Z=matrix(rnorm(m*n*N), m, n*N)
  for (i in 1:N){
    Z[,((i-1)*n+1):(i*n)] <- Srow.h %*% Z[,((i-1)*n+1):(i*n)] %*% Scol.h + M
  }
  Z
}


cor.est <- function(X,N){
  m <- dim(X)[1]
  n <- dim(X)[2] / N
  S.hat.row <- matrix(0,m,m)
  S.hat.col <- matrix(0,n,n)
  for (i in 1:N){
    S.hat.row <- S.hat.row + tcrossprod(X[,((i-1)*n+1):(i*n)])
    S.hat.col <- S.hat.col + crossprod(X[,((i-1)*n+1):(i*n)])
  }
  S.hat.row <- S.hat.row / N
  S.hat.col <- S.hat.col / N
  tmp <- diag(1/sqrt(diag(S.hat.row)))
  cor.row <- tmp %*% S.hat.row %*% tmp
  tmp <- diag(1/sqrt(diag(S.hat.col)))
  cor.col <- tmp %*% S.hat.col %*% tmp
  list(cor.row = cor.row, cor.col = cor.col,W.row = diag(S.hat.row),W.col = diag(S.hat.col))
}



soft.threshold.off.diag <- function(mat,lam){
  tmp.diag <- diag(mat)
  diag(mat) <- 0
  mat[abs(mat)<lam] <- 0
  idx.1 <- which(mat < -lam)
  idx.2 <- which(mat > lam)
  if ( length(idx.1)>0 ) mat[idx.1] <- mat[idx.1]+lam
  if ( length(idx.2)>0 ) mat[idx.2] <- mat[idx.2]-lam
  diag(mat) <- tmp.diag
  return ( mat )
}

soft.threshold <- function(vec,lam){
  if ( length(lam)>1 & length(lam)!=length(vec) ) {
    cat('\n ERROR: THE SIZE OF THE SECOND ARGUMENT SHOULD BE 1 OR THE SAME AS THE SIZE OF THE FIRST ARGUMENT.\n')
    return ( 0 )
  }
  vec[abs(vec)<lam] <- 0
  idx.1 <- which(vec < -lam)
  idx.2 <- which(vec > lam)
  if (length(lam)==1){
    if ( length(idx.1)>0 ) vec[idx.1]<- vec[idx.1]+lam
    if ( length(idx.2)>0 ) vec[idx.2]<- vec[idx.2]-lam
  } else if (length(lam)>1) {
    if ( length(idx.1)>0 ) vec[idx.1]<- vec[idx.1]+lam[idx.1]
    if ( length(idx.2)>0 ) vec[idx.2]<- vec[idx.2]-lam[idx.2]
  }
  return( vec )
}

max.singular.value <- function(A){
  norm_vec <- function(x) sum(abs(x))
  row_norm <- as.matrix(apply(A,1,norm_vec))
  col_norm <- as.matrix(apply(A,2,norm_vec))
  tmp1 <- abs(t(A)) %*% row_norm
  tmp2 <- abs(A) %*% col_norm
  max(c(tmp1,tmp2))
}

report.roc <- function(omega.hat, omega0){
  diag(omega0) = 0
  diag(omega.hat) = 0
  nonzero_idx = which(abs(omega0) > 1e-10)
  zero_idx = which(abs(omega0) < 1e-10)
  T1 = T2 = T1_rate = T2_rate = NULL
  if (length(nonzero_idx) == 0) { T1 = 0; T1_rate = 0; } else {
    T1 = length(which(abs(omega.hat[nonzero_idx]) < 1e-10))
    T1_rate = T1 / length(nonzero_idx)
  }
  if (length(zero_idx) == 0) { T2 = 0; T2_rate = 0; } else {
    T2 = length(which(abs(omega.hat[zero_idx]) > 1e-10))
    T2_rate = T2 / length(zero_idx)
  }
  return (list(T1=T1,T2=T2,T1_rate=T1_rate,T2_rate=T2_rate))
}


