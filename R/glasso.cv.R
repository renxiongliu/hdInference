loss.like <- function(omega.hat,sigma)
{
  return ( -log(det(omega.hat)) + sum(diag(sigma%*%omega.hat)) - dim(omega.hat)[1] )
}

glasso.pseudo <- function(omega,omega_path,num.of.solutions){
  sigma <- solve(omega)
  loss <- rep(0,num.of.solutions)
  for (i in 1:num.of.solutions){
    loss[i] <- loss.like(omega_path[,((i-1)*p+1):(i*p)],sigma)
  }
  i.best <- which.min(loss)
  list(omega.best=omega_path[,((i.best-1)*p+1):(i.best*p)],loss.best=loss[i.best],idx.best=i.best)
}

glasso.cv.tmp <- function(sigma.test,omega_path,num.of.grid){
  p <- dim(sigma.test)[1]
  loss <- rep(0,num.of.grid)
  for (i in 1:num.of.grid){
    loss[i] <- loss.like(omega_path[,((i-1)*p+1):(i*p)],sigma.test)
  }
  loss
}

glasso.cv <- function(sim, lam, tau, num.fold=5,method=1){
  quantile_alpha <- function(x){
    quantile(x,quan)
  }
  d = dim(sim$data)[2]
  n = dim(sim$data)[1]
  chunk.size = ceiling(n/num.fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  num.of.grid <- length(lam)*length(tau)
  omega.med <- matrix(0,d*d,length(partition))
  idx.best <- rep(0,num.fold)
  best.tau = rep(0,num.fold)
  for (i in 1:length(partition)){
    # cat("CV fold: ",i,"\n")
    val.idx <- partition[[i]]
    S.test <- cor(sim$data[val.idx,])
    S.train <- cor(sim$data[-val.idx,])
    omega.path <- glasso_path_new_Rwrapper(S.train,lam,tau,method=method)
    idx.best.i <- which.min(glasso.cv.tmp(S.test,omega.path$Omega_dc,num.of.grid))
    idx.best[i] <- idx.best.i
    tau.idx = idx.best.i %% length(tau)
    if (tau.idx == 0) { tau.idx = length(tau) }
    best.tau[i] = tau[tau.idx]
    omega.med[,i] <- as.vector(omega.path$Omega_dc[,((idx.best.i-1)*d+1):(idx.best.i*d)])
  }
  omega.median <- matrix(apply(omega.med,1,median),d,d)
  omega.mean <- matrix(apply(omega.med,1,mean),d,d)
  list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,tau.best=median(best.tau))
}


#' Constrained graphical lasso estimation using the truncated \eqn{\ell_1} penalty.
#'
#' @description
#' This function computes the constrained graphical lasso estimator using the truncated \eqn{\ell_1} penalty.
#'
#' @param sim A list containing a \eqn{n \times d} data matrix
#' @param bound a vector of upper bounds for the constraint
#' @param tau tuning parameter for the truncated \eqn{\ell_1} penalty
#' @param num.fold fold number for cross validation
#' @return cross validation score.
#' @export
glasso_nonconvex_constrained_cv = function(sim, bound=NULL, tau=NULL, num.fold=5){
  # rho=1.0,alpha=1.5, eps_abs=1e-5, eps_rel=1e-6, N.iter=1e3,dc.iter=1e1){

  quantile_alpha <- function(x)
  {
    quantile(x,quan)
  }

  if (is.null(bound))
  {
    bound = p*c(6,7,8,9,10,11,12)
  } else
  {
    bound = sort(bound)
  }
  if (is.null(tau))
  {
    tau = 10^(seq(0,-10,length.out=40)) ##### default warmstart strategy
  } else
  {
    if (length(tau) == 1) {
      tau <- 10^(seq(0,min(log10(tau),0),length.out=40)) ########## user-specified upper bound for tau
    } else
    {
      tau <- sort(tau,decreasing=TRUE) ########### user-specified sequence of tau for warmstart
    }
  }
  d = dim(sim$data)[2]
  n = dim(sim$data)[1]
  chunk.size = ceiling(n/num.fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  num.of.grid <- length(bound)
  idx.best <- rep(0,num.fold)
  cv.score.mat <- matrix(0,num.fold,length(bound))

  for (i in 1:length(partition))
  {
    cat("CV fold: ",i,"\n")
    val.idx <- partition[[i]]
    S.test <- cor(sim$data[val.idx,])
    S.train <- cor(sim$data[-val.idx,])
    omega.path <- glasso_nonconvex_constrained_path(S.train,bound,tau)
    small.tau.index <- rep(0,p*length(bound))
    for (j in seq_along(bound)){
      small.tau.index[((j-1)*p+1):(j*p)] = (j-1)*length(tau)*p + (p*(length(tau)-1)+1):(length(tau)*p)
    }
    omega.small.tau.path <- omega.path$Omega_dc[,small.tau.index]
    cv.score = glasso.cv.tmp(S.test,omega.small.tau.path,num.of.grid)
    idx.best.i <- which.min(cv.score)
    idx.best[i] <- idx.best.i
    cv.score.mat[i, ] = cv.score
  }
  return (cv.score.mat)
}


inference_LR <- function(sim, lam, tau, para_index, method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  sample.size = dim(sim$data)[1]
  out <- glasso.cv(sim, lam, tau, method=method)
  omega.best <- out$omega.mean
  best.index <- out$idx.best

  ########## get estimate of tau ############
  best.tau = out$tau.best
  ############ threshold at tau ##############
  omega.best.thred = omega.best
  omega.best.thred[abs(omega.best) < best.tau] = .0
  omega.null = omega.alt = NULL

  if (abs(omega.best.thred[para_index[1],para_index[2]]) > 1e-10){
    free_para_index <- which(abs(omega.best.thred)>1e-10) - 1
    out <- .C("glasso_cMLE",as.double(sim$sigmahat),
              as.integer(free_para_index),
              as.integer(length(free_para_index)),
              as.integer(p),
              as.integer(N.iter), as.double(eps_abs),
              as.double(eps_rel),
              as.double(rho), as.double(alpha), Omega=as.double(omega.best))
    omega.alt = matrix(out$Omega,p,p)
    omega.best.thred[para_index[1],para_index[2]] = omega.best.thred[para_index[2],para_index[1]] = 0
    free_para_index = which(abs(omega.best.thred) > 1e-10) - 1
    out <- .C("glasso_cMLE",as.double(sim$sigmahat),
              as.integer(free_para_index),
              as.integer(length(free_para_index)),
              as.integer(p),
              as.integer(N.iter), as.double(eps_abs),
              as.double(eps_rel),
              as.double(rho), as.double(alpha), Omega=as.double(omega.best))
    omega.null = matrix(out$Omega,p,p)
  } else {
    free_para_index = which(abs(omega.best.thred) > 1e-10) - 1
    out <- .C("glasso_cMLE",as.double(sim$sigmahat),
              as.integer(free_para_index),
              as.integer(length(free_para_index)),
              as.integer(p),
              as.integer(N.iter), as.double(eps_abs),
              as.double(eps_rel),
              as.double(rho), as.double(alpha), Omega=as.double(omega.best))
    omega.null = matrix(out$Omega,p,p)
    omega.best.thred[para_index[1],para_index[2]] = omega.best.thred[para_index[2],para_index[1]] = .1
    free_para_index = which(abs(omega.best.thred) > 1e-10) - 1
    out <- .C("glasso_cMLE",as.double(sim$sigmahat),
              as.integer(free_para_index),
              as.integer(length(free_para_index)),
              as.integer(p),
              as.integer(N.iter), as.double(eps_abs),
              as.double(eps_rel),
              as.double(rho), as.double(alpha), Omega=as.double(omega.best))
    omega.alt = matrix(out$Omega,p,p)
  }
  #cat(omega.alt,"\n \n")
  LR <- sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
  list(omega.null = omega.null, omega.alt = omega.alt, LR = LR,best.index=best.index,best.tau=best.tau)
}


glasso.cv.LR.inference <- function(sim, lam, tau_div_lam, para_index, method=1, num.fold=5){
  quantile_alpha <- function(x){
    quantile(x,quan)
  }
  d = dim(sim$data)[2]
  n = dim(sim$data)[1]
  chunk.size = ceiling(n/num.fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  num.of.grid <- length(lam)*length(tau_div_lam)
  omega.med <- matrix(0,d*d,length(partition))
  idx.best <- rep(0,num.fold)
  for (i in 1:length(partition)){
    # cat("CV fold: ",i,"\n")
    val.idx <- partition[[i]]
    S.test <- cor(sim$data[val.idx,])
    S.train <- cor(sim$data[-val.idx,])
    omega.path <- glasso_path_inference_LR_Rwrapper(S.train,lam,para_index,tau_div_lam,method=method)
    idx.best.i <- which.min(glasso.cv.tmp(S.test,omega.path$Omega_dc,num.of.grid))
    idx.best[i] <- idx.best.i
    omega.med[,i] <- as.vector(omega.path$Omega_dc[,((idx.best.i-1)*d+1):(idx.best.i*d)])
  }
  omega.median <- matrix(apply(omega.med,1,median),d,d)
  omega.mean <- matrix(apply(omega.med,1,mean),d,d)
  list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best)
}

glasso.cv.LR.inference_v2 <- function(sim, lam, tau, para_index, method=1, num.fold=5){
  quantile_alpha <- function(x){
    quantile(x,quan)
  }
  d = dim(sim$data)[2]
  n = dim(sim$data)[1]
  chunk.size = ceiling(n/num.fold)
  partition <- split(1:n, ceiling(1:n/chunk.size))
  num.of.grid <- length(lam)*length(tau)
  omega.med <- matrix(0,d*d,length(partition))
  idx.best <- rep(0,num.fold)
  best.tau = rep(0,num.fold)
  for (i in 1:length(partition)){
    # cat("CV fold: ",i,"\n")
    val.idx <- partition[[i]]
    S.test <- cor(sim$data[val.idx,])
    S.train <- cor(sim$data[-val.idx,])
    omega.path <- glasso_path_inference_LR_v2(S.train,lam,para_index,tau,method=method)
    idx.best.i <- which.min(glasso.cv.tmp(S.test,omega.path$Omega_dc,num.of.grid))
    idx.best[i] <- idx.best.i
    tau.idx = idx.best.i %% length(tau)
    if (tau.idx == 0) { tau.idx = length(tau) }
    best.tau[i] = tau[tau.idx]
    omega.med[,i] <- as.vector(omega.path$Omega_dc[,((idx.best.i-1)*d+1):(idx.best.i*d)])
  }
  omega.median <- matrix(apply(omega.med,1,median),d,d)
  omega.mean <- matrix(apply(omega.med,1,mean),d,d)
  list(omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,tau.best=median(best.tau))
}

glasso.subsample <- function(sim,lam,tau,num.subsample=20,val.ratio=.4,method=2,quan=.95,verbose = TRUE){
  quantile_alpha <- function(x){
    quantile(x,quan)
  }
  d = dim(sim$data)[2]
  n = dim(sim$data)[1]
  num.of.grid <- length(lam)*length(tau)
  omega.med <- matrix(0,d*d,num.subsample)
  idx.best <- rep(0,num.subsample)
  for (i in 1:num.subsample){
    if (verbose){
      cat("conducting subsampling ", i,"\n")
      flush.console()
    }
    val.idx <- sample(1:n,ceiling(val.ratio*n))
    S.test <- cor(sim$data[val.idx,])
    S.train <- cor(sim$data[-val.idx,])
    omega.path <- glasso_path_new_Rwrapper(S.train,lam,tau,method=method)
    idx.best.i <- which.min(glasso.cv.tmp(S.test,omega.path$Omega_dc,num.of.grid))
    idx.best[i] <- idx.best.i
    omega.med[,i] <- as.vector(omega.path$Omega_dc[,((idx.best.i-1)*d+1):(idx.best.i*d)])
  }
  omega.quan <- matrix(apply(omega.med,1,quantile_alpha),d,d)
  omega.median <- matrix(apply(omega.med,1,median),d,d)
  omega.mean <- matrix(apply(omega.med,1,mean),d,d)
  list(omega.quan=omega.quan,omega.median=omega.median,omega.mean=omega.mean,idx.best=idx.best,omega.path=omega.med)
}

# inference based on likelihood ratio #
glasso_path_test_LR <- function(sim, lam, tau_div_lam, sample.size, para_index, method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  out <- glasso.cv.LR.inference(sim, lam, tau_div_lam, para_index,method=method)
  omega.best <- out$omega.mean
  best.index <- out$idx.best

  ########## get estimate of tau ############
  # tau_div_lam_tmp = tau_div_lam[best.index %% length(tau_div_lam)]
  # best.lam = lam[(best.index - best.index %% length(tau_div_lam))/length(tau_div_lam) + 1]
  # best.tau = mean(tau_div_lam_tmp * best.lam)
  best.tau = sqrt(log(p) / sample.size)
  ############ threshold at tau ##############
  free_para_index <- which(abs(omega.best)>best.tau) - 1
  out <- .C("glasso_cMLE",as.double(sim$sigmahat),
            as.integer(free_para_index),
            as.integer(length(free_para_index)),
            as.integer(p),
            as.integer(N.iter), as.double(eps_abs),
            as.double(eps_rel),
            as.double(rho), as.double(alpha), Omega=as.double(omega.best))
  omega.null <- matrix(out$Omega,p,p)
  #cat(omega.null,"\n \n")
  #cat(omega.best,"\n \n")
  #omega.null <- omega.best
  para_index <- para_index - 1
  free_para_index <- c(free_para_index,para_index[1]+para_index[2]*p,para_index[2]+para_index[1]*p)
  out <- .C("glasso_cMLE",as.double(sim$sigmahat),
            as.integer(free_para_index),
            as.integer(length(free_para_index)),
            as.integer(p),
            as.integer(N.iter), as.double(eps_abs),
            as.double(eps_rel),
            as.double(rho), as.double(alpha), Omega=as.double(omega.best))
  omega.alt <- matrix(out$Omega,p,p)
  #cat(omega.alt,"\n \n")
  LR <- sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
  list(omega.null = omega.null, omega.alt = omega.alt, LR = LR,best.index=best.index)
}

# inference based on likelihood ratio (tau as input istead of tau_div_lam) #
glasso_path_test_LR_v2 <- function(sim, lam, tau, sample.size, para_index, method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  out <- glasso.cv.LR.inference_v2(sim, lam, tau, para_index,method=method)
  omega.best <- out$omega.mean
  best.index <- out$idx.best

  ########## get estimate of tau ############
  # tau_div_lam_tmp = tau_div_lam[best.index %% length(tau_div_lam)]
  # best.lam = lam[(best.index - best.index %% length(tau_div_lam))/length(tau_div_lam) + 1]
  # best.tau = mean(tau_div_lam_tmp * best.lam)
  best.tau = out$tau.best
  ############ threshold at tau ##############
  free_para_index <- which(abs(omega.best)>best.tau) - 1
  out <- .C("glasso_cMLE",as.double(sim$sigmahat),
            as.integer(free_para_index),
            as.integer(length(free_para_index)),
            as.integer(p),
            as.integer(N.iter), as.double(eps_abs),
            as.double(eps_rel),
            as.double(rho), as.double(alpha), Omega=as.double(omega.best))
  omega.null <- matrix(out$Omega,p,p)
  #cat(omega.null,"\n \n")
  #cat(omega.best,"\n \n")
  #omega.null <- omega.best
  para_index <- para_index - 1
  free_para_index <- c(free_para_index,para_index[1]+para_index[2]*p,para_index[2]+para_index[1]*p)
  out <- .C("glasso_cMLE",as.double(sim$sigmahat),
            as.integer(free_para_index),
            as.integer(length(free_para_index)),
            as.integer(p),
            as.integer(N.iter), as.double(eps_abs),
            as.double(eps_rel),
            as.double(rho), as.double(alpha), Omega=as.double(omega.best))
  omega.alt <- matrix(out$Omega,p,p)
  #cat(omega.alt,"\n \n")
  LR <- sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
  list(omega.null = omega.null, omega.alt = omega.alt, LR = LR,best.index=best.index,best.tau=best.tau)
}


