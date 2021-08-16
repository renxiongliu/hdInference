# source("util_funcs.R")
# dyn.load("./src/glasso.so")

glasso <- function(Sigma.hat, lam, rho=1.0,alpha=1.5,eps_abs = 1e-4,eps_rel=1e-5,N.iter=1e3)
{
  p <- dim(Sigma.hat)[1]
  Omega <- diag(1/diag(Sigma.hat))
  Gamma <- matrix(0,p,p)
  Delta_old <- Omega
  Delta <- soft.threshold.off.diag(Omega+Gamma,lam/rho)
  Lambda <- rep(0,p)
  for (j in 1:N.iter){
    # Omega update #
    tmp <- eigen(rho*(Delta - Gamma) - Sigma.hat,symmetric=TRUE)
    for (i in 1:p){
      Lambda[i] <- (tmp$values[i] + sqrt(tmp$values[i]^2+4*rho))/(2*rho)
    }
    Omega <- tmp$vectors %*% diag(Lambda) %*% t(tmp$vectors)

    # Delta update #
    Delta_old <- Delta
    Delta <- soft.threshold.off.diag(alpha*Omega+(1-alpha)*Delta+Gamma,lam/rho)
    # Gamma update #
    Gamma <- Gamma + alpha*Omega+(1-alpha)*Delta_old - Delta
    # check stopping rule #
    dual_res <- rho*norm(Delta-Delta_old,"F")
    primal_res <- norm(Omega - Delta,"F")
    eps_primal <- eps_abs*p + eps_rel*max(norm(Omega,"F"),norm(Delta,"F"))
    eps_dual <- eps_abs*p + eps_rel*norm(Gamma,"F")*rho
    # cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
    if ((primal_res < eps_primal) && (dual_res <- eps_dual)){
      cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
      break;
    }
  }
  Delta
}

glasso_dc <- function(Sigma.hat, lam, rho=1.0,alpha=1.5,tau=1e-2,eps_abs=1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=5)
{
  p <- dim(Sigma.hat)[1]
  Omega <- diag(1/diag(Sigma.hat))
  Gamma <- matrix(0,p,p)
  Lambda <- matrix(0,p,p)
  Delta_old <- Omega
  Delta <- soft.threshold.off.diag(Omega+Gamma,lam/rho)
  D <- rep(0,p)
  for (j in 1:N.iter){
    # Omega update #
    tmp <- eigen(rho*(Delta - Gamma) - Sigma.hat,symmetric=TRUE)
    for (i in 1:p){
      D[i] <- (tmp$values[i] + sqrt(tmp$values[i]^2+4*rho))/(2*rho)
    }
    Omega <- tmp$vectors %*% diag(D) %*% t(tmp$vectors)

    # Delta update #
    Delta_old <- Delta
    Delta <- soft.threshold.off.diag(alpha*Omega+(1-alpha)*Delta_old+Gamma,lam/rho)

    # Gamma update #
    Gamma <- Gamma + alpha*Omega+(1-alpha)*Delta - Delta
    # check stopping rule #
    dual_res <- rho*norm(Delta-Delta_old,"F")
    primal_res <- norm(Omega - Delta,"F")
    eps_primal <- eps_abs*p + eps_rel*max(norm(Omega,"F"),norm(Delta,"F"))
    eps_dual <- eps_abs*p + eps_rel*norm(Gamma,"F")*rho
    # cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
    if ((primal_res < eps_primal) && (dual_res <- eps_dual)){
      cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
      break;
    }
  }
  for (k in 1:dc.iter){
    Lambda[] <- 0
    Lambda[abs(Delta)<tau] <- 1
    for (j in 1:N.iter){
      # Omega update #
      tmp <- eigen(rho*(Delta - Gamma) - Sigma.hat,symmetric=TRUE)
      for (i in 1:p){
        D[i] <- (tmp$values[i] + sqrt(tmp$values[i]^2+4*rho))/(2*rho)
      }
      Omega <- tmp$vectors %*% diag(D) %*% t(tmp$vectors)

      # Delta update #
      Delta_old <- Delta
      Delta <- soft.threshold(Omega+Gamma,Lambda/rho)

      # Gamma update #
      Gamma <- Gamma + Omega - Delta
      # check stopping rule #
      dual_res <- rho*norm(Delta-Delta_old,"F")
      primal_res <- norm(Omega - Delta,"F")
      eps_primal <- eps_abs*p + eps_rel*max(norm(Omega,"F"),norm(Delta,"F"))
      eps_dual <- eps_abs*p + eps_rel*norm(Gamma,"F")*rho
      # cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
      if ((primal_res < eps_primal) && (dual_res <- eps_dual)){
        cat("#iter",j,"primal residual:", primal_res," dual residual:",dual_res,"\n")
        break;
      }
    }
  }
  Delta
}


glasso_path_Rwrapper <- function(Sigma.hat, lam, method=1,rho=1.0,alpha=1.5,tau=1e-2,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
  Omega_dc <- Omega
  out_dc <- .C("glasso_nonconvex_path",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(tau),as.integer(p),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega=as.double(Omega),Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len))
}

glasso_path_new_Rwrapper <- function(Sigma.hat, lam, tau,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  tau_len <- length(tau)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_len),p,p*lam_len*tau_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_len),p,p*lam_len*tau_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_new",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(tau),as.integer(tau_len),as.integer(p),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega=as.double(Omega),Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_len))
}

glasso_nonconvex_constrained_path_R = function(Sigma.hat,bound,tau,
  free_para_index=NULL,zero_para_index=NULL,rho=1.0,alpha=1.5,
  eps_abs=1e-5,eps_rel=1e-6, N.iter=2e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  bound_len <- length(bound)
  tau_len <- length(tau)
  if (is.null(free_para_index))
  {
    num_of_free_para = 0
  } else
  {
    num_of_free_para <- length(free_para_index) / 2
    free_para_index <- free_para_index - 1 ## convert to C indexing ##
  }

  if (is.null(zero_para_index))
  {
    num_of_zero_para = 0
  } else
  {
    num_of_zero_para = length(zero_para_index) / 2
    zero_para_index = zero_para_index - 1
  }
  Omega <- matrix(rep(diag(1/diag(Sigma.hat)),bound_len),p,p*bound_len)
  Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),bound_len*tau_len),p,p*bound_len*tau_len)
  out_dc <- .C("glasso_nonconvex_constrained_path",as.double(Sigma.hat),
        as.integer(free_para_index), as.integer(num_of_free_para),
        as.integer(zero_para_index), as.integer(num_of_zero_para),
        as.integer(bound),as.integer(bound_len),as.double(tau),
        as.integer(tau_len),as.integer(p),as.integer(N.iter),
        as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),
        as.double(rho),as.double(alpha),Omega=as.double(Omega),Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*bound_len),Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*bound_len*tau_len))
}



########### seems not work ###########
glasso_nonconvex_constrained_single_warmstart = function(Sigma.hat,bound,tau=1e-8,free_para_index=NULL,
                                                         zero_para_index=NULL,rho=1.0,alpha=1.5,
                                                         eps_abs=1e-5,eps_rel=1e-6, N.iter=1e3,
                                                         dc.iter=1e1, num.of.warmstart=20)
{
  if (is.null(free_para_index))
  {
    num_of_free_para = 0
  } else
  {
    num_of_free_para <- length(free_para_index) / 2
    free_para_index <- free_para_index - 1 ## convert to C indexing ##
  }

  if (is.null(zero_para_index))
  {
    num_of_zero_para = 0
  } else
  {
    num_of_zero_para = length(zero_para_index) / 2
    zero_para_index = zero_para_index - 1
  }
  tau_max = max(tau, 1)
  omega = glasso_nonconvex_constrained_path_R(Sigma.hat,bound,tau_max,
                                               free_para_index=free_para_index,
                                               zero_para_index=zero_para_index,
                                               rho=rho,alpha=alpha, eps_abs=eps_abs,
                                               eps_rel=eps_rel, N.iter=N.iter,dc.iter=dc.iter)$Omega_dc
  if (tau_max != tau)
  {
      tau_seq = seq(tau_max,tau,length.out=num.of.warmstart)
      omega = diag(1 / diag(Sigma.hat))
      for (tau.i in tau_seq)
      {
        out = .C("glasso_nonconvex_constrained_single",as.double(Sigma.hat),
               as.integer(free_para_index), as.integer(num_of_free_para),
               as.integer(zero_para_index), as.integer(num_of_zero_para),
               as.double(bound),as.double(tau.i),
               as.integer(p),as.integer(N.iter),
               as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),
               as.double(rho),as.double(alpha),Omega_dc=as.double(omega))
        omega = matrix(out$Omega_dc,p,p)
      }
  }
  omega
}

glasso_path_v2_Rwrapper <- function(Sigma.hat, lam, tau_div_lam,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  tau_div_lam_len <- length(tau_div_lam)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_v2",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),as.double(tau_div_lam),as.integer(tau_div_lam_len),as.integer(p),as.integer(N.iter),as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),as.double(rho),as.double(alpha),as.integer(method),Omega=as.double(Omega),Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_div_lam_len))
}


inference_general_LR <- function(sim, para_index, lam=NULL, tau=NULL,method=2,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  sample.size = dim(sim$data)[1]
  if (is.null(lam)) {
    lam = sqrt(2*log(p)/sample.size)
    # s.hat = sqrt(sample.size)/log(p)
    # B = qt(1-s.hat/(2*p),df=sample.size-1)
    # lam = B / sqrt(sample.size-1+B^2)
  }
  if (is.null(tau)) { tau = .5*sqrt(log(p)/sample.size) }
  out = glasso_path_new_Rwrapper(sim$sigmahat, lam, tau, method=method,rho=rho,alpha=alpha, eps_abs=eps_abs,eps_rel=eps_rel,N.iter=N.iter,dc.iter=dc.iter)
  omega.est = out$Omega_dc
  thred_level = .5*sqrt(omega.est^2+diag(omega.est)%*% t(diag(omega.est))) * sqrt(log(p)/sample.size)
  #thred_level = 1e-10
  omega.est.thred = soft.threshold(omega.est, thred_level)
  out.tmp = report.roc(omega.est.thred,sim$omega)
  cat("T1 error: ", out.tmp$T1, out.tmp$T1_rate,"\t","T2 error: ", out.tmp$T2, out.tmp$T2_rate,"\n")
  diag(omega.est.thred) = 1.0 # any nonzero number will do #
  S_union_B = omega.est.thred
  S_minus_B = omega.est.thred
  for (j in 1:(length(para_index)/2))
  {
    if (abs(omega.est.thred[para_index[2*j-1],para_index[2*j]]) > 1e-10) {
      S_minus_B[para_index[2*j-1],para_index[2*j]] = 0
      S_minus_B[para_index[2*j],para_index[2*j-1]] = 0
    } else {
      S_union_B[para_index[2*j-1],para_index[2*j]] = 1
      S_union_B[para_index[2*j],para_index[2*j-1]] = 1
    }
  }
  free_para_index = which(abs(S_union_B)>1e-10) - 1
  omega.alt = matrix(.C("glasso_cMLE",as.double(sim$sigmahat),
                 as.integer(free_para_index),
                 as.integer(length(free_para_index)),
                 as.integer(p),
                 as.integer(N.iter), as.double(eps_abs),
                 as.double(eps_rel),
                 as.double(rho), as.double(alpha), Omega=as.double(omega.est))$Omega,p,p)
  free_para_index = which(abs(S_minus_B)>1e-10) - 1
  omega.null = matrix(.C("glasso_cMLE",as.double(sim$sigmahat),
                        as.integer(free_para_index),
                        as.integer(length(free_para_index)),
                        as.integer(p),
                        as.integer(N.iter), as.double(eps_abs),
                        as.double(eps_rel),
                        as.double(rho), as.double(alpha), Omega=as.double(omega.est))$Omega,p,p)
  LR <- sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
  LR = (LR - length(para_index)/2) / sqrt(length(para_index))
  return (list(omega.null = omega.null, omega.alt = omega.alt, LR = LR))
}

#' Likelihood ratio test for a subset of parameters in Gaussian graphical model
#'
#' @description
#' This function performs likelihood ratio test for a constrained graphical lasso estimator using the truncated \eqn{\ell_1} penalty.
#'
#' @param sim A list containing a \eqn{n \times d} data matrix and a covariance matrix estimator
#' @param para_index indices of the parameters of interest
#' @param bound a vector of upper bounds for the constraint
#' @param tau tuning parameter for the truncated \eqn{\ell_1} penalty
#' @param inference_type either 'LR' or 'LR_gen' with the former based on chi-square test and the later based on normal approximation
#' @return omega.null:
#' @return omega.alt:
#' @return LR:
#' @export
#'
inference_constrained <- function(sim, para_index, bound=NULL, tau=NULL, inference_type="LR",
 rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=2e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  sample.size = dim(sim$data)[1]
  if (is.null(bound)) { bound = 4*p } else { if (length(bound) > 1) stop("only support single bound! \n") }
  if (is.null(tau))
  {
    tau = 10^(seq(0,-10,length.out=40)) ##### default warmstart strategy
    # tau = 10^(seq(1,-10,length.out=40)) ##### added later on Feb 27, 2017 because solutions become infeasible when tau = 1.0
  }
  else
  {
    if (length(tau) == 1) {
      tau <- 10^(seq(0,min(log10(tau),0),length.out=40)) ########## user-specified upper bound for tau
    } else
    {
      tau <- sort(tau,decreasing=TRUE) ########### user-specified sequence of tau for warmstart
    }
  }
  omega_warmstart_path = glasso_nonconvex_constrained_path_R(sim$sigmahat,bound,tau,free_para_index=para_index,zero_para_index=para_index)$Omega_dc
  omega.null = omega_warmstart_path[,(p*(length(tau)-1)+1):(p*length(tau))]

  ## compute alternative
  out = .C("glasso_nonconvex_constrained_single",as.double(sim$sigmahat),
           as.integer(para_index-1), as.integer(length(para_index)/2),
           as.integer(NULL), as.integer(0),
           as.integer(bound),as.double(min(tau)),
           as.integer(p),as.integer(N.iter),
           as.integer(dc.iter),as.double(eps_abs),as.double(eps_rel),
           as.double(rho),as.double(alpha),Omega_dc=as.double(omega.null))

  omega.alt = matrix(out$Omega_dc,p,p)
  LR = sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
  if (inference_type == "LR_gen")
  {
    LR = (LR - length(para_index)/2) / sqrt(length(para_index))
  }
  return (list(omega.null = omega.null, omega.alt = omega.alt, LR = LR))
}

inference <- function(sim, para_index, lam=NULL, tau=NULL, inference_type="Normal", method=2,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(sim$sigmahat)[1]
  sample.size = dim(sim$data)[1]
  if (is.null(lam)) {
    lam = sqrt(2*log(p)/sample.size)
    # s.hat = sqrt(sample.size)/log(p)
    # B = qt(1-s.hat/(2*p),df=sample.size-1)
    # lam = B / sqrt(sample.size-1+B^2)
  }
  if (is.null(tau)) { tau = .5*sqrt(log(p)/sample.size) }
  out = glasso_path_new_Rwrapper(sim$sigmahat, lam, tau,
        method=method,rho=rho,alpha=alpha,
        eps_abs=eps_abs,eps_rel=eps_rel,
        N.iter=N.iter,dc.iter=dc.iter)
  omega.est = out$Omega_dc
  thred_level = .5*sqrt(omega.est^2+diag(omega.est)%*% t(diag(omega.est))) * sqrt(log(p)/sample.size)
  #thred_level = 1e-10
  omega.est.thred = soft.threshold(omega.est, thred_level)
  out.tmp = report.roc(omega.est.thred,sim$omega)
  cat("T1 error: ", out.tmp$T1, out.tmp$T1_rate,"\t","T2 error: ", out.tmp$T2, out.tmp$T2_rate,"\n")
  diag(omega.est.thred) = 1.0 # any nonzero number will do #
  if (inference_type=="LR") {
    omega.null = omega.alt = NULL
    if (abs(omega.est.thred[para_index[1],para_index[2]]) > 1e-10){
      free_para_index <- which(abs(omega.est.thred)>1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.alt = matrix(out$Omega,p,p)
      omega.est.thred[para_index[1],para_index[2]] = omega.est.thred[para_index[2],para_index[1]] = 0
      free_para_index = which(abs(omega.est.thred) > 1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.null = matrix(out$Omega,p,p)
    } else {
      free_para_index = which(abs(omega.est.thred) > 1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.null = matrix(out$Omega,p,p)
      omega.est.thred[para_index[1],para_index[2]] = omega.est.thred[para_index[2],para_index[1]] = .1
      free_para_index = which(abs(omega.est.thred) > 1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.alt = matrix(out$Omega,p,p)
    }
    LR <- sample.size*(log(det(omega.alt)) - sum(diag(sim$sigmahat%*%omega.alt)) - log(det(omega.null)) + sum(diag(sim$sigmahat%*%omega.null)))
    return (list(omega.null = omega.null, omega.alt = omega.alt, LR = LR))
  } else if ( inference_type=="Normal")
  {
    omega.MLE = NULL
    if (abs(omega.est.thred[para_index[1],para_index[2]]) > 1e-10){
      free_para_index <- which(abs(omega.est.thred)>1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.MLE = matrix(out$Omega,p,p)
    } else {
      omega.est.thred[para_index[1],para_index[2]] = omega.est.thred[para_index[2],para_index[1]] = .1
      free_para_index = which(abs(omega.est.thred) > 1e-10) - 1
      out <- .C("glasso_cMLE",as.double(sim$sigmahat),
                as.integer(free_para_index),
                as.integer(length(free_para_index)),
                as.integer(p),
                as.integer(N.iter), as.double(eps_abs),
                as.double(eps_rel),
                as.double(rho), as.double(alpha), Omega=as.double(omega.est))
      omega.MLE = matrix(out$Omega,p,p)
    }
    idx = which(abs(omega.MLE) > 1e-10, arr.ind=TRUE)
    fisher.dim = dim(idx)[1]
    I_oracle_inv = matrix(0,fisher.dim,fisher.dim)
    for (k in 1:fisher.dim){
      for (kk in 1:fisher.dim){
        I_oracle_inv[k,kk] = sim$sigma[idx[k,2],idx[kk,1]] * sim$sigma[idx[kk,2],idx[k,1]]
      }
    }
    I_oracle = solve(I_oracle_inv)
    idx.half = which(((abs(omega.MLE) > 1e-10) & lower.tri(sim$omega,diag=TRUE)), arr.ind=TRUE)
    fisher.dim.half = dim(idx.half)[1]
    fisher.variances = rep(0,fisher.dim.half)
    for (k in 1:fisher.dim.half){
      i = idx.half[k,1]
      j = idx.half[k,2]
      index = NULL
      for (kk in 1:fisher.dim){
        if ((idx[kk,1] == i) && ((idx[kk,2] == j))) index = c(index,kk)
        if ((idx[kk,2] == i) && ((idx[kk,1] == j))) index = c(index,kk)
      }
      if (i == j) {
        fisher.variances[k] = 2 * I_oracle[index[1],index[1]]
      } else {
        fisher.variances[k] = I_oracle[index[1],index[2]]+.5*(I_oracle[index[1],index[1]] + I_oracle[index[2],index[2]])
      }
    }

    ########## get center and variance of CI ############
    omega.std = NULL
    omega.center = omega.MLE[para_index[1],para_index[2]]
    for (k in 1:fisher.dim.half){
      i = idx.half[k,1]
      j = idx.half[k,2]
      if ((para_index[1] == i) && ((para_index[2] == j))) omega.std = sqrt(fisher.variances[k]/sample.size)
      if ((para_index[1] == j) && ((para_index[2] == i))) omega.std = sqrt(fisher.variances[k]/sample.size)
    }
    return (list(omega.MLE = omega.MLE,omega.center=omega.center,omega.std=omega.std))
  } else {
    stop("inference type not supported! please specify inference type to be either LR or Normal. \n")
  }
}


# inference based on asymptotic normality #
glasso_path_inference_Rwrapper <- function(Sigma.hat, lam, free_para_index, tau_div_lam,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  num_of_free_para <- length(free_para_index) / 2
  free_para_index <- free_para_index - 1 ## convert to C indexing ##
  tau_div_lam_len <- length(tau_div_lam)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_inference",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),
               as.integer(free_para_index), as.integer(num_of_free_para), as.double(tau_div_lam),as.integer(tau_div_lam_len),
               as.integer(p),as.integer(N.iter),as.integer(dc.iter),
               as.double(eps_abs),as.double(eps_rel),as.double(rho),
               as.double(alpha),as.integer(method),Omega=as.double(Omega),
               Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),
       Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_div_lam_len))
}

# inference based on likelihood ratio #
glasso_path_inference_LR_Rwrapper <- function(Sigma.hat, lam, zero_para_index, tau_div_lam,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  zero_para_index <- zero_para_index - 1 ## convert to C indexing ##
  tau_div_lam_len <- length(tau_div_lam)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_div_lam_len),p,p*lam_len*tau_div_lam_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_inference_LR",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),
               as.integer(zero_para_index), as.integer(length(zero_para_index)/2), as.double(tau_div_lam),as.integer(tau_div_lam_len),
               as.integer(p),as.integer(N.iter),as.integer(dc.iter),
               as.double(eps_abs),as.double(eps_rel),as.double(rho),
               as.double(alpha),as.integer(method),Omega=as.double(Omega),
               Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),
       Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_div_lam_len))
}


# inference based on likelihood ratio #
glasso_path_inference_LR_v2 <- function(Sigma.hat, lam, zero_para_index, tau,method=1,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1)
{
  p <- dim(Sigma.hat)[1]
  lam_len <- length(lam)
  zero_para_index <- zero_para_index - 1 ## convert to C indexing ##
  tau_len <- length(tau)
  if (lam[1] >= lam[lam_len]) {
    Omega <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(1/diag(Sigma.hat)),lam_len*tau_len),p,p*lam_len*tau_len)
  } else {
    Omega <- matrix(rep(diag(p),lam_len),p,p*lam_len)
    Omega_dc <- matrix(rep(diag(p),lam_len*tau_len),p,p*lam_len*tau_len)
    Omega[,1:p] <- solve(Sigma.hat + 1e-8*diag(p))
    Omega_dc[,1:p] <- Omega[,1:p]
  }
  out_dc <- .C("glasso_nonconvex_path_inference_LR_v2",as.double(Sigma.hat),as.double(lam),as.integer(lam_len),
               as.integer(zero_para_index), as.integer(length(zero_para_index)/2), as.double(tau),as.integer(tau_len),
               as.integer(p),as.integer(N.iter),as.integer(dc.iter),
               as.double(eps_abs),as.double(eps_rel),as.double(rho),
               as.double(alpha),as.integer(method),Omega=as.double(Omega),
               Omega_dc=as.double(Omega_dc))
  list(Omega=matrix(out_dc$Omega,nrow=p,ncol=p*lam_len),
       Omega_dc=matrix(out_dc$Omega_dc,nrow=p,ncol=p*lam_len*tau_len))
}

## LR might be negative ... ##
inference_v2 <- function(sim, para_index, lam=NULL, tau=NULL, inference_type="Normal", method=2,rho=1.0,alpha=1.5,eps_abs = 1e-5,eps_rel=1e-6,N.iter=1e3,dc.iter=1e1){
  p <- dim(sim$sigmahat)[1]
  sample.size = dim(sim$data)[1]
  if (is.null(lam)) {
    lam = 2*sqrt(log(p)/sample.size)
    # s.hat = sqrt(sample.size)/log(p)
    # B = qt(1-s.hat/(2*p),df=sample.size-1)
    # lam = B / sqrt(sample.size-1+B^2)
  }
  if (is.null(tau)) { tau = .5*sqrt(log(p)/sample.size) }
  omega_L1 = omega_trL1 = diag(p)
  gamma_L1 = gamma_trL1 = matrix(0,p,p)
  free_index = para_index - 1
  out_free <- .C("glasso_nonconvex_free",
                 as.double(sim$sigmahat),
                 as.double(lam),
                 as.integer(free_index),
                 as.integer(length(free_index)/2),
                 as.double(tau),
                 as.integer(p),
                 as.integer(N.iter), as.integer(dc.iter),
                 as.double(eps_abs),
                 as.double(eps_rel),
                 as.double(rho), as.double(alpha),
                 as.integer(method),
                 gamma_L1 = as.double(gamma_L1),
                 gamma_trL1 = as.double(gamma_trL1),
                 omega_L1 = as.double(omega_L1),
                 omega_trL1=as.double(omega_trL1))
  omega_L1_free = matrix(out_free$omega_L1,p,p)
  omega_trL1_free = matrix(out_free$omega_trL1,p,p)
  #thred_level = .5*sqrt(omega.est^2+diag(omega.est)%*% t(diag(omega.est))) * sqrt(log(p)/sample.size)
  #thred_level = 1e-10
  #omega.est.thred = soft.threshold(omega.est, thred_level)
  out.tmp = report.roc(omega_trL1_free,sim$omega)
  cat("T1 error: ", out.tmp$T1, out.tmp$T1_rate,"\t","T2 error: ", out.tmp$T2, out.tmp$T2_rate,"\n")
  if (inference_type=="LR") {
    zero_index = para_index - 1
    out_zero <- .C("glasso_nonconvex_zero",
                   as.double(sim$sigmahat),
                   as.double(lam),
                   as.integer(zero_index),
                   as.integer(length(zero_index)/2),
                   as.double(tau),
                   as.integer(p),
                   as.integer(N.iter), as.integer(dc.iter),
                   as.double(eps_abs),
                   as.double(eps_rel),
                   as.double(rho), as.double(alpha),
                   as.integer(method),
                   gamma_L1 = as.double(out_free$gamma_L1),
                   gamma_trL1 = as.double(out_free$gamma_trL1),
                   omega_L1 = as.double(out_free$omega_L1),
                   omega_trL1=as.double(out_free$omega_trL1))
    omega_L1_zero = matrix(out_zero$omega_L1,p,p)
    omega_trL1_zero = matrix(out_zero$omega_trL1,p,p)
    LR <- sample.size*(log(det(omega_trL1_free)) - sum(diag(sim$sigmahat%*%omega_trL1_free)) - log(det(omega_trL1_zero)) + sum(diag(sim$sigmahat%*%omega_trL1_zero)))
    return (list(omega.null = omega_trL1_zero, omega.alt = omega_trL1_free, LR = LR))
  } else if ( inference_type=="Normal") {
    omega.MLE = omega_trL1_free
    idx = which(abs(omega.MLE) > 1e-10, arr.ind=TRUE)
    fisher.dim = dim(idx)[1]
    I_oracle_inv = matrix(0,fisher.dim,fisher.dim)
    for (k in 1:fisher.dim){
      for (kk in 1:fisher.dim){
        I_oracle_inv[k,kk] = sim$sigma[idx[k,2],idx[kk,1]] * sim$sigma[idx[kk,2],idx[k,1]]
      }
    }
    I_oracle = solve(I_oracle_inv)
    idx.half = which(((abs(omega.MLE) > 1e-10) & lower.tri(sim$omega,diag=TRUE)), arr.ind=TRUE)
    fisher.dim.half = dim(idx.half)[1]
    fisher.variances = rep(0,fisher.dim.half)
    for (k in 1:fisher.dim.half){
      i = idx.half[k,1]
      j = idx.half[k,2]
      index = NULL
      for (kk in 1:fisher.dim){
        if ((idx[kk,1] == i) && ((idx[kk,2] == j))) index = c(index,kk)
        if ((idx[kk,2] == i) && ((idx[kk,1] == j))) index = c(index,kk)
      }
      if (i == j) {
        fisher.variances[k] = 2 * I_oracle[index[1],index[1]]
      } else {
        fisher.variances[k] = I_oracle[index[1],index[2]]+.5*(I_oracle[index[1],index[1]] + I_oracle[index[2],index[2]])
      }
    }

    ########## get center and variance of CI ############
    omega.std = NULL
    omega.center = omega.MLE[para_index[1],para_index[2]]
    for (k in 1:fisher.dim.half){
      i = idx.half[k,1]
      j = idx.half[k,2]
      if ((para_index[1] == i) && ((para_index[2] == j))) omega.std = sqrt(fisher.variances[k]/sample.size)
      if ((para_index[1] == j) && ((para_index[2] == i))) omega.std = sqrt(fisher.variances[k]/sample.size)
    }
    return (list(omega.MLE = omega.MLE,omega.center=omega.center,omega.std=omega.std))
  } else {
    stop("inference type not supported! please specify inference type to be either LR or Normal. \n")
  }
}

plot.omega = function(omega,omega_true){
  gcinfo(FALSE)
  d <- dim(omega)[1]
  theta <- matrix(0,d,d)
  theta_true <- theta
  theta[abs(matrix(omega))>1e-10] = 1
  theta_true[abs(matrix(omega_true))>1e-10] = 1
  par = par(mfrow = c(1, 2), pty = "s", omi=c(.3,.3,.3,.3), mai = c(.3,.3,.3,.3))
  image(as.matrix(theta), col = gray.colors(256),  main = "Estimated")
  image(as.matrix(theta_true), col = gray.colors(256),  main = "Truth")
}

