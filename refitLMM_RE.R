##### refit and calculate the loglikelihood  
### refitmodel for RE ###

library(Matrix)

refitLMM_RE <- function(indexy, y, x, group, par0 = NULL)
{
  
  nt <- length(y)
  np <- ncol(x)
  
  uindex <- unique(indexy)
  n <- length(uindex)
  
  ng <- length(unique(group))
  
  nJ <- rep(0,n)
  for(i in 1:n)
  {
    nJ[i] <- sum(indexy==uindex[i])
  }
  ## REML 
  
  loglike_reml <- function(par)
  {
    sigmav2 = exp(par[1])
    sigmae2 = exp(par[2:(ng+1)])
    gammah <- sigmav2/sigmae2[group]
    detV <- sum(nJ*log(sigmae2[group])) + sum(log(1+gammah*nJ))
    xvx <- 0 
    Vinvlist <- list()
    
    for(i in 1:n)
    {
      indexi <- indexy == uindex[i]
      ni <- nJ[i]
      mat1 <- matrix(1,ni,ni)
      diag1 <- diag(1,ni)
      
      S1 <- mat1 *sigmav2 + diag1 * sigmae2[group[i]]
      S1inv <- (diag1 - mat1*sigmav2/(sigmae2[group[i]] + ni*sigmav2))/sigmae2[group[i]]
      Vinvlist[[i]] <- S1inv
      
    }
    
    Vinv <- bdiag(Vinvlist)
    
    detXV <- determinant(t(x) %*% Vinv %*% x, logarithm = TRUE)$modulus
    
    betahat <- solve(t(x)%*%Vinv %*% x) %*% t(x) %*% Vinv %*% y
    
    value <- as.numeric(detV + detXV + t(y - x %*% betahat) %*% Vinv %*%  (y - x %*% betahat))
    return(value)
    
  }
  
  if(is.null(par0))
  {
    res <- optim(loglike_reml,par = rep(0,1+ng), method = "L-BFGS-B", 
                 lower = rep(-10,1+ng), upper = rep(10,1+ng))
  }
  else{
    res <- optim(loglike_reml,par = par0, method = "L-BFGS-B", 
                 lower = rep(-10,1+ng), upper = rep(10,1+ng))
  }
  
  
  sig2est <- exp(res$par)
  
  sigmav2 = sig2est[1]
  sigmae2 = sig2est[2:(ng+1)]
  gammah <- sigmav2/sigmae2[group]
  detV <- sum(nJ*log(sigmae2[group])) + sum(log(1+gammah*nJ))
  xvx <- 0 
  Vinvlist <- list()
  
  for(i in 1:n)
  {
    indexi <- indexy == uindex[i]
    ni <- nJ[i]
    mat1 <- matrix(1,ni,ni)
    diag1 <- diag(1,ni)
    
    S1 <- mat1 *sigmav2 + diag1 * sigmae2[group[i]]
    S1inv <- (diag1 - mat1*sigmav2/(sigmae2[group[i]] + ni*sigmav2))/sigmae2[group[i]]
    Vinvlist[[i]] <- S1inv
    
  }
  
  Vinv <- bdiag(Vinvlist)
  
  betahat <- solve(t(x)%*%Vinv %*% x) %*% t(x) %*% Vinv %*% y
  
  loglikevalue <- -0.5*detV - 0.5*t(y - x %*% betahat) %*% Vinv %*% (y - x %*% betahat)
  
  
  covest <- solve(t(x)%*%Vinv %*% x) 
  se <- sqrt(diag(covest))
  
  return(list(est = betahat, 
              se = se,
              covest = covest,
              sig2est = sig2est,
              loglikevalue = loglikevalue))
}









