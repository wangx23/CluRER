library(igraph)
library(Matrix)
library(dplyr)
CluRER <- function(indexy, y, x, betam0, gamma0, varu0,lambda = 0.1,
                   maxiter = 200, nu =1, tolabs = 1e-4, tolrel = 1e-2,
                   tolr = 1e-3)
{
  nt <- length(y)
  np <- ncol(x)
  
  uindex <- unique(indexy)
  n <- length(uindex)
  
  
  npair <- n*(n-1)/2
  
  nJ <- as.numeric(table(indexy))
  
  
  ## the A matrix ####
  
  fullg <- make_full_graph(n)
  
  linkMatrix <- get.data.frame(fullg) %>%
    dplyr::mutate(edge = row_number()) %>% dplyr::select(edge, from, to)
  linkMatrix$from <- linkMatrix$from %>% as.integer()
  linkMatrix$to <- linkMatrix$to %>% as.integer()
  Dmat <- as(matrix(0, npair, n), "dgCMatrix")
  Dmat[as.matrix(linkMatrix[, c(1, 2)])] = 1
  Dmat[as.matrix(linkMatrix[, c(1, 3)])] = -1
  
  
  
  #### initial values
  sigmav2 <- varu0
  sigmae2 <- sigmav2/gamma0
  tau0 <- log(gamma0)
  deltam0 <- Dmat %*% tau0
  
  vm <- rep(0, npair)
  
  deltam <- deltam0
  deltamold <- as.numeric(deltam0)
  tauv <- tau0
  gammav <- gamma0
  
  for(m in 1:maxiter)
  {
    
    ### update beta 
    
    # the V matrix ### (covariance)
    
    Vinvlist <- list()
    
    for(i in 1:n)
    {
      indexi <- indexy == uindex[i]
      ni <- nJ[i]
      mat1 <- matrix(1,ni,ni)
      diag1 <- diag(1,ni)
      
      # S1inv <- (diag1 - mat1*sigmav2/(sigmae2 + ni*sigmav2))/sigmae2
      S1inv <- gammav[i]*(diag1 - mat1*sigmav2/(sigmae2[i] + ni*sigmav2))
      Vinvlist[[i]] <- S1inv
      
    }
    
    Vinv <- bdiag(Vinvlist)
    
    betam <- as.numeric(solve(t(x) %*% Vinv %*% x)%*% t(x) %*% Vinv %*% y)
    
    
    ### update sigmav2 ####
    resid <- y - x %*% betam
    residu <- as.numeric(by(resid, indexy, sum))
    resid2 <- as.numeric(by(resid^2, indexy, sum))
    
    #(sum(gammav * resid2) - sum(gammav^2/(1+nJ*gammav)*residu^2))/nt
    
    sigmav2_new <- as.numeric(t(resid) %*% Vinv %*% (resid)/nt)
    
    ### update gammav ###
    
    
    tau_fun <- function(x)
    {
      gx <- exp(x)
      value11 <- -sum(nJ*log(gx)) + sum(log(1+nJ*gx)) 
      value12 <- sum(resid2*gx)/sigmav2_new - sum(residu^2*gx^2/(1+nJ*gx))/sigmav2_new
      value2 <- nu *sum((Dmat %*% x  - deltam  + vm/nu)^2)
      value11 + value12 + value2
    }
    
    
    
    grad_fun <- function(x)
    {
      gx <- exp(x)
      value11 <- -nJ/gx + nJ/(1+nJ*gx)
      value12 <- resid2/sigmav2_new - residu^2*(2*gx + nJ*gx^2)/(1+nJ*gx)^2/sigmav2_new
      value2 <- 2*nu * t(Dmat) %*% (Dmat %*% x  - deltam  + vm/nu)
      as.numeric(value11*exp(x) + value12*exp(x)+value2)
    }
    
    
    
    tau_new <- optim(tau_fun, par = tauv, gr = grad_fun, method ="BFGS", control = list(maxit =1000))$par
    
    
    ### update deltam 
    
    tau_diff = as.numeric(Dmat %*% tau_new)
    
    deltam <- tau_diff + (1/nu) * vm
    
    deltam <- unlist(lapply(deltam, FUN = mcp, lam = lambda, gam = 3, nu = 1))
    
    ### update v ###
    
    vm <-  vm + nu * (tau_diff - deltam)
    
    gamma_new <- exp(tau_new)
    
    sigmae2 <- sigmav2_new/gamma_new
    sigmav2 <- sigmav2_new
    
    tauv <- tau_new
    gammav <- gamma_new
    
    
    norm1 = norm(tau_diff,type = "2")
    norm2 = norm(deltam,type = "2")
    
    tolpri = tolabs*sqrt(npair) + tolrel*max(norm1, norm2);
    toldual = tolabs*sqrt(n) + tolrel * norm(vm %*% Dmat, type = "2")
    
    rm = norm(tau_diff- deltam,type = "2")
    sm = nu * norm((deltam - deltamold)%*%Dmat, "2")
    
    deltamold <- deltam
    
    if(rm <= tolr){
      break
    }
  }
  
  niteration <- m
  
  groupest <- getgroup(matrix(deltam, nrow = 1),n)
  
  #### get the loglikelihood 
  
  loglikevalue <- 0
  detV <- 0
  detV1 <- 0
  
  resid <- y - x %*% betam
  
  
  for(i in 1:n)
  {
    indexi <- indexy == uindex[i]
    ni <- nJ[i]
    mat1 <- matrix(1,ni,ni)
    diag1 <- diag(1,ni)
    
    S1inv <- (diag1 - mat1*sigmav2/(sigmae2[i] + ni*sigmav2))/sigmae2[i]
    
    residi <- resid[indexi]
    
    loglikevalue <- loglikevalue + t(residi) %*% S1inv %*%(residi)
    
    
  }
  
  
  
  
  loglikevalue <- - 0.5*(sum(nJ*log(sigmae2)) + sum(log(1+nJ*gammav))) - 0.5*loglikevalue
  
  
  outls <- list(betaest = betam,
                sigmae2 = sigmae2,
                sigmav2 = sigmav2,
                gammav = gammav,
                tauv = tauv,
                loglikevalue = loglikevalue,
                groupest = groupest,
                niteration = niteration,
                rm = rm)
  return(outls)
  
  
}





#####  mcp function #####


sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  if(xn==0){
    return(x)
  }else{
    thval <- 1 - th/xn
    return(thval*((thval) >0)*x)
  }
}
mcp <- function(x,lam,gam,nu)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  return(z)
}

#### getgroup 

getgroup = function(deltam, n, tol = 1e-2)
{
  p = nrow(deltam)
  b2value =sqrt(colMeans(deltam^2))
  b2value[b2value <= tol] = 0
  
  d2 = matrix(0, n, n)
  for(j in 1:(n-1))
  {
    indexj1 = (2*n -j)*(j-1)/2 + 1
    indexj2 = indexj1 + n - j - 1
    d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
  }
  d2 = t(d2) + d2
  
  
  ngj = 1:n
  groupest = rep(0,n)
  j = 0
  
  while (length(ngj) >0)
  {
    j = j + 1
    gj = (1:n)[d2[ngj[1],] ==0]
    indexj = ngj %in% gj
    gj = ngj[indexj]
    ngj = ngj[!indexj]
    groupest[gj] = j * rep(1,length(gj))
  }
  
  
  return(groupest)
  
}




