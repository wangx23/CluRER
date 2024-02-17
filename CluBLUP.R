#### calculate BLUP based on group and estimate information#####
### sig2est first is sigmav, later are sigmae for different groups ###
CluBLUP <- function(indexy, y, x, group, betahat, sig2est )
{
  Zm <- model.matrix(~ -1 + as.factor(indexy)) 
  nind <- length(unique(indexy))
  ngest <- length(unique(group))
  
  sig2v <- sig2est[1]
  sig2e <- sig2est[2:length(sig2est)]
  
  Gm <- sig2v*diag(nind)
  Sigm <- Zm %*%Gm %*% t(Zm)
  
  Rv <- sig2e[rep(group,as.numeric(table(indexy)))]
  diag(Sigm) <- diag(Sigm) + Rv
  
  rande <- Gm %*% t(Zm) %*% solve(Sigm) %*%(y - x %*% betahat)
  resid_value <- as.numeric(y - x %*% betahat - rep(rande,  as.numeric(table(indexy))))
  
  return(list(rande = rande, residuals = resid_value))
}