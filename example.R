#### examples #####

library(flexclust)
library(cluster)
library(lme4)
library(flexmix)

source("CluRER.R")
source("refitLMM_RE.R")




### simulated data
n <- 15

betab <- c(3,0.2,-0.2)
sdu <- 0.18
tauv <- c(-0.7,0.5)
probvalue <- c(4/5,1/5) 

p <- length(betab)-1
gammat <- exp(tauv)
sde <- sdu/sqrt(gammat)
sdes <- 0.05


set.seed(1)
m <- sample(30:85, 15)
gammat <- sdu^2/sde^2
id_info <- data.frame(id = 1:n,
                      group = c(3,sample(1:2, size = n-1,prob = probvalue,
                                         replace = TRUE)))

ylist <- list()
xlist <- list()

u <- rnorm(n)*sdu

xmati <- cbind(1,matrix(rnorm(m[1]*p),ncol = p))
yi <- xmati %*% betab + u[1] + sdes*rnorm(m[1])
xlist[[1]] <- xmati
ylist[[1]] <- cbind(id = id_info$id[1],y = yi)


for(i in 2:n)
{
  xmati <- cbind(1,matrix(rnorm(m[i]*p),ncol = p))
  yi <- xmati %*% betab + u[i] + sde[id_info$group[i]]*rnorm(m[i])
  xlist[[i]] <- xmati
  ylist[[i]] <- cbind(id = id_info$id[i],y = yi)
}

xmat <- do.call("rbind",xlist)
ydf <- as.data.frame(do.call("rbind",ylist))
colnames(ydf) <- c("id","y")

## end of simulated data

# initial values
reslmm <- lmer(y~xmat[,-1] + (1|id), data = ydf)
residmm <- residuals(reslmm)
vare0 <- as.numeric(by(residmm, ydf$id, FUN = var))
varu0 <- as.numeric(as.data.frame(VarCorr(reslmm))["vcov"][1,])
gamma0 <- varu0/vare0

lamvec <- c(exp(seq(-1,-0.65,length = 50)),
            exp(seq(-0.65,-0.6,length = 50)),
            exp(seq(-0.6,1, length = 50)))
nlam <- length(lamvec)
nll2 <- ngest <- rep(0, nlam)
groupest <- matrix(0, n, nlam)
sigv2est <- rep(0, nlam)
sige2est <- matrix(0,n,nlam)
for(j in 1:length(lamvec))
{
  resj <- CluRER(indexy = ydf$id, y = ydf$y, x = xmat,betam0 = betam0,
                 gamma0 = gamma0,varu0 = varu0,
                 lambda = lamvec[j],maxiter = 500)
  
  refitj <- refitLMM_RE(indexy = ydf$id, y = ydf$y, x = xmat, group = resj$groupest)
  nll2[j] <- -refitj$loglikevalue
  ngest[j] <- length(unique(resj$groupest))
  groupest[,j] <- resj$groupest
  sigv2est[j] <- resj$sigmav2
  sige2est[,j] <- resj$sigmae2
}

lamdf <- data.frame(lam = lamvec, nll = nll2, ng = ngest)
lamdf <- lamdf %>% mutate(bic = 2*nll + log(log(n))*log(sum(m))*ng)

index1 <- which.min(lamdf$bic)

ari <- randIndex(groupest[,index1], id_info$group)
ngest <- length(unique(groupest[,index1]))
ari


# kmeans based on initial two Gaps ####
gapsim <- clusGap(matrix(gamma0,ncol=1),FUNcluster = kmeans,K.max=10, B=1000)
gapres <- maxSE(gapsim$Tab[,3],gapsim$Tab[,4],method = "globalSEmax")
ari_kmeans_gap <- randIndex(kmeans(gamma0, gapres, iter.max = 1000, nstart = 100)$cluster,id_info$group)
ari_kmeans_gap


gapres <- maxSE(gapsim$Tab[,3],gapsim$Tab[,4],method = "Tibs2001SEmax")
ari_kmeans_gap <- randIndex(kmeans(gamma0, gapres, iter.max = 1000, nstart = 100)$cluster,id_info$group)
ari_kmeans_gap


# kmeans based on initial log two Gaps ####
gapsim <- clusGap(matrix(log(gamma0),ncol=1),FUNcluster = kmeans,K.max=10, B=1000)
gapres_log <- maxSE(gapsim$Tab[,3],gapsim$Tab[,4],method = "globalSEmax")
ari_kmeans_gap_log <- randIndex(kmeans(log(gamma0), gapres_log, iter.max = 1000, nstart = 100)$cluster,id_info$group)
ari_kmeans_gap_log



gapres_log <- maxSE(gapsim$Tab[,3],gapsim$Tab[,4],method = "Tibs2001SEmax")
ari_kmeans_gap_log <- randIndex(kmeans(log(gamma0), gapres_log, iter.max = 1000, nstart = 100)$cluster,id_info$group)
ari_kmeans_gap_log


# flexmix  ####

resfb <- stepFlexmix(y~xmat[,-1]|id,k = 1:8,nrep = 20,verbose = FALSE, data = ydf)
resf <- getModel(resfb, "BIC")
fgroup <- unique(cbind(ydf$id, flexmix::clusters(resf)))
arif <- randIndex(fgroup[,2], id_info$group)
ngf <- length(unique(fgroup[,2]))
arif





