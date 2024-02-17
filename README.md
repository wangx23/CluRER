# CluRER
Clustered variances in LMM

find clusters of variances in linear mixed models,
$y_{ih}=x_{ih}^{T}\beta+v_{i}+\epsilon_{ih},$
  where $v_{i} ~ {\sim}N\left(0,\sigma_{v}^{2}\right)$ is the random effect of lab $i$,
$\sigma_{i,\epsilon}^2$, that is, $\epsilon_{ih} \overset{iid}{\sim} N\left(0,\sigma_{i,\epsilon}^{2}\right)$ for $h=1,\dots,m_{i}$ and $i=1,\dots,n$. 

The goal is to find clusters of $\sigma_{i,\epsilon}$. 

