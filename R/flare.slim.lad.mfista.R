#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim.lad.mfista(): Regression with LAD Lasso()                             #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Nov 5th, 2012                                                              #
# Version: 0.9.7                                                                   #
#----------------------------------------------------------------------------------#

flare.slim.lad.mfista <- function(Y, X, lambda, nlambda, n, d, maxdf, mu, max.ite, prec,intercept, verbose)
{
  if(verbose==TRUE)
    cat("LAD Lasso regression via MFISTA.\n")
  XX = t(X)%*%X
  L = eigen(XX)$values[1]
  beta = array(0,dim=c(d,nlambda))
  ite.ext = rep(0,nlambda)
  obj = array(0,dim=c(max.ite,nlambda))
  runt = array(0,dim=c(max.ite,nlambda))
  if(intercept) intercept=1
  else intercept=0
  str=.C("slim_lad_mfista", as.double(Y), as.double(X), as.double(XX), 
         as.double(beta), as.integer(n), as.integer(d), as.double(mu),
         as.integer(ite.ext), as.double(lambda), as.integer(nlambda), 
         as.integer(max.ite), as.double(prec), as.double(L), 
         as.double(obj),as.double(runt),as.integer(intercept),PACKAGE="flare")
  
  beta.list = vector("list", nlambda)
  obj.list = vector("list", nlambda)
  runt.list = vector("list", nlambda)
  ite.ext = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
    obj.i = unlist(str[14])[((i-1)*max.ite+1):((i-1)*max.ite+ite.ext[i])]
    obj.list[[i]] = obj.i
    runt.i = unlist(str[15])[((i-1)*max.ite+1):((i-1)*max.ite+ite.ext[i])]
    runt.list[[i]] = runt.i
  }
  
  return(list(beta=beta.list, ite=ite.ext, obj=obj.list,runt=runt.list))
}
