#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim.lq(): Regression with Lq norm Lasso()                                 #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Nov 5th, 2012                                                              #
# Version: 0.9.7                                                                   #
#----------------------------------------------------------------------------------#

flare.slim.lq <- function(Y, X, q, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept,verbose)
{
  if(verbose==TRUE)
    cat("LQ norm regrelarization regression (q =", q, ") .\n")
  XX = t(X)%*%X
  beta = array(0,dim=c(d,nlambda))
  ite.ext = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  ite.int2 = rep(0,nlambda)
  obj = array(0,dim=c(max.ite,nlambda))
  runt = array(0,dim=c(max.ite,nlambda))
  if(intercept) intercept=1
  else intercept=0
  str=.C("slim_lq", as.double(Y), as.double(X), as.double(XX), 
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.ext), as.integer(ite.int1), as.integer(ite.int2), 
         as.double(q), as.double(lambda), as.integer(nlambda), as.integer(max.ite), 
         as.double(prec), as.integer(intercept),as.double(obj),as.double(runt),
         PACKAGE="flare")
  ite.ext = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.int1 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite.int2 = matrix(unlist(str[10]), byrow = FALSE, ncol = nlambda)
  ite = vector("list", 3)
  ite[[1]] = ite.ext
  ite[[2]] = ite.int1
  ite[[3]] = ite.int2
  beta.list = vector("list", nlambda)
  obj.list = vector("list", nlambda)
  runt.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
    obj.i = unlist(str[17])[((i-1)*max.ite+1):((i-1)*max.ite+ite.ext[i])]
    obj.list[[i]] = obj.i
    runt.i = unlist(str[18])[((i-1)*max.ite+1):((i-1)*max.ite+ite.ext[i])]
    runt.list[[i]] = runt.i
  }
  return(list(beta=beta.list, ite=ite, obj=obj.list,runt=runt.list))
}
