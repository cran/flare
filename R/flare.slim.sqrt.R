#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim.sqrt(): Regression with Square Root Lasso()                           #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Nov 5th, 2012                                                              #
# Version: 0.9.7                                                                   #
#----------------------------------------------------------------------------------#

flare.slim.sqrt <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
{
  cat("SQRT Lasso regression.\n")
  XX = t(X)%*%X
  beta = array(0,dim=c(d,nlambda))
  ite.ext = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  ite.int2 = rep(0,nlambda)
  beta.list = vector("list", nlambda)
  intercept=0
  str=.C("slim_sqrt", as.double(Y), as.double(X), as.double(XX), 
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.ext), as.integer(ite.int1), as.integer(ite.int2), 
         as.double(lambda), as.integer(nlambda), as.integer(max.ite), 
         as.double(prec), as.integer(intercept),PACKAGE="flare")
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  ite.ext = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.int1 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite.int2 = matrix(unlist(str[10]), byrow = FALSE, ncol = nlambda)
  ite = vector("list", 3)
  ite[[1]] = ite.ext
  ite[[2]] = ite.int1
  ite[[3]] = ite.int2
  return(list(beta=beta.list, ite=ite))
}
