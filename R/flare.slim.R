#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim(): The user interface for slim()                                      #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Nov 23th, 2012                                                             #
# Version: 0.9.5                                                                   #
#----------------------------------------------------------------------------------#

flare.slim <- function(X, 
                       Y, 
                       lambda = NULL,
                       nlambda = NULL,
                       lambda.min.ratio = NULL,
                       rho = NULL,
                       method="lq",
                       q = 2,
                       prec = 1e-3,
                       max.ite = 1e3,
                       verbose = TRUE)
{
  if(method!="dantzig" && method!="lq"){
    cat("\"method\" must be dantzig, lasso or lq.\n")
    return(NULL)
  }
  if(method=="lq"){
    if(q<1 || q>2){
      cat("q must be in [1, 2] when method = \"lq\".\n")
      return(NULL)
    }
  }
  if(verbose) {
    cat("Sparse Linear Regression with L1 Regularization.\n")
  }
  n = nrow(X)
  d = ncol(X)
  maxdf = max(n,d)
#   if(intercept){
#     X = cbind(rep(1, n), X)
#     d = d+1
#   }
  xm=matrix(rep(.colMeans(X,n,d),n),nrow=n,ncol=d,byrow=T)
  x1=X-xm
  sdx=sqrt(diag(t(x1)%*%x1)/(n-1))
  Cxinv=diag(1/sdx)
  xx=x1%*%Cxinv
  ym=mean(Y)
  y1=Y-ym
  sdy=sqrt(sum(y1^2)/(n-1))
  yy=y1/sdy
  corr=cor(xx,yy)
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 5
    if(is.null(lambda.min.ratio)){
      if(method=="dantzig")
        lambda.min.ratio = 0.8
      else
        lambda.min.ratio = 0.25
    }
    if(method=="dantzig")
      lambda.max = max(corr)
    else 
      lambda.max = pi*sqrt(log(d)/n) #min(max(S-diag(diag(S))),-min(S-diag(diag(S))))
    lambda.min = lambda.min.ratio*lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    rm(lambda.max,lambda.min,lambda.min.ratio)
    gc()
  }
  if(is.null(rho))
    rho = sqrt(d)
  begt=Sys.time()
  if(method=="dantzig") # dantzig
    out = flare.slim.dantzig(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
  if(method=="lq" && q>=1) {#  && q<1e5 && q!=2 && q!="lasso"
    if(is.null(q)) q=2;
    if(q==1) # lad lasso
      out = flare.slim.lad(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
    if(q==2) # sqrt lasso
      out = flare.slim.sqrt(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
    if(q>1 && q<2)
      out = flare.slim.lq(yy, xx, q, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
  }
  runt=Sys.time()-begt
  
  sparsity=rep(0,nlambda)
  for(i in 1:nlambda)
    sparsity[i] = sum(out$beta[[i]]!=0)/d
    
  est = list()    
  beta1=matrix(0,nrow=d,ncol=nlambda)
  intercept=matrix(0,nrow=1,ncol=nlambda)
  for(k in 1:nlambda){
    intercept[k]=ym-xm[1,]%*%Cxinv%*%out$beta[[k]]*sdy
    beta1[,k]=Cxinv%*%out$beta[[k]]*sdy
  }
  est$beta = beta1
  est$intercept = intercept
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$sparsity = sparsity
  est$method = method
  est$q = q
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  return(est)
}

print.slim <- function(x, ...)
{  
  cat("\n slim options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method=", x$method, "\n")
  if(x$method=="lq"){
    if(x$q==1){
      cat("q=",x$q," loss, LAD Lasso\n")
    } else {
      if(x$q==2)
        cat("q=",x$q," loss, SQRT Lasso\n")
      else
        cat("q=",x$q," loss\n")
    }
  }
  cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
  if(units.difftime(x$runtime)=="secs") unit="secs"
  if(units.difftime(x$runtime)=="mins") unit="mins"
  if(units.difftime(x$runtime)=="hours") unit="hours"
  cat("Runtime:",x$runtime," ",unit,"\n")
}

plot.slim <- function(x, ...)
{
  matplot(x$lambda, t(x$beta), type="l", main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient")
}
