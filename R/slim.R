#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim(): The user interface for slim()                                      #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Mar 16th 2014                                                              #
# Version: 1.2.0                                                                   #
#----------------------------------------------------------------------------------#

slim <- function(X, 
                 Y, 
                 lambda = NULL,
                 nlambda = NULL,
                 lambda.min.value = NULL,
                 lambda.min.ratio = NULL,
                 rho = 1,
                 method="lq",
                 q = 2,
                 res.sd = FALSE,
                 prec = 1e-5,
                 max.ite = 1e5,
                 verbose = TRUE)
{
  if(method!="dantzig" && method!="lq" && method!="lasso"){
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
  if(is.null(X)||is.null(Y)) {
    cat("No data input.\n")
    return(NULL)
  }
  X = as.matrix(X)
  Y = as.matrix(Y)
  if(ncol(Y) != 1) {
    cat("Y must be a vector or one-column matrix.\n")
    return(NULL)
  }
  if(nrow(X) != nrow(Y)) {
    cat("X and Y must have the same number of rows.\n")
    return(NULL)
  }
  if(anyNA(X) || anyNA(Y)){
    cat("The input has missing values.\n")
    Xrow.na = rowSums(is.na(X)) > 0
    Yrow.na = rowSums(is.na(Y)) > 0
    idx.na = which(Xrow.na | Yrow.na)
    X = X[-idx.na,,drop=FALSE]
    Y = Y[-idx.na,,drop=FALSE]
    n = nrow(X)
    if(n==0) {
      cat("Too many missing values.\n")
      return(NULL)
    }
  }
  n = nrow(X)
  d = ncol(X)
  if(n==0 || d==0) {
    cat("No data input.\n")
    return(NULL)
  }
  if(n <= 1) {
    cat("At least two observations are required.\n")
    return(NULL)
  }
  maxdf = max(n,d)
  x.mean = colMeans(X)
  x1 = sweep(X, 2, x.mean, FUN = "-")
  x.sd = sqrt(colSums(x1^2)/(n-1))
  if(any(x.sd <= 0)) {
    cat("X contains at least one constant column.\n")
    return(NULL)
  }
  sdxinv = 1/x.sd
  xx = sweep(x1, 2, sdxinv, FUN = "*")
  ym=mean(Y)
  y1=Y-ym
  if(res.sd == TRUE){
    sdy=sqrt(sum(y1^2)/(n-1))
    if(sdy <= 0){
      cat("Y has zero variance after centering.\n")
      return(NULL)
    }
    yy=y1/sdy
  }else{
    sdy = 1
    yy = y1
  }
  intercept = FALSE
  
  if(intercept){
    xx = cbind(rep(1, nrow(xx)), xx)
    X = cbind(rep(1, nrow(X)), X)
    d = d+1
  }
  
  if(!is.null(lambda)) {
    lambda = as.double(lambda)
    nlambda = length(lambda)
  }
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 5
    if(method=="dantzig"){
      if(intercept)
        lambda.max = max(abs(crossprod(xx[,2:d],yy/n)))
      else
        lambda.max = max(abs(crossprod(xx,yy/n)))
    }
    if(method=="lq"){
      if(q==2){
        if(intercept)
          lambda.max = max(abs(crossprod(xx[,2:d],yy/sqrt(sum(yy^2))/sqrt(n))))
        else
          lambda.max = max(abs(crossprod(xx,yy/sqrt(sum(yy^2))/sqrt(n))))
      }else{
        if(q==1){
          if(intercept)
            lambda.max = max(abs(crossprod(xx[,2:d],sign(yy)/n)))
          else
            lambda.max = max(abs(crossprod(xx,sign(yy)/n)))
        }else{
          if(intercept){
            lambda.max = max(abs(crossprod(xx[,2:d],sign(yy)*(abs(yy)^(q-1))/(sum(abs(yy)^q)^((q-1)/q))/n^(1/q))))# 1<=q<=2
          }else{
            lambda.max = max(abs(crossprod(xx,sign(yy)*(abs(yy)^(q-1))/(sum(abs(yy)^q)^((q-1)/q))/n^(1/q)))) # 1<=q<=2
          }
        }
      }
    }
    if(method=="lasso"){
      if(intercept)
        lambda.max = max(abs(crossprod(xx[,2:d],yy/n)))
      else
        lambda.max = max(abs(crossprod(xx,yy/n)))
    }
    if(method=="dantzig"){
      if(is.null(lambda.min.ratio)){
        lambda.min.ratio = 0.5
      }
      if(is.null(lambda.min.value)){
        lambda.min.value = lambda.min.ratio*lambda.max
      }
    }else{
      if(is.null(lambda.min.value)){
        lambda.min.value = sqrt(log(d)/n)
      }else{
        if(is.null(lambda.min.ratio)){
          lambda.min.ratio = lambda.min.value/lambda.max
        }
      }
    }
    if(lambda.max<lambda.min.value){
      lambda.max = 1
      lambda.min.value = 0.4
    }
    lambda = exp(seq(log(lambda.max), log(lambda.min.value), length.out = nlambda))
    rm(lambda.max,lambda.min.value,lambda.min.ratio)
  }
  nlambda = length(lambda)
  if(nlambda == 0){
    cat("At least one lambda value is required.\n")
    return(NULL)
  }
  if(any(!is.finite(lambda)) || any(lambda <= 0)){
    cat("lambda must contain positive finite values.\n")
    return(NULL)
  }
  if(is.null(rho))
    rho = 1
  begt=Sys.time()
  if(method=="dantzig"){ # dantzig
    if(d>=n)
      out = slim.dantzig.ladm.scr(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    else
      out = slim.dantzig.ladm.scr2(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    q = "infty"
  }
  if(method=="lq") {#  && q!=2 && q!="lasso"
    if(q==1) # lad lasso
      out = slim.lad.ladm.scr.btr(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    if(q==2) # sqrt lasso
      out = slim.sqrt.ladm.scr(yy, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    if(q>1 && q<2) # lq lasso
      out = slim.lq.ladm.scr.btr(yy, xx, q, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
  }
  if(method=="lasso")
    out = slim.lasso.ladm.scr(yy, xx, lambda, nlambda, n, d, maxdf, max.ite, prec, intercept, verbose)
  runt=Sys.time()-begt
  
  df=rep(0,nlambda)
  if(intercept){
    for(i in seq_len(nlambda))
      df[i] = sum(out$beta[[i]][2:d]!=0)
  }else{
    for(i in seq_len(nlambda))
      df[i] = sum(out$beta[[i]]!=0)
  }
  
  est = list()
  intcpt0 = numeric(nlambda)
  intcpt = numeric(nlambda)
  if(intercept){
    beta1=matrix(0,nrow=d-1,ncol=nlambda)
    for(k in seq_len(nlambda)){
      tmp.beta = out$beta[[k]][2:d]
      beta1[,k]=sdxinv*tmp.beta*sdy
      intcpt[k] = ym-as.numeric(x.mean%*%beta1[,k])+out$beta[[k]][1]*sdy
      intcpt0[k] = intcpt[k]
    }
  }else{
    beta1=matrix(0,nrow=d,ncol=nlambda)
    for(k in seq_len(nlambda)){
      tmp.beta = out$beta[[k]]
      intcpt0[k] = 0
      beta1[,k] = sdxinv*tmp.beta*sdy
      intcpt[k] = ym-as.numeric(x.mean%*%beta1[,k])
    }
  }
  
  est$beta0 = out$beta
  est$beta = beta1
  est$intercept = matrix(intcpt, nrow=1)
  est$intercept0 = matrix(intcpt0, nrow=1)
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$df = df
  est$method = method
  est$q = q
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  if(verbose) print(est)
  return(est)
}

print.slim <- function(x, ...)
{  
  cat("\n")
  cat("slim options summary: \n")
  cat(x$nlambda, "lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method =", x$method, "\n")
  if(x$method=="lq"){
    if(x$q==1){
      cat("q =",x$q," loss, LAD Lasso\n")
    } else {
      if(x$q==2)
        cat("q =",x$q,"loss, SQRT Lasso\n")
      else
        cat("q =",x$q,"loss\n")
    }
  }
  cat("Degree of freedom:",min(x$df),"----->",max(x$df),"\n")
  unit = as.character(units(x$runtime))
  cat("Runtime:",x$runtime,unit,"\n")
}

plot.slim <- function(x, ...)
{
  matplot(x$lambda, t(x$beta), type="l", main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient")
}

coef.slim <- function(object, lambda.idx = c(1:3), beta.idx = c(1:3), ...)
{
  lambda.n = length(lambda.idx)
  beta.n = length(beta.idx)
  cat("\n Values of estimated coefficients: \n")
  cat(" index     ")
  for(i in seq_len(lambda.n)){
    cat("",formatC(lambda.idx[i],digits=5,width=10),"")
  }
  cat("\n")
  cat(" lambda    ")
  for(i in seq_len(lambda.n)){
    cat("",formatC(object$lambda[lambda.idx[i]],digits=4,width=10),"")
  }
  cat("\n")
  cat(" intercept ")
  for(i in seq_len(lambda.n)){
    cat("",formatC(object$intercept[lambda.idx[i]],digits=4,width=10),"")
  }
  cat("\n")
  for(i in seq_len(beta.n)){
    cat(" beta",formatC(beta.idx[i],digits=5,width=-5))
    for(j in seq_len(lambda.n)){
      cat("",formatC(object$beta[beta.idx[i],lambda.idx[j]],digits=4,width=10),"")
    }
    cat("\n")
  }
}
predict.slim <- function(object, newdata, lambda.idx = c(1:3), Y.pred.idx = c(1:5), ...)
{
  newdata = as.matrix(newdata)
  if(ncol(newdata) != nrow(object$beta)){
    cat("newdata has incompatible column size.\n")
    return(NULL)
  }
  pred.n = nrow(newdata)
  lambda.n = length(lambda.idx)
  Y.pred.n = length(Y.pred.idx)
  intcpt = matrix(object$intercept[,lambda.idx], nrow=pred.n,
                  ncol=lambda.n,byrow=TRUE)
  Y.pred = newdata%*%object$beta[,lambda.idx] + intcpt
  cat("\n Values of predicted responses: \n")
  cat("   index   ")
  for(i in seq_len(lambda.n)){
    cat("",formatC(lambda.idx[i],digits=5,width=10),"")
  }
  cat("\n")
  cat("   lambda  ")
  for(i in seq_len(lambda.n)){
    cat("",formatC(object$lambda[lambda.idx[i]],digits=4,width=10),"")
  }
  cat("\n")
  for(i in seq_len(Y.pred.n)){
    cat("    Y",formatC(Y.pred.idx[i],digits=5,width=-5))
    for(j in seq_len(lambda.n)){
      cat("",formatC(Y.pred[Y.pred.idx[i],j],digits=4,width=10),"")
    }
    cat("\n")
  }
  return(list(Y.pred = Y.pred))
}
