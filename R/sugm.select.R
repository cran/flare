#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# sugm.select(): Model selection using:                                            #
#                     1.cross validation (cv)                                      #
#                     2.stability approach to regularization selection (stars)     #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Dec 2nd 2013                                                               #
# Version: 1.1.0                                                                   #
#----------------------------------------------------------------------------------#

## Main Function
sugm.select <- function(est, 
                        criterion = "stars", 
                        stars.subsample.ratio = NULL,
                        stars.thresh = 0.1, 
                        rep.num = 20,
                        fold = 5,
                        loss="likelihood", 
                        verbose = TRUE)
{
  if(est$method!="clime" && est$method!="tiger") {
    cat("\"method\" must be either \"clime\" or \"tiger\" \n")
    return(NULL)
  }
  gcinfo(FALSE)
  
  if(est$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
    class(est) = "select"
    return(est)
  }
  
  if(!est$cov.input)
  {
    if(is.null(criterion))
      criterion = "stars" # stars cv
    criterion = match.arg(criterion, c("stars", "cv"))
    if(rep.num < 1)
      stop("rep.num must be >= 1.")
    if(stars.thresh < 0 || stars.thresh > 1)
      stop("stars.thresh must be in [0, 1].")
    n = nrow(est$data)
    d = ncol(est$data)
    nlambda = length(est$lambda)
    
    if(criterion == "cv"){
      if(verbose)
      {
        cat("Conducting cross validation (cv) selection....\n")
        flush.console()
      }
      
      out = sugm.cv(est, loss=loss, fold = fold)
      est$opt.lambda = out$lambda.opt
      est$opt.index = out$opt.idx
      est$loss = out$loss
      rm(out)
      gc()
      
      if(verbose){
        cat("done\n")
        flush.console()
      }
      
      if(verbose)
      {
        cat("Computing the optimal graph....\n")
        flush.console()
      }
      
      est$refit = est$path[[est$opt.index]]
      #est$refit[abs(est$icov[[est$opt.index]])<5e-2]=0
      est$opt.sparsity=sum(est$refit)/d/(d-1)
      est$opt.icov = est$icov[[est$opt.index]]
      
      if(verbose){
        cat("done\n")
        flush.console()
      }
    }
    
    if(criterion == "stars"){
      if(is.null(stars.subsample.ratio))
      {
        if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
        if(n<=144) stars.subsample.ratio = 0.8
      } 
      if(stars.subsample.ratio <= 0 || stars.subsample.ratio >= 1)
        stop("stars.subsample.ratio must be in (0, 1).")
      
      est$merge = list()
      for(i in seq_len(nlambda)) est$merge[[i]] = Matrix(0,d,d)
      
      for(rep.idx in seq_len(rep.num))
      {
        if(verbose)
        {
          mes <- paste(c("Conducting Subsampling....in progress:", floor(100*rep.idx/rep.num), "%"), collapse="")
          cat(mes, "\r")
          flush.console()	
        }
        ind.sample = sample(seq_len(n), floor(n*stars.subsample.ratio), replace=FALSE)
        
        if(est$method == "clime")
          tmp = sugm(est$data[ind.sample,,drop=FALSE], lambda = est$lambda, method = "clime", sym = est$sym, verbose = FALSE,
                      standardize=est$standardize)$path
        if(est$method == "tiger")
          tmp = sugm(est$data[ind.sample,,drop=FALSE], lambda = est$lambda, method = "tiger", sym = est$sym, verbose = FALSE,
                     standardize=est$standardize)$path
        
        for(lambda.idx in seq_len(nlambda))
          est$merge[[lambda.idx]] = est$merge[[lambda.idx]] + tmp[[lambda.idx]]
        
        rm(ind.sample,tmp)
        gc()
      }
      
      if(verbose){
        cat("\n")
        mes = "Conducting Subsampling....done. "
        cat(mes, "\r")
        cat("\n")
        flush.console()
      }
      
      est$variability = rep(0,nlambda)
      for(i in seq_len(nlambda)){
        est$merge[[i]] = est$merge[[i]]/rep.num
        est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
      }
      
      stars.cross = which(est$variability >= stars.thresh)
      if(length(stars.cross) == 0){
        est$opt.index = nlambda
      }else{
        est$opt.index = max(stars.cross[1]-1,1)
      }
      est$refit = est$path[[est$opt.index]]
      est$opt.lambda = est$lambda[est$opt.index]
      est$opt.sparsity = est$sparsity[est$opt.index]
      est$opt.icov = est$icov[[est$opt.index]]

    }
    est$criterion = criterion
    class(est) = "select"
    return(est)
  }
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

print.select = function(x, ...)
{
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.\n")
    return(NULL)
  }
  if(!x$cov.input)
  {
    if(x$method == "clime")
      cat("Method: CLIME\n")
    else
      cat("Method: SLasso\n")
    cat("selection criterion:",x$criterion,"\n")
    cat("Graph dimension:",ncol(x$data),"\n")
    cat("sparsity level:", x$opt.sparsity,"\n")
    cat("optimal paramter:", x$opt.lambda,"\n")
  }
  if(x$criterion == "cv")
    cat("cross validation loss used:",x$loss)
}

plot.select = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.\n")
    return(NULL)
  }
  if(!x$cov.input)
  {
    par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    g = graph_from_adjacency_matrix(as.matrix(x$refit), mode="undirected", diag=FALSE)
    layout.grid = layout_with_fr(g)
    
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA)	  
    plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    lines(x$opt.lambda,x$opt.sparsity,type = "p")  
  }
}
