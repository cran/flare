\name{flare-internal}
\alias{sugm.likelihood}
\alias{sugm.tracel2}
\alias{sugm.cv}
\alias{part.cv}
\alias{sugm.clime.ladm.scr}
\alias{sugm.tiger.ladm.scr}
\alias{slim.lad.ladm.scr.btr}
\alias{slim.sqrt.ladm.scr}
\alias{slim.dantzig.ladm.scr}
\alias{slim.lq.ladm.scr.btr}
\alias{slim.lasso.ladm.scr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal flare functions}
\description{
  Internal flare functions
}
\usage{
sugm.likelihood(Sigma, Omega)
sugm.tracel2(Sigma, Omega)
sugm.cv(obj, loss=c("likelihood", "tracel2"), fold=5)
part.cv(n, fold)
sugm.clime.ladm.scr(Sigma, lambda, nlambda, n, d, maxdf, rho, shrink, prec, 
                    max.ite, verbose)
sugm.tiger.ladm.scr(data, n, d, maxdf, rho, lambda, shrink, prec, 
                    max.ite, verbose)
slim.lad.ladm.scr.btr(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, 
                      intercept, verbose)
slim.sqrt.ladm.scr(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, 
                   intercept, verbose)
slim.dantzig.ladm.scr(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, 
                      intercept, verbose)
slim.lq.ladm.scr.btr(Y, X, q, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, 
                     intercept, verbose)
slim.lasso.ladm.scr(Y, X, lambda, nlambda, n, d, maxdf, max.ite, prec, 
                    intercept, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{Covariance matrix.}
  \item{Omega}{Inverse covariance matrix.}
  \item{obj}{An object with S3 class returned from \code{"sugm"}.}
  \item{loss}{Type of loss function for cross validation.}
  \item{fold}{The number of fold for cross validatio.}
  \item{n}{The number of observations (sample size).}
  \item{d}{Dimension of data.}
  \item{maxdf}{Maximal degree of freedom.}
  \item{lambda}{Grid of non-negative values for the regularization parameter lambda.}
  \item{nlambda}{The number of the regularization parameter lambda.}
  \item{shrink}{Shrinkage of regularization parameter based on precision of estimation.}
  \item{rho}{Value of augmented Lagrangian multipiler.}
  \item{prec}{Stopping criterion.}
  \item{max.ite}{Maximal value of iterations.}
  \item{data}{\code{n} by \code{d} data matrix. }
  \item{Y}{Dependent variables in linear regression. }
  \item{X}{Design matrix in linear regression. }
  \item{q}{The vector norm used for the loss term. }
  \item{intercept}{The indicator of whether including intercepts specifically. }
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
}
\details{
  These are not intended for use by users.   
}

%\references{
%}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sugm}}, \code{\link{slim}} and \code{\link{flare-package}}.
}
