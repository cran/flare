\name{flare-internal}
\alias{tiger.likelihood}
\alias{tiger.tracel2}
\alias{flare.tiger.cv}
\alias{part.cv}
\alias{flare.tiger.clime.cdadm}
\alias{flare.tiger.clime.ladm}
\alias{flare.tiger.clime.hadm}
\alias{flare.tiger.aclime.cdadm}
\alias{flare.tiger.slasso.cdadm}
\alias{flare.slim.lad}
\alias{flare.slim.sqrt}
\alias{flare.slim.dantzig}
\alias{flare.slim.lq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal flare functions}
\description{
  Internal flare functions
}
\usage{
tiger.likelihood(Sigma, Omega)
tiger.tracel2(Sigma, Omega)
flare.tiger.cv(obj, loss=c("likelihood", "tracel2"), fold=5)
part.cv(n, fold)
flare.tiger.clime.cdadm(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite)
flare.tiger.clime.ladm(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite)
flare.tiger.clime.hadm(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite)
flare.tiger.aclime.cdadm(Sigma, wmat, n, d, maxdf, lambda, rho, shrink, prec, max.ite)
flare.tiger.slasso.cdadm(data, n, d, maxdf, lambda, shrink, prec, max.ite)
flare.slim.lad(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
flare.slim.sqrt(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
flare.slim.dantzig(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
flare.slim.lq(Y, X, q, lambda, nlambda, n, d, maxdf, rho, max.ite, prec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{Covariance matrix.}
  \item{Omega}{Inverse covariance matrix.}
  \item{obj}{An object with S3 class returned from \code{"tiger"}.}
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
  \item{wmat}{Weight matrix for adaptive Lasso. It is only used when \code{method="aclime"}.}
  \item{data}{\code{n} by \code{d} data matrix. }
  \item{Y}{Dependent variables in linear regression. }
  \item{X}{Design matrix in linear regression. }
  \item{q}{The vector norm used for the loss term. }
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
\code{\link{flare.tiger}}, \code{\link{flare.slim}} and \code{\link{flare-package}}.
}
