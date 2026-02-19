\name{sugm.select}
\alias{sugm.select}

\title{
Model selection for high-dimensional undirected graphical models
}

\description{
Implements regularization parameter selection for high-dimensional undirected graphical models. Supported approaches are Stability Approach to Regularization Selection (\code{"stars"}) and cross-validation (\code{"cv"}).
}

\usage{
sugm.select(est, criterion = "stars", stars.subsample.ratio = NULL,
            stars.thresh = 0.1, rep.num = 20, fold = 5,
            loss="likelihood", verbose = TRUE)
}

\arguments{
  \item{est}{
An object with S3 class \code{"sugm"}
}
  \item{criterion}{
Model selection criterion. \code{"stars"} and \code{"cv"} are available for both graph estimation methods. The default value is \code{"stars"}.
}
  \item{stars.subsample.ratio}{
The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n > 144} and \code{0.8} when \code{n <= 144}, where \code{n} is the sample size. Must be in \code{(0,1)}. Only applicable when \code{criterion = "stars"}.
}
  \item{stars.thresh}{
The variability threshold in STARS. Must be in \code{[0,1]}. The default value is \code{0.1}. Only applicable when \code{criterion = "stars"}.
}
  \item{rep.num}{
The number of subsamples. Must be at least 1. The default value is \code{20}.
}
  \item{fold}{
The number of folds used in cross-validation. Must be between 2 and \code{n}. The default value is \code{5}. Only applicable when \code{criterion = "cv"}.
}
  \item{loss}{
  Loss used in cross-validation. Two losses are available: \code{"likelihood"} and \code{"tracel2"}. Default is \code{"likelihood"}. Only applicable when \code{criterion = "cv"}.
  }
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
Stability Approach to Regularization Selection (STARS) selects an optimal regularization parameter by variability across subsamples, and tends to over-select edges in Gaussian graphical models. In addition to selecting regularization parameters, STARS can provide an additional merged graph estimate based on edge frequencies across subsamples. K-fold cross-validation is also available for selecting \code{lambda}, using the following losses:
\deqn{
likelihood: Tr(\Sigma \Omega) - \log|\Omega|
}
\deqn{
tracel2: Tr(diag(\Sigma \Omega - I)^2).
}
}

\value{
An object with S3 class "select" is returned:
  \item{refit}{
The optimal graph selected from the graph path
}
  \item{opt.icov}{
The optimal precision matrix selected.
}
  \item{merge}{
The graph path estimated by merging the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
}
  \item{variability}{
The variability along the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
}
  \item{opt.index}{
The index of the selected regularization parameter.
}
  \item{opt.lambda}{
The selected regularization/thresholding parameter.
}
  \item{opt.sparsity}{
The sparsity level of \code{"refit"}.
}
\item{loss}{
Cross-validation loss used for selection. Only applicable when \code{criterion = "cv"}.
}

and anything else included in the input \code{est}.
}

\note{
Model selection is not available when the data input is a sample covariance matrix.
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Tuo Zhao <tourzhao@gatech.edu>
}

\references{
1. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011.\cr
2. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{sugm}} and \code{\link{flare-package}}.
}

\examples{
## load package required
library(flare)

#generate data
L = sugm.generator(d = 10, graph="hub")
out1 = sugm(L$data)

#model selection using stars
#out1.select1 = sugm.select(out1, criterion = "stars", stars.thresh = 0.1)
#plot(out1.select1)

#model selection using cross validation
out1.select2 = sugm.select(out1, criterion = "cv")
plot(out1.select2)
}
