\name{flare-package}
\alias{flare-package}
\alias{flare}
\docType{package}
\title{
flare: Family of Lasso Regression
}
\description{
The package "flare" provides the implementation of a family of high-dimensional Lasso regression based machine learning toolkits, including a family of various Lasso regression and sparse Gaussian graphical model estimation. Lasso variants including Dantzig Selector, LAD Lasso, SQRT Lasso for estimating high dimensional sparse linear model. The sparse Gaussian graphical model estimation includes TIGER and CLIME using L1 penalty. We adopt the combination of the dual smoothing and monotone fast iterative soft-thresholding algorithm (MFISTA). The computation is memory-optimized using the sparse matrix output. 
}
\details{
\tabular{ll}{
Package: \tab flare\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2013-08-25\cr
License: \tab GPL-2\cr
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. A. Belloni, V. Chernozhukov and L. Wang. Pivotal recovery of sparse signals via conic programming. \emph{Biometrika}, 2012.\cr
2. L. Wang. L1 penalized LAD estimator for high dimensional linear regression. \emph{Journal of Multivariate Analysis}, 2013.\cr
3. E. Candes and T. Tao. The Dantzig selector: Statistical estimation when p is much larger than n. \emph{Annals of Statistics}, 2007.\cr
4. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011.\cr
5. H. Liu and L. Wang. TIGER: A tuning-insensitive approach for optimally estimating large undirected graphs. \emph{Technical Report}, 2012. \cr
6. A. Beck and M. Teboulle. Fast gradient-based algorithms for constrained total variation image denoising and deblurring problems. \emph{IEEE Transactions on Image Processing}, 2009.\cr
7. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012. \cr
8. J. Liu and J. Ye. Efficient L1/Lq Norm Regularization. \emph{Technical Report}, 2010.
}

\seealso{
\code{\link{flare.tiger}} and \code{\link{flare.slim}}.
}
