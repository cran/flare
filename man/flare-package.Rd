\name{flare-package}
\alias{flare-package}
\alias{flare}
\docType{package}
\title{
flare: Family of Lasso Regression
}
\description{
The package "flare" provides the implementation of a family of new regression methods (Dantzig Selector, LAD Lasso, SQRT Lasso, Lq Lasso) and their extensions to sparse precision matrix estimation (TIGER and CLIME using either L1 or adaptive penalty) in high dimensions. We adopt the alternating direction method of multipliers and convert the original optimization problem into a sequential L1 penalized least square minimization problem, which can be efficiently solved by combining the linearization and efficient coordinate descent algorithm. The computation is memory-optimized using the sparse matrix output.
}
\details{
\tabular{ll}{
Package: \tab flare\cr
Type: \tab Package\cr
Version: \tab 0.9.6\cr
Date: \tab 2012-11-25\cr
License: \tab GPL-2\cr
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang , Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. E. Candes and T. Tao. The Dantzig selector: Statistical estimation when p is much larger than n. \emph{Annals of Statistics}, 2007.\cr
2. A. Belloni, V. Chernozhukov and L. Wang. Pivotal recovery of sparse signals via conic programming. \emph{Biometrika}, 2012.\cr
3. L. Wang. L1 penalized LAD estimator for high dimensional linear regression. \emph{Technical Report}, 2012.\cr
4. J. Liu and J. Ye. Efficient L1/Lq Norm Regularization. \emph{Technical Report}, 2010.
5. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011.\cr
6. H. Liu and L. Wang. TIGER: A tuning-insensitive approach for optimally estimating large undirected graphs. \emph{Technical Report}, 2012. \cr
7. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{flare.tiger}} and \code{\link{flare.slim}}.
}
