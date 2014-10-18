\name{flare-package}
\alias{flare-package}
\alias{flare}
\docType{package}
\title{
flare: a new Family of Lasso Regression
}
\description{
The package "flare" provides the implementation of a family of novel regression methods (Lasso, Dantzig Selector, LAD Lasso, SQRT Lasso, Lq Lasso) and their extensions to sparse precision matrix estimation (TIGER and CLIME using L1) in high dimensions. We adopt the alternating direction method of multipliers and convert the original optimization problem into a sequence of L1-penalized least square minimization problems with the linearization method and multi-stage screening of variables. Missing values can be tolerated for Dantzig selector in the design matrix and response vector, and CLIME in the data matrix. The computation is memory-optimized using the sparse matrix output. In addition, we also provide several convenient regularization parameter selection and visulaization tools.
}
\details{
\tabular{ll}{
Package: \tab flare\cr
Type: \tab Package\cr
Version: \tab 1.5.0\cr
Date: \tab 2014-10-18\cr
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
3. L. Wang. L1 penalized LAD estimator for high dimensional linear regression. \emph{Journal of Multivariate Analysis}, 2012.\cr
4. J. Liu and J. Ye. Efficient L1/Lq Norm Regularization. \emph{Technical Report}, 2010.
5. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011.\cr
6. S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein, Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. \emph{Foundations and Trends in Machine Learning}, 2011.
7. H. Liu and L. Wang. TIGER: A tuning-insensitive approach for optimally estimating large undirected graphs. \emph{Technical Report}, 2012. \cr
8. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{sugm}} and \code{\link{slim}}.
}
