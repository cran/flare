\name{sugm}
\alias{sugm}

\title{
High-deimensional Sparse Undirected Graphical Models.
}

\description{
The function "sugm" estimates sparse undirected graphical models, i.e. Gaussian precision matrix, in high dimensions. We adopt two estimation procedures based on column by column regression scheme: (1) Tuning-Insensitive Graph Estimation and Regression based on square root Lasso (tiger); (2) The Constrained L1 Minimization for Sparse Precision Matrix Estimation using either L1 penalty (clime). The optimization algorithm for all three methods are implemented based on the alternating direction method of multipliers (ADMM) with the linearization method and multi-stage screening of variables. Missing values can be tolerated for CLIME in the data matrix. The computation is memory-optimized using the sparse matrix output.
}

\usage{
sugm(data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
     rho = NULL, method = "tiger", sym = "or", shrink=NULL, 
     prec = 1e-4, max.ite = 1e4, standardize = FALSE, 
     perturb = TRUE, verbose = TRUE)
}

\arguments{
  \item{data}{There are 2 options for \code{"clime"}: (1) \code{data} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension). For \code{"tiger"}, covariance input is not supported and \code{d}\eqn{ \ge }3 is required. For \code{"clime"}, \code{d}\eqn{ \ge }2 is required.}
  \item{lambda}{A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is from \code{lambda.max} to \code{lambda.min.ratio*lambda.max}. For \code{"tiger"}, the default value of \code{lambda.max} is \eqn{\pi\sqrt{\log(d)/n}}. For \code{"clime"} , the default value of \code{lambda.max} is the minimum regularization parameter, which yields an all-zero off-diagonal estiamtes.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{lambda.max}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{lambda.max} to \code{lambda.min.ratio*lambda.max} in log scale. The default value is 0.25 for \code{"tiger"} and 0.5 for \code{"clime"}.}
  \item{rho}{Penalty parameter used in the optimization algorithm for \code{clime}. The default value is \eqn{\sqrt{d}}.}
  \item{method}{\code{"tiger"} is applied if \code{method = "tiger"} and \code{"clime"} is applied if \code{method="clime"}. Default value is \code{"tiger"}.}
  \item{sym}{Symmetrization of output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}.}
  \item{shrink}{Shrinkage of regularization parameter based on precision of estimation. The default value is 1.5 if \code{method = "clime"} and the default value is 0 if \code{method="tiger"}.}
  \item{prec}{Stopping criterion. The default value is 1e-4.}
  \item{max.ite}{The iteration limit. The default value is 1e4.}
  \item{standardize}{Variables are standardized to have mean zero and unit standard deviation if \code{standardize = TRUE}. The default value is \code{FALSE}.}
  \item{perturb}{The diagonal of \code{Sigma} is added by a positive value to guarantee that \code{Sigma} is positive definite if \code{perturb = TRUE}. User can specify a numeric value for \code{perturbe}. The default value is \code{perturb = TRUE}.}
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{
CLIME solves the following minimization problem\cr
\deqn{
    \min || \Omega ||_1 \quad \textrm{s.t. } || S \Omega - I ||_\infty \le \lambda,
  }\cr
where \eqn{||\cdot||_1} and \eqn{||\cdot||_\infty} are element-wise 1-norm and \eqn{\infty}-norm respectively. 

\code{"tiger"} solves the following minimization problem\cr
\deqn{
    \min ||X-XB||_{2,1} + \lambda ||B||_1 \quad \textrm{s.t. } B_{jj} = 0,
  }\cr
where \eqn{||\cdot||_{1}} and \eqn{||\cdot||_{2,1}} are element-wise 1-norm and \eqn{L_{2,1}}-norm respectively.
}

\value{
An object with S3 class \code{"sugm"} is returned:  
  \item{data}{
  The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input.
  }
  \item{cov.input}{
  An indicator of the sample covariance. 
  }
  \item{lambda}{
  The sequence of regularization parameters \code{lambda} used in the program.
  }
  \item{nlambda}{
  The number of values used in \code{lambda}.
  }
  \item{icov}{
  A list of \code{d} by \code{d} precision matrices corresponding to regularization parameters.
  }
  \item{sym}{
  The \code{sym} from the input.
  }
  \item{method}{
  The \code{method} from the input. 
  }
  \item{path}{
  A list of \code{d} by \code{d} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
  }
  \item{sparsity}{
  The sparsity levels of the graph path.
  }
  \item{ite}{
  If \code{method = "clime"}, it is a list of two matrices where ite[[1]] is the number of external iterations and ite[[2]] is the number of internal iterations with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. If \code{method="tiger"}, it is a matrix of iteration with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. 
  }
  \item{df}{
  It is a \code{d} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the lasso solution path. 
  }
  \item{standardize}{The \code{standardize} from the input.}
  \item{perturb}{The \code{perturb} from the input.}
  \item{verbose}{The \code{verbose} from the input.}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. T. Cai, W. Liu and X. Luo. A constrained L1 minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011. \cr
2. H. Liu, L. Wang. TIGER: A tuning-insensitive approach for optimally estimating large undirected graphs. \emph{Technical Report}, 2012. \cr
3. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{flare-package}}, \code{\link{sugm.generator}}, \code{\link{sugm.select}}, \code{\link{sugm.plot}}, \code{\link{sugm.roc}}, \code{\link{plot.sugm}}, \code{\link{plot.select}}, \code{\link{plot.roc}}, \code{\link{plot.sim}}, \code{\link{print.sugm}}, \code{\link{print.select}}, \code{\link{print.roc}} and \code{\link{print.sim}}.
}

\examples{

## load package required
library(flare)

## generating data
n = 50
d = 50
D = sugm.generator(n=n,d=d,graph="band",g=1)
plot(D)

## sparse precision matrix estimation with method "clime"
out1 = sugm(D$data, method = "clime")
plot(out1)
sugm.plot(out1$path[[4]])

## sparse precision matrix estimation with method "tiger"
out2 = sugm(D$data, method = "tiger")
plot(out2)
sugm.plot(out2$path[[5]])
}
