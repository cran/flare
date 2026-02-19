\name{sugm}
\alias{sugm}

\title{
High-dimensional Sparse Undirected Graphical Models
}

\description{
The function "sugm" estimates sparse undirected graphical models (Gaussian precision matrices) in high dimensions. Two procedures are implemented using a column-wise regression scheme: (1) Tuning-Insensitive Graph Estimation and Regression based on square-root Lasso (\code{"tiger"}); and (2) The Constrained L1 Minimization for Sparse Precision Matrix Estimation (\code{"clime"}). The optimization algorithm is based on the alternating direction method of multipliers (ADMM), linearization, and multi-stage screening. Missing values can be tolerated for CLIME when the input is a data matrix. Computation is memory-optimized using sparse matrix output.
}

\usage{
sugm(data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
     rho = NULL, method = "tiger", sym = "or", shrink = NULL, 
     prec = 1e-4, max.ite = 1e4, standardize = FALSE, 
     perturb = TRUE, verbose = TRUE)
}

\arguments{
  \item{data}{There are two options for \code{"clime"}: (1) an \code{n} by \code{d} data matrix, or (2) a \code{d} by \code{d} sample covariance matrix. The input type is identified by symmetry. For \code{"tiger"}, covariance input is not supported and \code{d}\eqn{ \ge }3 is required. For \code{"clime"}, \code{d}\eqn{ \ge }2 is required.}
  \item{lambda}{A sequence of decreasing, positive, finite numbers controlling regularization. Typical usage is \code{lambda = NULL}, in which case the sequence is generated from \code{nlambda} and \code{lambda.min.ratio}.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.ratio}{The minimum value of generated \code{lambda} as a fraction of \code{lambda.max}. The default value is 0.4 for both \code{"tiger"} and \code{"clime"}.}
  \item{rho}{Penalty parameter used in the optimization algorithm. The default value is 1.}
  \item{method}{\code{"tiger"} is applied if \code{method = "tiger"} and \code{"clime"} is applied if \code{method="clime"}. Default value is \code{"tiger"}.}
  \item{sym}{Symmetrization of output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}.}
  \item{shrink}{Shrinkage of the regularization parameter based on estimation precision. The default value is 0.}
  \item{prec}{Stopping criterion. The default value is 1e-4.}
  \item{max.ite}{The iteration limit. The default value is 1e4.}
  \item{standardize}{Variables are standardized to have mean zero and unit standard deviation if \code{standardize = TRUE}. The default value is \code{FALSE}.}
  \item{perturb}{For \code{"clime"}, if \code{TRUE}, adds \eqn{1/\sqrt{n}} to the diagonal of \code{Sigma}; if \code{FALSE}, no perturbation is added; a numeric value can also be supplied directly. The default value is \code{TRUE}.}
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
  Iteration counts returned by the underlying optimization solver.
  }
  \item{df}{
  A \code{d} by \code{nlambda} matrix containing nonzero counts along the estimated path.
  }
  \item{standardize}{The \code{standardize} from the input.}
  \item{perturb}{The \code{perturb} from the input.}
  \item{verbose}{The \code{verbose} from the input.}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Tuo Zhao <tourzhao@gatech.edu>
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
