\name{flare.slim}
\alias{flare.slim}

\title{
Sparse Linear Regression using Non-smooth Loss Functions and L1 Regularization
}

\description{
The function "flare.slim" implements a family of Lasso variants for estimating high dimensional sparse linear models including Dantzig Selector, LAD Lasso, SQRT Lasso, Lq Lasso for estimating high dimensional sparse linear model. We adopt the alternating direction method of multipliers (ADMM) and convert the original optimization problem into a sequential L1-penalized least square minimization problem, which can be efficiently solved by combining the linearization and the efficient coordinate descent algorithm. The computation is memory-optimized using the sparse matrix output. 
}

\usage{
flare.slim(X, Y,  lambda = NULL, nlambda = NULL,
       lambda.min.ratio = NULL, rho = NULL, method="lq", 
       q = 2, prec = 1e-3, max.ite = 1e3, verbose = TRUE)
}

\arguments{
  \item{Y}{The \eqn{n} dimensional response vector.}
  \item{X}{The \eqn{n} by \eqn{d} design matrix.}
  \item{lambda}{A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is from \eqn{lambda.max} to \code{lambda.min.ratio*lambda.max}. For Lq regression, the default value of \eqn{lambda.max} is \eqn{\pi\sqrt{\log(d)/n}}. For Dantzig selector, the default value of \eqn{lambda.max} is the minimum regularization parameter, which yields an all-zero estiamtes.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.25} for Lq Lasso and \code{0.5} for Dantzig selector.}
  \item{rho}{The penalty parameter used in \code{ADMM}. The default value is \eqn{\sqrt{d}}.}
  \item{method}{Dantzig selector is applied if \code{method = "dantzig"} and \eqn{L_q} Lasso is applied if \code{method = "lq"}. The default value is \code{"lq"}.}
  \item{q}{The loss function used in Lq Lasso. It is only applicable when \code{method = "lq"} and must be in \eqn{[1,2]}. The default value is 2.}
  \item{prec}{Stopping criterion. The default value is 1e-3.}
  \item{max.ite}{The iteration limit. The default value is 1e3.}
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{ 
Dantzig selector solves the following optimization problem\cr
\deqn{
    \min || \beta ||_1, \quad \textrm{s.t. } || X'(Y - X \beta) ||_{\infty} < \lambda
    }\cr
\eqn{L_q} loss Lasso solves the following optimization problem\cr
\deqn{
    \min n^{-\frac{1}{q}}|| Y - X \beta ||_q + \lambda || \beta ||_1
    }\cr
    where \eqn{1<= q <=2}. Lq Lasso is equivalent to LAD Lasso and SQR Lasso when \eqn{q=1} and \eqn{q=2} respectively.
}
\value{
An object with S3 class \code{"flare.slim"} is returned: 
  \item{beta}{
  A matrix of regression estimates whose columns correspond to regularization parameters.
  }
  \item{intercept}{The value of intercepts corresponding to regularization parameters. }
  \item{Y}{The value of \code{Y} used in the program. }
  \item{X}{The value of \code{X} used in the program. }
  \item{lambda}{
  The sequence of regularization parameters \code{lambda} used in the program.
  }
  \item{nlambda}{
  The number of values used in \code{lambda}.
  }
  \item{method}{
  The \code{method} from the input. 
  }
  \item{sparsity}{
  The sparsity levels of the solution path.
  }
  \item{ite}{
  A list of vectors where ite[[1]] is the number of external iteration and ite[[2]] is the number of internal iteration with the i-th entry corresponding to the i-th regularization parameter. 
  }
  \item{verbose}{The \code{verbose} from the input.}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. E. Candes and T. Tao. The Dantzig selector: Statistical estimation when p is much larger than n. \emph{Annals of Statistics}, 2007.\cr
2. A. Belloni, V. Chernozhukov and L. Wang. Pivotal recovery of sparse signals via conic programming. \emph{Biometrika}, 2012.\cr
3. L. Wang. L1 penalized LAD estimator for high dimensional linear regression. \emph{Technical Report}, 2012.\cr
4. J. Liu and J. Ye. Efficient L1/Lq Norm Regularization. \emph{Technical Report}, 2010.\cr
5. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{flare-package}}.
}

\examples{
## generate data
set.seed(1)
n=100
d=200
d1=10
rho0=0.3
lambda=c(3:1)*sqrt(log(d)/n)
Sigma=matrix(0,nrow=d,ncol=d)
Sigma[1:d1,1:d1]=rho0
diag(Sigma)=1
mu=rep(0,d)
X=mvrnorm(n=n,mu=mu,Sigma=Sigma)
eps=rt(n=n,df=n-1)
beta=c(rep(sqrt(1/3),3),rep(0,d-3))
Y=X\%*\%beta+eps

## Regression with "dantzig" and general "lq" respectively
out1=flare.slim(X=X,Y=Y,lambda=lambda,method = "dantzig")
out2=flare.slim(X=X,Y=Y,lambda=lambda,method = "lq",q=1.5)

## Print results
print(out1)
print(out2)
plot(out1)
}
