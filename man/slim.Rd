\name{slim}
\alias{slim}

\title{
Sparse Linear Regression using Nonsmooth Loss Functions and L1 Regularization
}

\description{
The function "slim" implements a family of Lasso variants for estimating high dimensional sparse linear models including Dantzig Selector, LAD Lasso, SQRT Lasso, Lq Lasso for estimating high dimensional sparse linear model. We adopt the alternating direction method of multipliers (ADMM) and convert the original optimization problem into a sequential L1-penalized least square minimization problem, which can be efficiently solved by combining the linearization and multi-stage screening of varialbes. Missing values can be tolerated for Dantzig selector in the design matrix and response vector. 
}

\usage{
slim(X, Y, lambda = NULL, nlambda = NULL, 
     lambda.min.value = NULL,lambda.min.ratio = NULL, 
     rho = 1, method="lq", q = 2, res.sd = FALSE, 
     prec = 1e-5, max.ite = 1e5, verbose = TRUE)
}

\arguments{
  \item{Y}{The \eqn{n}-dimensional response vector.}
  \item{X}{The \eqn{n} by \eqn{d} design matrix. \code{d}\eqn{ \ge }2 is required.}
  \item{lambda}{A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is from \code{lambda.max} to \code{lambda.min.ratio*lambda.max}. For Lq regression, the default value of \code{lambda.max} is \eqn{\pi\sqrt{\log(d)/n}}. For Dantzig selector, the default value of \code{lambda.max} is the minimum regularization parameter, which yields an all-zero estiamtes.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 5.}
  \item{lambda.min.value}{The smallest value for \code{lambda}, as a fraction of the uppperbound (\code{lambda.max}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{lambda.max} to \code{lambda.min.ratio*lambda.max} in log scale. The default value is \eqn{\log(d)/n} for for Dantzig selector 0.3\code{*lambda.max} for Lq Lasso.}
  \item{lambda.min.ratio}{The smallest ratio of the value for \code{lambda}. The default value is 0.3 for Lq Lasso and 0.5 for Dantzig selector.}
  \item{rho}{The penalty parameter used in \code{ADMM}. The default value is \eqn{\sqrt{d}}.}
  \item{method}{Dantzig selector is applied if \code{method = "dantzig"} and \eqn{L_q} Lasso is applied if \code{method = "lq"}. Standard Lasso is provided if \code{method = "lasso"}. The default value is \code{"lq"}.}
  \item{q}{The loss function used in Lq Lasso. It is only applicable when \code{method = "lq"} and must be in [1,2]. The default value is 2.}
  \item{res.sd}{Flag of whether the response varialbles are standardized. The default value is \code{FALSE}.}
  \item{prec}{Stopping criterion. The default value is 1e-5.}
  \item{max.ite}{The iteration limit. The default value is 1e5.}
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{ 
Standard Lasso\cr
\deqn{
    \min {\frac{1}{2n}}|| Y - X \beta ||_2^2 + \lambda || \beta ||_1
    }\cr
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
An object with S3 class \code{"slim"} is returned: 
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
3. L. Wang. L1 penalized LAD estimator for high dimensional linear regression. \emph{Journal of Multivariate Analysis}, 2012.\cr
4. J. Liu and J. Ye. Efficient L1/Lq Norm Regularization. \emph{Technical Report}, 2010.
5. S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein, Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. \emph{Foundations and Trends in Machine Learning}, 2011.
6. B. He and X. Yuan. On non-ergodic convergence rate of Douglas-Rachford alternating direction method of multipliers. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{flare-package}}, \code{\link{print.slim}}, \code{\link{plot.slim}}, \code{\link{coef.slim}} and \code{\link{predict.slim}}.
}

\examples{
## load library
library(flare)
## generate data
n = 50
d = 100
X = matrix(rnorm(n*d), n, d)
beta = c(3,2,0,1.5,rep(0,d-4))
eps = rnorm(n)
Y = X\%*\%beta + eps
nlamb = 5
ratio = 0.3

## Regression with "dantzig", general "lq" and "lasso" respectively
out1 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="dantzig")
out2 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="lq",q=1)
out3 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="lq",q=1.5)
out4 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="lq",q=2)
out5 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="lasso")

## Display results
print(out4)
plot(out4)
coef(out4)
}
