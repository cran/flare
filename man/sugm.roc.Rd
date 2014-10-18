\name{sugm.roc}
\alias{sugm.roc}

\title{
Draw ROC Curve for an object with S3 class \code{"sugm"}
}

\description{
Draws ROC curve for a graph path according to the true graph structure.
}

\usage{
sugm.roc(path, theta, verbose = TRUE)
}

\arguments{
  \item{path}{
A graph path.
}
  \item{theta}{
The true graph structure.
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
To avoid the horizontal oscillation, false positive rates is automatically sorted in the ascent oder and true positive rates also follow the same order.
}

\note{
For a lasso regression,  the number of nonzero coefficients is at most \code{n-1}. If \code{d>>n}, even when regularization parameter is very small, the estimated graph may still be sparse. In this case, the AUC may not be a good choice to evaluate the performance.
}

\value{
An object with S3 class "roc" is returned:
 \item{F1}{
The F1 scores along the graph path.
}
  \item{tp}{
The true positive rates along the graph path
}
  \item{fp}{
The false positive rates along the graph paths
}
  \item{AUC}{
Area under the ROC curve
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\seealso{
\code{\link{sugm}} and \code{\link{flare-package}}
}

\examples{
## load package required
library(flare)

#generate data
L = sugm.generator(d = 30, graph = "random", prob = 0.1)
out1 = sugm(L$data, lambda=10^(seq(log10(.4), log10(0.03), length.out=20)))

#draw ROC curve
Z1 = sugm.roc(out1$path,L$theta)

#Maximum F1 score
max(Z1$F1)
}
