\name{coef.slim}
\alias{coef.slim}
\title{
Extract Model Coefficients for an object with S3 class \code{"slim"}
}
\description{
Extract estimated regression coefficient vectors from the solution path.
}
\usage{
\method{coef}{slim}(object, lambda.idx = c(1:3), beta.idx = c(1:3), \dots)
}
\arguments{
  \item{object}{
  An object with S3 class \code{"slim"} 
  }
  \item{lambda.idx}{
  The indices of the regularization parameters in the solution path to display. The default values are \code{c(1:3)}.
  }
  \item{beta.idx}{
  The indices of estimated regression coefficient vectors in the solution path to display. The default values are \code{c(1:3)}.
  }
  \item{\dots}{
  Arguments to be passed to methods.
  }
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Tuo Zhao <tourzhao@gatech.edu>
}

\seealso{
\code{\link{slim}} and \code{\link{flare-package}}.
}
