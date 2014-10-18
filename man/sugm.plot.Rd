\name{sugm.plot}
\alias{sugm.plot}

\title{
Graph visualization for an object with S3 class \code{"sugm"}
}

\description{
Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
}

\usage{
sugm.plot(G, epsflag = FALSE, graph.name = "default", cur.num = 1, 
          location)
}

\arguments{
  \item{G}{
The adjacency matrix corresponding to the graph.
}
  \item{epsflag}{
If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.
}
  \item{graph.name}{
The name of the output eps files. The default value is "default".
}
  \item{cur.num}{
The number of plots saved as eps files. Only applicable when \code{epsflag = TRUE}. The default value is 1.
}
  \item{location}{
Target directory. The default value is the current working directory.
}
}

\details{
The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\seealso{
\code{\link{flare}} and \code{\link{flare-package}}
}

\examples{
## load package required
library(flare)

## visualize the hub graph
L = sugm.generator(graph = "hub")
sugm.plot(L$theta)

## visualize the band graph
L = sugm.generator(graph = "band",g=5)
sugm.plot(L$theta)

## visualize the cluster graph
L = sugm.generator(graph = "cluster")
sugm.plot(L$theta)

#show working directory
getwd()
#plot 5 graphs and save the plots as eps files in the working directory  
sugm.plot(L$theta, epsflag = TRUE, cur.num = 5)
}
